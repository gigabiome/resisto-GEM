#!/usr/bin/env python
"""Prepare a curated reference GEM for reference-guided CarveMe reconstruction.

Runs once per reference model. Outputs (prefix = --reference-id):
  <prefix>.harmonised.xml      reference with compartments renamed to CarveMe's C_<x> convention
                               (passed to `carve --reference`)
  <prefix>.biomass_module.xml  the selected biomass reaction plus the reference reactions that
                               build its lumped precursors or drain its non-currency products;
                               added unchanged to every strain model
  <prefix>.pool.xml            gap-fill pool: CarveMe universe (minus its generic biomass) plus
                               reference reactions that are not in the universe
  <prefix>.pool_sources.tsv    origin of every pool reaction and any ID conflicts
  <prefix>.biomass.json        biomass ID, module reactions and ATPM bounds
  <prefix>.validation.json     all checks below; the run fails if gap-filling would be infeasible

Rules (all recorded in the outputs):
  * IDs are not changed, except compartments (x -> C_x) and the reference variants described next.
  * Reaction ID in both universe and reference with a different stoichiometry:
      - the universe version stays in the pool under the original ID (CarveMe drafts use it);
      - if the reference version uses metabolites that are not in the universe (organism-specific
        chemistry, e.g. P. aeruginosa lipid X in R_LPADSS), it is also added as R_<id>__<ref>;
      - otherwise (same metabolites, different coefficients) it is not added, to avoid near-duplicate
        reactions that could form energy-generating cycles.
  * Biomass module: starting from the selected biomass reaction, repeatedly add "biomass-like"
    reference reactions (objective, or 'BIOMASS' in the ID, or non-exchange boundary reactions) that
    produce a non-currency reactant, or consume a non-currency product, of a module reaction.
    Other biomass-like reactions (alternative biomass equations) are left out.
  * Not copied into the pool: module reactions (added directly), other biomass-like reactions,
    R_ATPM (the universe has it) and non-exchange boundary reactions that can supply a metabolite
    ("sinks", lb < 0, including mislabelled exchanges such as iDR478's R_Ex_lac_L__e), because
    CarveMe's medium handling only controls R_EX_ reactions.
  * Required sinks: if pool + module cannot grow, the smallest set of reference sinks that restores
    growth is found greedily (start with all, drop each in alphabetical order if growth survives).
    These sinks join the biomass module, i.e. are added to every strain model; e.g. iDR478 needs
    SK_apoACP_c because its lipoteichoic acid / cardiolipin pathways consume ACP.
  * Checks: (1) harmonised reference grows with its own biomass on the medium;
            (2) pool + biomass module grows on the medium. A CarveMe draft is a subset of the
                universe, so this guarantees gap-filling is feasible for every genome.
    If (2) fails, per-precursor producibility is reported as a diagnostic.
"""
import argparse
import csv
import math
import re
import sys

from reframed import FBA, Environment, load_cbmodel
from reframed.core.cbmodel import CBModel, CBReaction

from resisto_common import (CARVEME_BIOMASS, copy_reaction_into, harmonise_compartments,
                            load_medium, load_universe, read_fbc_fields, rid, save_model,
                            set_solver, software_versions, universe_path, write_json)

CURRENCY = {'atp', 'adp', 'amp', 'gtp', 'gdp', 'pi', 'ppi', 'h', 'h2o'}
ELEMENT = re.compile(r'([A-Z][a-z]*)(\d*\.?\d*)')


def elements(formula):
    """Element counts of a formula string; None if missing or not parseable (e.g. 'R' groups)."""
    if not formula or not isinstance(formula, str):
        return None
    counts = {}
    for el, n in ELEMENT.findall(formula):
        counts[el] = counts.get(el, 0) + (float(n) if n else 1.0)
    return counts


def use_fbc_fields(model, sbml_path):
    """Make each reference metabolite's FORMULA/CHARGE equal to its SBML fbc:chemicalFormula /
    fbc:charge where those are set. reframed otherwise takes them from the notes, which can hold
    several ';'-joined formulas (iSD1509: 158 metabolites, silently dropped on writing) or
    charges that disagree with fbc. Returns counts of changed formulas and charges."""
    fields = read_fbc_fields(sbml_path)
    changed = {'formula': 0, 'charge': 0, 'multi_formula_unresolved': 0}
    for m_id, met in model.metabolites.items():
        formula, charge = fields.get(m_id, (None, None))
        old_f, old_c = met.metadata.get('FORMULA'), met.metadata.get('CHARGE')
        if formula and formula != old_f:
            met.metadata['FORMULA'] = formula
            changed['formula'] += 1
        elif isinstance(old_f, str) and ';' in old_f:
            met.metadata.pop('FORMULA')
            changed['multi_formula_unresolved'] += 1
        if charge is not None and str(charge) != str(old_c):
            met.metadata['CHARGE'] = str(charge)
            changed['charge'] += 1
    return changed


def balance(model, rxn):
    """'balanced', 'imbalanced' or 'unknown' (missing/unparseable formula or charge)."""
    total, charge = {}, 0.0
    for m_id, coeff in rxn.stoichiometry.items():
        meta = model.metabolites[m_id].metadata
        el = elements(meta.get('FORMULA'))
        ch = meta.get('CHARGE')
        if el is None or ch in (None, ''):
            return 'unknown'
        for k, v in el.items():
            total[k] = total.get(k, 0.0) + coeff * v
        charge += coeff * float(ch)
    ok = all(abs(v) < 1e-6 for v in total.values()) and abs(charge) < 1e-6
    return 'balanced' if ok else 'imbalanced'


def base_id(m_id):
    """M_atp_c -> atp (strip prefix and compartment suffix)."""
    core = m_id[2:] if m_id.startswith('M_') else m_id
    return core.rsplit('_', 1)[0] if '_' in core else core


def file_sha256(path):
    import hashlib
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def finite(x):
    return None if x is None or math.isinf(x) else x


def fba(model, objective, constraints):
    sol = FBA(model, objective=objective, constraints=constraints, get_values=False)
    ok = sol.status.name in ('OPTIMAL', 'SUBOPTIMAL')
    return (sol.fobj if ok and sol.fobj is not None else 0.0), sol.status.name


def biomass_module(ref, biomass_id):
    def biomass_like(r_id, rxn):
        boundary = len(rxn.stoichiometry) == 1 and not r_id.startswith('R_EX_')
        return rxn.objective != 0 or 'BIOMASS' in r_id.upper() or boundary

    candidates = {r for r, x in ref.reactions.items() if r != biomass_id and biomass_like(r, x)}
    module = [biomass_id]
    changed = True
    while changed:
        changed = False
        needed, drained = set(), set()
        for r in module:
            for m, c in ref.reactions[r].stoichiometry.items():
                if base_id(m) in CURRENCY:
                    continue
                (needed if c < 0 else drained).add(m)
        for r in sorted(candidates - set(module)):
            st = ref.reactions[r].stoichiometry
            if any(st.get(m, 0) > 0 for m in needed) or any(st.get(m, 0) < 0 for m in drained):
                module.append(r)
                changed = True
    excluded = sorted(r for r in candidates - set(module)
                      if ref.reactions[r].objective != 0 or 'BIOMASS' in r.upper())
    return module, excluded


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--reference', required=True, help='Reference GEM (SBML)')
    p.add_argument('--reference-id', required=True)
    p.add_argument('--sha256', default=None, help='Expected SHA-256 of the reference file (from the registry)')
    p.add_argument('--biomass', required=True, help='Biomass reaction ID in the reference (with or without R_)')
    p.add_argument('--universe', required=True, choices=['gramneg', 'grampos', 'bacteria'])
    p.add_argument('--medium', required=True, help='Medium ID in the media database')
    p.add_argument('--mediadb', default=None, help='Media database (CarveMe format); default: CarveMe built-in')
    p.add_argument('--max-uptake', type=float, default=100.0,
                   help='Uptake bound for medium compounds in the checks (CarveMe gap-fill uses 100)')
    p.add_argument('--atpm-source', choices=['reference', 'carveme'], default='reference')
    p.add_argument('--min-growth', type=float, default=0.1)
    p.add_argument('--solver', default='gurobi')
    p.add_argument('--threads', type=int, default=1)
    args = p.parse_args()

    sha = file_sha256(args.reference)
    if args.sha256 and sha != args.sha256.lower():
        sys.exit(f'{args.reference}: SHA-256 {sha} does not match the registry ({args.sha256}). '
                 'Re-download the reference or update assets/references.csv deliberately.')
    set_solver(args.solver, args.threads)
    prefix = args.reference_id
    biomass_id = rid(args.biomass)

    universe = load_universe(args.universe)
    if universe.biomass_reaction != CARVEME_BIOMASS:
        sys.exit(f'Universe biomass is {universe.biomass_reaction!r}, expected {CARVEME_BIOMASS!r}')
    if 'R_ATPM' not in universe.reactions:
        sys.exit('Universe has no R_ATPM reaction')

    ref = load_cbmodel(args.reference, flavor='bigg')
    fbc_changes = use_fbc_fields(ref, args.reference)
    if biomass_id not in ref.reactions:
        sys.exit(f'Biomass reaction {biomass_id} not found in {args.reference}')
    comp_map = harmonise_compartments(ref, universe)
    ref.id = f'{prefix}_harmonised'
    for r_id, rxn in ref.reactions.items():
        rxn.objective = 1 if r_id == biomass_id else 0
    ref.biomass_reaction = biomass_id

    module, excluded_biomass = biomass_module(ref, biomass_id)

    # --- biomass module model ------------------------------------------------------------------
    mod = CBModel(f'{prefix}_biomass_module')
    for r_id in module:
        copy_reaction_into(mod, ref, r_id)
    for r_id, rxn in mod.reactions.items():
        rxn.objective = 1 if r_id == biomass_id else 0
        rxn.gpr = None

    # --- gap-fill pool -------------------------------------------------------------------------
    pool = universe.copy()
    pool.id = f'{prefix}_pool'
    pool.remove_reactions([CARVEME_BIOMASS])
    universe_mets = set(universe.metabolites)

    rows, skipped = [], []
    for r_id, rxn in ref.reactions.items():
        if r_id in module:
            continue
        if r_id in excluded_biomass:
            skipped.append((r_id, 'alternative biomass reaction'))
            continue
        if r_id == 'R_ATPM':
            continue
        if len(rxn.stoichiometry) == 1 and not r_id.startswith('R_EX_') and (rxn.lb is None or rxn.lb < 0):
            skipped.append((r_id, 'sink (can supply a metabolite)'))
            continue
        if r_id in universe.reactions:
            su, sr = universe.reactions[r_id].stoichiometry, rxn.stoichiometry
            same = set(su) == set(sr) and all(abs(su[m] - sr[m]) <= 1e-9 for m in su)
            if same:
                rows.append({'reaction': r_id, 'reference_reaction': r_id, 'in_universe': 1,
                             'in_reference': 1, 'note': ''})
                continue
            new_mets = sorted(set(sr) - universe_mets)
            if new_mets:
                variant = f'{r_id}__{prefix}'
                copy_reaction_into(pool, ref, r_id, new_id=variant)
                rows.append({'reaction': r_id, 'reference_reaction': variant, 'in_universe': 1,
                             'in_reference': 1,
                             'note': f'stoichiometry differs; reference variant added as {variant} '
                                     f'(reference-specific metabolites: {",".join(new_mets)})'})
                rows.append({'reaction': variant, 'reference_reaction': r_id, 'in_universe': 0,
                             'in_reference': 1, 'note': 'reference variant of a universe reaction ID'})
            else:
                rows.append({'reaction': r_id, 'reference_reaction': r_id, 'in_universe': 1,
                             'in_reference': 1,
                             'note': 'stoichiometry differs (coefficients only); universe version used'})
            continue
        copy_reaction_into(pool, ref, r_id)
        rows.append({'reaction': r_id, 'reference_reaction': r_id, 'in_universe': 0,
                     'in_reference': 1, 'note': ''})
    listed = {r['reaction'] for r in rows}
    for r_id in pool.reactions:
        if r_id not in listed:
            rows.append({'reaction': r_id, 'reference_reaction': '', 'in_universe': 1,
                         'in_reference': 0, 'note': ''})
    # Pool reactions carry no gene rules: reference GPRs name the reference strain's genes, which do
    # not exist in a new strain (gap-filled reactions get no GPR, as in CarveMe).
    for rxn in pool.reactions.values():
        rxn.objective = 0
        rxn.gpr = None
    pool.genes.clear()

    # Mass/charge balance of reference reactions, in the reference and in the pool. Shared BiGG
    # metabolite IDs can carry different formulas/charges (protonation states) in the universe,
    # so a reaction balanced in the reference may be imbalanced on universe metabolites.
    for row in rows:
        row['balance_in_reference'] = balance(ref, ref.reactions[row['reference_reaction']]) \
            if row['in_reference'] and row['reference_reaction'] in ref.reactions else ''
        row['balance_in_pool'] = balance(pool, pool.reactions[row['reaction']])
    degraded = sorted(r['reaction'] for r in rows if r['in_reference'] and not r['in_universe']
                      and r['balance_in_reference'] == 'balanced' and r['balance_in_pool'] == 'imbalanced')

    # --- ATPM ----------------------------------------------------------------------------------
    ref_atpm = ref.reactions.get('R_ATPM')
    atpm = {'reference': (finite(ref_atpm.lb), finite(ref_atpm.ub)) if ref_atpm else None,
            'carveme': (finite(universe.reactions['R_ATPM'].lb), finite(universe.reactions['R_ATPM'].ub))}
    if args.atpm_source == 'reference' and atpm['reference'] is None:
        sys.exit('Reference has no R_ATPM; use --atpm-source carveme')
    atpm_bounds = atpm[args.atpm_source]

    # --- checks --------------------------------------------------------------------------------
    compounds, mediadb = load_medium(args.medium, args.mediadb)
    env = Environment.from_compounds(compounds, max_uptake=args.max_uptake)

    ref_check = ref.copy()
    if args.atpm_source == 'carveme' and ref_atpm is not None:
        ref_check.set_flux_bounds('R_ATPM', *atpm['carveme'])
    ref_growth, ref_status = fba(ref_check, {biomass_id: 1},
                                 env.apply(ref_check, inplace=False, warning=False))

    def pool_with(extra_sinks):
        model = pool.copy()
        for r_id in module:
            copy_reaction_into(model, mod, r_id)
        for r_id in extra_sinks:
            copy_reaction_into(model, ref, r_id)
        model.set_flux_bounds('R_ATPM', *atpm_bounds)
        return model

    def pool_growth_with(extra_sinks):
        model = pool_with(extra_sinks)
        return fba(model, {biomass_id: 1}, env.apply(model, inplace=False, warning=False))

    pool_growth, pool_status = pool_growth_with([])
    required_sinks = []
    sinks = sorted(r for r, why in skipped if why.startswith('sink'))
    if pool_growth < args.min_growth and sinks:
        g_all, _ = pool_growth_with(sinks)
        if g_all >= args.min_growth:
            # Greedy minimisation in fixed (alphabetical) order: drop each sink if growth survives.
            required_sinks = list(sinks)
            for r_id in sinks:
                trial = [s for s in required_sinks if s != r_id]
                if pool_growth_with(trial)[0] >= args.min_growth:
                    required_sinks = trial
            for r_id in required_sinks:
                copy_reaction_into(mod, ref, r_id)
                mod.reactions[r_id].gpr = None
                module.append(r_id)
            skipped = [(r, why) for r, why in skipped if r not in required_sinks]
            pool_growth, pool_status = pool_growth_with([])  # module now contains the sinks
    mod.genes.clear()

    # Per-precursor producibility is only a diagnostic when growth fails: a demand test gives
    # false negatives for carrier metabolites (e.g. charged tRNAs, acyl-ACP) that need recycling.
    precursors = {}
    if pool_growth < args.min_growth:
        check = pool_with([])
        constraints = env.apply(check, inplace=False, warning=False)
        for r_id in module:
            for m_id, coeff in mod.reactions[r_id].stoichiometry.items():
                if coeff >= 0 or m_id in precursors or base_id(m_id) in CURRENCY:
                    continue
                check.add_reaction(CBReaction('R_DMTEST', stoichiometry={m_id: -1}, reversible=False,
                                              lb=0, ub=1000))
                precursors[m_id], _ = fba(check, {'R_DMTEST': 1}, constraints)
                check.remove_reactions(['R_DMTEST'])
    not_producible = sorted(m for m, v in precursors.items() if v < 1e-6)

    validation = {
        'reference_id': prefix, 'reference_file': args.reference, 'reference_sha256': sha, 'biomass': biomass_id,
        'universe': args.universe, 'universe_file': universe_path(args.universe),
        'medium': args.medium, 'media_db': mediadb, 'medium_compounds': compounds,
        'max_uptake': args.max_uptake, 'atpm_source': args.atpm_source, 'atpm_bounds': atpm_bounds,
        'compartment_map': comp_map, 'reference_fields_set_from_fbc': fbc_changes,
        'biomass_module': module, 'alternative_biomass_excluded': excluded_biomass,
        'reference': {'reactions': len(ref.reactions), 'metabolites': len(ref.metabolites),
                      'genes': len(ref.genes)},
        'pool': {'reactions': len(pool.reactions),
                 'reference_only': sum(1 for r in rows if r['in_reference'] and not r['in_universe']),
                 'reference_variants': sum(1 for r in rows if r['note'].startswith('reference variant')),
                 'conflicts_universe_kept': sum(1 for r in rows if 'coefficients only' in r['note']),
                 'reference_reactions_imbalanced_only_in_pool': degraded},
        'reference_reactions_skipped': skipped,
        'checks': {
            'reference_growth_on_medium': ref_growth, 'reference_growth_status': ref_status,
            'pool_growth_on_medium': pool_growth, 'pool_growth_status': pool_status,
            'required_reference_sinks': required_sinks,
            'precursor_diagnostics_not_producible': not_producible,
        },
        'software': software_versions(args.solver),
    }
    write_json(validation, f'{prefix}.validation.json')
    write_json({'biomass': biomass_id, 'module': module, 'atpm_source': args.atpm_source,
                'atpm_bounds': atpm_bounds, 'atpm_reference': atpm['reference'],
                'atpm_carveme': atpm['carveme'],
                'stoichiometry': ref.reactions[biomass_id].stoichiometry}, f'{prefix}.biomass.json')
    with open(f'{prefix}.pool_sources.tsv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['reaction', 'reference_reaction', 'in_universe', 'in_reference',
                                           'balance_in_reference', 'balance_in_pool', 'note'],
                           delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: r['reaction']))
    save_model(ref, f'{prefix}.harmonised.xml')
    save_model(mod, f'{prefix}.biomass_module.xml')
    save_model(pool, f'{prefix}.pool.xml')

    print(f'{prefix}: module={len(module)} rxns; reference growth {ref_growth:.4g} ({ref_status}); '
          f'pool growth {pool_growth:.4g} ({pool_status}); required reference sinks: {required_sinks or "none"}'
          + (f'; non-producible precursors: {not_producible}' if not_producible else ''))
    failures = []
    if ref_growth < args.min_growth:
        failures.append(f'reference does not grow on {args.medium} (growth {ref_growth:.4g})')
    if pool_growth < args.min_growth:
        failures.append(f'pool + biomass module does not grow on {args.medium} (growth {pool_growth:.4g}); '
                        'gap-filling would be infeasible')
    if failures:
        sys.exit('; '.join(failures) + '. See validation JSON.')


if __name__ == '__main__':
    main()
