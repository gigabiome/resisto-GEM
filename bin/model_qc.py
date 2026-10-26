#!/usr/bin/env python
"""Independent quality checks of a strain model with COBRApy (a different SBML parser, and LP
solver GLPK, rather than the reframed/CarveMe build with Gurobi/CPLEX, so it also catches I/O problems).

Checks (all reported; nothing is filtered here):
  sbml_*                 libSBML/COBRApy validation messages per category
  objective              must be exactly the expected reference biomass reaction
  growth_default         FBA growth with the bounds stored in the model (medium applied at build)
  growth_build_reframed  growth the build step computed with reframed (from build.json), for comparison
  growth_without_carbon  growth after closing uptake of every carbon-containing exchange compound;
                         > 0 means biomass can be made without a carbon source ("free lunch")
  egc_atp                max ATP hydrolysis (R_ATPM, bounds 0..1000) with every exchange uptake
                         closed; > 0 means an energy-generating cycle (Fritzemeier et al. 2017)
  mass/charge imbalance  non-boundary reactions (excluding the biomass module) whose elements or
                         charges do not balance, overall and among gap-filled reactions; reactions
                         with a metabolite lacking formula/charge are counted as unknown instead
  blocked_reactions      reactions that cannot carry flux under the stored bounds (optional; FVA)
  gpr_coverage           fraction of non-boundary reactions with a gene rule
"""
import argparse
import csv
import json
import logging

import cobra
from cobra.flux_analysis import find_blocked_reactions

from resisto_common import read_fbc_fields

logging.getLogger('cobra').setLevel(logging.ERROR)
TOL = 1e-6
QC_SOLVER = 'glpk'


def strip(r_id):
    return r_id[2:] if r_id.startswith('R_') else r_id


def has_carbon(met):
    if not met.formula:
        return None
    return met.elements.get('C', 0) > 0


def slim(model):
    value = model.slim_optimize(error_value=float('nan'))
    return 0.0 if value != value else max(value, 0.0)  # NaN (infeasible) -> 0


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--model', required=True)
    p.add_argument('--build-json', required=True)
    p.add_argument('--gapfill-tsv', required=True)
    p.add_argument('--sample', required=True)
    p.add_argument('--blocked', action='store_true', help='Also count blocked reactions (FVA; slower)')
    p.add_argument('--processes', type=int, default=1)
    args = p.parse_args()

    with open(args.build_json) as f:
        build = json.load(f)
    with open(args.gapfill_tsv) as f:
        gapfilled = {strip(r['reaction']) for r in csv.DictReader(f, delimiter='\t')}
    module = {strip(r) for r in build['biomass_module_added']}
    biomass = strip(build['biomass'])

    model, errors = cobra.io.validate_sbml_model(args.model)
    qc = {'sample': args.sample, 'sbml_loaded': model is not None}
    for key, msgs in errors.items():
        qc[f'sbml_{key.lower()}'] = len(msgs)
    if model is None:
        print(json.dumps(qc))
        with open(f'{args.sample}.qc.json', 'w') as f:
            json.dump(qc, f, indent=2)
        raise SystemExit('Model could not be loaded by COBRApy')

    # Pin the LP solver: COBRApy otherwise picks Gurobi/CPLEX automatically when installed (it did,
    # once gurobipy was added), making the QC solver depend on the environment. GLPK keeps the
    # QC independent of the solver used to build the model.
    model.solver = QC_SOLVER
    qc['qc_solver'] = QC_SOLVER

    objective = [r.id for r in model.reactions if r.objective_coefficient != 0]
    qc['objective'] = ';'.join(objective)
    qc['objective_ok'] = objective == [biomass]
    qc['reactions'] = len(model.reactions)
    qc['metabolites'] = len(model.metabolites)
    qc['genes'] = len(model.genes)
    qc['exchanges'] = len(model.exchanges)

    qc['growth_default'] = slim(model)
    qc['growth_build_reframed'] = build['growth']['default_bounds']
    qc['growth_matches_build'] = abs(qc['growth_default'] - qc['growth_build_reframed']) <= \
        1e-4 * max(1.0, abs(qc['growth_build_reframed']))

    with model:
        unknown = 0
        for rxn in model.exchanges:
            met = next(iter(rxn.metabolites))
            carbon = has_carbon(met)
            if carbon is None:
                unknown += 1
            if carbon is not False and rxn.lower_bound < 0:   # unknown formula treated as carbon
                rxn.lower_bound = 0
        qc['growth_without_carbon'] = slim(model)
        qc['exchanges_without_formula'] = unknown

    with model:
        for rxn in model.exchanges:
            rxn.lower_bound = 0
        if 'ATPM' in model.reactions:
            atpm = model.reactions.ATPM
            atpm.bounds = (0, 1000)
            model.objective = atpm
            qc['egc_atp'] = slim(model)
        else:
            qc['egc_atp'] = None

    fbc = read_fbc_fields(args.model)          # species id (M_ prefix) -> (formula, charge)
    known_formula = {k[2:] if k.startswith('M_') else k for k, (f, _) in fbc.items() if f}
    known_charge = {k[2:] if k.startswith('M_') else k for k, (_, c) in fbc.items() if c is not None}
    boundary = {r.id for r in model.boundary}
    internal = [r for r in model.reactions if r.id not in boundary and r.id not in module]
    # A reaction is only judged when every metabolite has an fbc formula (mass) / fbc charge (charge).
    # COBRApy reads a missing fbc:charge as 0, so set-ness is taken from libSBML directly. CarveMe's
    # universe lacks fbc:charge for ~1,000 metabolites, so "unknown" is expected to be large.
    mass_bad, charge_bad, mass_unknown, charge_unknown = [], [], 0, 0
    for rxn in internal:
        mets = list(rxn.metabolites)
        if all(m.id in known_formula for m in mets):
            totals = {}
            for m, coeff in rxn.metabolites.items():
                for el, n in m.elements.items():
                    totals[el] = totals.get(el, 0) + coeff * n
            if any(abs(v) > TOL for v in totals.values()):
                mass_bad.append(rxn.id)
        else:
            mass_unknown += 1
        if all(m.id in known_charge for m in mets):
            if abs(sum(coeff * m.charge for m, coeff in rxn.metabolites.items())) > TOL:
                charge_bad.append(rxn.id)
        else:
            charge_unknown += 1
    qc['mass_imbalanced'] = len(mass_bad)
    qc['mass_unknown'] = mass_unknown
    qc['charge_imbalanced'] = len(charge_bad)
    qc['charge_unknown'] = charge_unknown
    qc['mass_imbalanced_gapfilled'] = sorted(set(mass_bad) & gapfilled)
    qc['charge_imbalanced_gapfilled'] = sorted(set(charge_bad) & gapfilled)
    qc['gpr_coverage'] = sum(1 for r in internal if r.gene_reaction_rule) / max(1, len(internal))
    qc['gapfilled'] = len(gapfilled)

    if args.blocked:
        blocked = find_blocked_reactions(model, processes=args.processes)
        qc['blocked_reactions'] = len(blocked)
        qc['blocked_fraction'] = len(blocked) / max(1, len(model.reactions))

    qc['flags'] = ';'.join(name for name, bad in [
        ('objective', not qc['objective_ok']),
        ('no_growth', qc['growth_default'] < TOL),
        ('growth_mismatch', not qc['growth_matches_build']),
        ('free_lunch', qc['growth_without_carbon'] > TOL),
        ('egc', (qc['egc_atp'] or 0) > TOL),
        ('imbalanced_gapfill', bool(qc['mass_imbalanced_gapfilled'] or qc['charge_imbalanced_gapfilled'])),
    ] if bad)

    with open(f'{args.sample}.qc.json', 'w') as f:
        json.dump(qc, f, indent=2, default=str)
    print(f"{args.sample}: growth {qc['growth_default']:.4g}; flags: {qc['flags'] or 'none'}")


if __name__ == '__main__':
    main()
