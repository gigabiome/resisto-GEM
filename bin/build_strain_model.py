#!/usr/bin/env python
"""Turn a CarveMe draft into a strain model with the species reference biomass, gap-filled
against the prepared reference pool.

Steps:
  1. Load the CarveMe draft; remove CarveMe's generic biomass ('Growth').
  2. Add the reference biomass module (biomass reaction, its lumping/drain reactions and any
     required reference sinks; see prepare_reference.py) and make the biomass the only objective.
  3. Set R_ATPM bounds from <ref>.biomass.json (reference NGAM by default).
  4. Gap-fill with CarveMe's own gapFill() MILP (same algorithm and defaults as `carve -g`), using
     the pool (CarveMe universe + reference reactions) instead of the universe alone:
       cost of adding reaction r = 1 / (1 + score_r)
       score_r = CarveMe homology score of r for this genome (from CarveMe's DIAMOND output and
                 BiGG GPRs, i.e. reaction_scoring()), 0 if unscored
                 + reference_bonus if r is a reference reaction.
     Medium = --medium, uptake bound --gapfill-max-uptake (CarveMe CLI uses 100), growth >= --min-growth.
  5. Apply the medium with --init-max-uptake (CarveMe `-i` default 10) as the model's default bounds.
  6. Write <sample>.xml.gz (SBML L3 FBC, BiGG IDs; gzip-compressed, readable by COBRApy/reframed).
Gap-filled reactions get no GPR (as in CarveMe); the homology-based candidate GPR, if any, is
listed in <sample>.gapfill.tsv together with source, score and cost. Run details: <sample>.build.json.
"""
import argparse
import csv
import json
import os
import sys

import pandas as pd
from carveme import config, project_dir
from carveme.cli.carve import build_model_id
from carveme.reconstruction.diamond import load_diamond_results
from carveme.reconstruction.gapfilling import gapFill, merge_models
from carveme.reconstruction.scoring import reaction_scoring
from reframed import FBA, Environment, load_cbmodel
from reframed.core.transformation import disconnected_metabolites
from reframed.solvers import solver_instance

from resisto_common import (CARVEME_BIOMASS, SOURCE_TAG, copy_reaction_into, load_medium,
                            save_model, set_solver, software_versions, write_json)


def growth(model, biomass_id, constraints=None):
    sol = FBA(model, objective={biomass_id: 1}, constraints=constraints, get_values=False)
    ok = sol.status.name in ('OPTIMAL', 'SUBOPTIMAL')
    return (sol.fobj if ok and sol.fobj is not None else 0.0), sol.status.name


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--sample', required=True)
    p.add_argument('--species', required=True)
    p.add_argument('--reference-id', required=True)
    p.add_argument('--draft', required=True, help='CarveMe draft model (SBML)')
    p.add_argument('--diamond', required=True, help='CarveMe DIAMOND output (.tsv) for this genome')
    p.add_argument('--biomass-module', required=True)
    p.add_argument('--biomass-json', required=True)
    p.add_argument('--pool', required=True)
    p.add_argument('--pool-sources', required=True)
    p.add_argument('--medium', required=True)
    p.add_argument('--mediadb', default=None)
    p.add_argument('--min-growth', type=float, default=0.1)
    p.add_argument('--gapfill-max-uptake', type=float, default=float(config.getint('gapfill', 'max_uptake')))
    p.add_argument('--init-max-uptake', type=float, default=10.0)
    p.add_argument('--reference-bonus', type=float, default=1.0)
    p.add_argument('--solver', default='cplex')
    p.add_argument('--threads', type=int, default=1)
    args = p.parse_args()

    set_solver(args.solver, args.threads)
    out = args.sample

    # 1-2. draft and biomass ----------------------------------------------------------------------
    model = load_cbmodel(args.draft, flavor='bigg')
    if model.biomass_reaction != CARVEME_BIOMASS:
        sys.exit(f'Draft objective is {model.biomass_reaction!r}, expected CarveMe {CARVEME_BIOMASS!r}')
    if 'R_ATPM' not in model.reactions:
        sys.exit('Draft has no R_ATPM')
    draft_stats = {'reactions': len(model.reactions), 'metabolites': len(model.metabolites),
                   'genes': len(model.genes)}
    model.remove_reactions([CARVEME_BIOMASS])

    with open(args.biomass_json) as f:
        biomass_info = json.load(f)
    biomass_id = biomass_info['biomass']
    module = load_cbmodel(args.biomass_module, flavor='bigg')
    added_module = []
    for r_id in biomass_info['module']:
        if r_id in model.reactions:
            if model.reactions[r_id].stoichiometry != module.reactions[r_id].stoichiometry:
                sys.exit(f'Biomass-module reaction {r_id} is in the draft with a different equation')
            continue
        rxn = copy_reaction_into(model, module, r_id)
        rxn.metadata[SOURCE_TAG] = 'biomass_module'
        added_module.append(r_id)
    for r_id, rxn in model.reactions.items():
        rxn.objective = 1 if r_id == biomass_id else 0
    model.biomass_reaction = biomass_id

    # 3. ATPM -------------------------------------------------------------------------------------
    lb, ub = biomass_info['atpm_bounds']
    model.set_flux_bounds('R_ATPM', lb, ub)

    # 4. scores and gap-filling -----------------------------------------------------------------
    pool = load_cbmodel(args.pool, flavor='bigg')
    sources = {}
    with open(args.pool_sources) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            sources[row['reaction']] = row

    gprs = pd.read_csv(project_dir + config.get('generated', 'bigg_gprs'))
    gprs = gprs[gprs.reaction.isin(list(pool.reactions))]
    annotations = load_diamond_results(args.diamond)
    reaction_scores, _ = reaction_scoring(annotations, gprs)
    homology = {} if reaction_scores is None else dict(reaction_scores[['reaction', 'normalized_score']].values)
    homology_gpr = {} if reaction_scores is None else dict(reaction_scores[['reaction', 'GPR']].values)

    def is_reference(r_id):
        return sources.get(r_id, {}).get('in_reference') == '1'

    scores = {r_id: homology.get(r_id, 0.0) + (args.reference_bonus if is_reference(r_id) else 0.0)
              for r_id in pool.reactions}

    compounds, mediadb = load_medium(args.medium, args.mediadb)
    before = set(model.reactions)

    # Same preparation as carveme.reconstruction.gapfilling.multiGapFill (CarveMe 1.6.6), with a
    # single medium and the solver handle kept so that the MILP status can be recorded.
    for r_id in set(pool.reactions) - before:
        if r_id.startswith('R_EX'):
            pool.set_flux_bounds(r_id, lb=0)
    merged = merge_models(model, pool, inplace=False)
    solver = solver_instance(merged)
    milp = {}
    original_solve = solver.solve

    def solve(*a, **kw):
        sol = original_solve(*a, **kw)
        milp.update(status=sol.status.name, message=getattr(sol, 'message', None), objective=sol.fobj)
        return sol
    solver.solve = solve
    constraints = Environment.from_compounds(compounds, max_uptake=args.gapfill_max_uptake) \
        .apply(merged, inplace=False, warning=False)
    try:
        gapFill(model, pool, constraints=constraints, min_growth=args.min_growth, scores=scores,
                inplace=True, solver=solver, tag=args.medium)
    except RuntimeError as e:
        write_json({'sample': out, 'error': str(e), 'milp': milp}, f'{out}.build.json')
        sys.exit(f'Gap-filling failed: {e}')

    added = sorted(set(model.reactions) - before)
    gap_rows = []
    for r_id in added:
        rxn = model.reactions[r_id]
        src = sources.get(r_id, {})
        if src.get('in_reference') == '1' and src.get('in_universe') == '1':
            origin = 'reference+universe'
        elif src.get('in_reference') == '1':
            origin = 'reference'
        else:
            origin = 'universe'
        rxn.metadata[SOURCE_TAG] = f'gapfill:{origin}'
        gap_rows.append({'reaction': r_id, 'source': origin,
                         'homology_score': homology.get(r_id, 0.0), 'score': scores.get(r_id, 0.0),
                         'cost': 1.0 / (1.0 + scores.get(r_id, 0.0)),
                         'homology_gpr': homology_gpr.get(r_id) or '',
                         'equation': str(rxn).split(': ', 1)[-1]})

    gf_growth, gf_status = growth(model, biomass_id, Environment.from_compounds(
        compounds, max_uptake=args.gapfill_max_uptake).apply(model, inplace=False, warning=False))
    if gf_growth < args.min_growth - 1e-9:
        sys.exit(f'Gap-filled model grows {gf_growth:.4g} < {args.min_growth} on {args.medium}')

    # 5. default medium ---------------------------------------------------------------------------
    Environment.from_compounds(compounds, max_uptake=args.init_max_uptake).apply(model, inplace=True, warning=False)
    model.remove_metabolites(disconnected_metabolites(model))
    final_growth, final_status = growth(model, biomass_id)

    model.id = build_model_id(out)
    save_model(model, f'{out}.xml.gz')

    with open(f'{out}.gapfill.tsv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['reaction', 'source', 'homology_score', 'score', 'cost', 'homology_gpr',
                                           'equation'],
                           delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(gap_rows)

    count = lambda s: sum(1 for r in gap_rows if r['source'] == s)
    write_json({
        'sample': out, 'species': args.species, 'reference_id': args.reference_id, 'model_id': model.id,
        'draft': draft_stats, 'biomass': biomass_id, 'biomass_module_added': added_module,
        'atpm_bounds': [lb, ub], 'medium': args.medium, 'media_db': mediadb,
        'min_growth': args.min_growth, 'gapfill_max_uptake': args.gapfill_max_uptake,
        'init_max_uptake': args.init_max_uptake, 'reference_bonus': args.reference_bonus,
        'scored_reactions': len(homology),
        'gapfill': {'milp': milp, 'added': len(gap_rows), 'from_reference': count('reference'),
                    'from_reference+universe': count('reference+universe'),
                    'from_universe': count('universe')},
        'growth': {'gapfill_medium': gf_growth, 'gapfill_medium_status': gf_status,
                   'default_bounds': final_growth, 'default_bounds_status': final_status},
        'final': {'reactions': len(model.reactions), 'metabolites': len(model.metabolites),
                  'genes': len(model.genes)},
        'inputs': {k: os.path.basename(v) for k, v in vars(args).items()
                   if k in ('draft', 'diamond', 'biomass_module', 'biomass_json', 'pool', 'pool_sources', 'mediadb') and v},
        'software': software_versions(args.solver),
    }, f'{out}.build.json')
    print(f'{out}: added {len(gap_rows)} reactions ({count("reference")} reference-only, '
          f'{count("reference+universe")} shared, {count("universe")} universe-only); '
          f'growth {gf_growth:.4g} on {args.medium}, {final_growth:.4g} with default bounds; MILP {milp.get("status")}')


if __name__ == '__main__':
    main()
