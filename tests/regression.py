#!/usr/bin/env python
"""Regression test of the model-building scripts (no genomes, no Nextflow).

Run inside the resisto environment from the repository root, after
scripts/download_reference_gems.sh:

    python tests/regression.py --solver cplex      # or scip (no licence needed)

Checks the behaviours established when the scripts were written (2026-10-22, SCIP):
  * every registry reference prepares, grows on LB_btn, matches its SHA-256, and needs exactly the
    expected reference sinks; pool files carry no genes and no boundaryCondition species
  * a CarveMe 1.6.6 E. coli draft (CarveMe's packaged K-12 benchmark model + DIAMOND table) is
    gap-filled with iML1515 reactions only, and the expected Fe-S/lipid reactions are added
  * COBRApy and reframed give the same growth for the result; no free lunch, no ATP cycle
Alternative optimal gap-fill sets are possible with another solver; those differences are
reported as warnings, everything else is an error.
"""
import argparse
import csv
import json
import os
import subprocess
import sys
import tempfile

import carveme

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, 'bin')
EXPECTED_SINKS = {'iML1515': [], 'KpSC_pan_v2': [], 'iCN718': [], 'iSD1509': [],
                  'iYS854': ['R_SK_2fe1s_c', 'R_SK_apoACP_c'], 'iDR478': ['R_SK_apoACP_c']}
EXPECTED_ECOLI_GAPFILL = {'R_EX_ni2_e', 'R_K2L4Aabctex', 'R_S2FE2SS', 'R_S2FE2SS2', 'R_S2FE2ST',
                          'R_S4FE4SR', 'R_S4FE4ST', 'R_SCYSDS', 'R_UDCPDPS'}
errors, warnings = [], []


def run(cmd, cwd):
    r = subprocess.run([sys.executable] + cmd, cwd=cwd, capture_output=True, text=True,
                       env={**os.environ, 'PYTHONPATH': BIN})
    if r.returncode != 0:
        errors.append(f'{os.path.basename(cmd[0])} failed: {r.stderr.strip()[-500:]}')
    return r


def check(cond, msg):
    if not cond:
        errors.append(msg)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--solver', default='scip')
    p.add_argument('--reference-dir', default=os.path.join(ROOT, 'reference_gems'))
    p.add_argument('--keep', action='store_true', help='keep the temporary directory')
    args = p.parse_args()

    tmp = tempfile.mkdtemp(prefix='resisto_regression_')
    media = os.path.join(ROOT, 'assets', 'media.tsv')
    with open(os.path.join(ROOT, 'assets', 'references.csv')) as f:
        registry = {r['reference_id']: r for r in csv.DictReader(f)}

    for rid, r in registry.items():
        run([os.path.join(BIN, 'prepare_reference.py'), '--reference', os.path.join(args.reference_dir, r['reference_file']),
             '--reference-id', rid, '--sha256', r['sha256'], '--biomass', r['biomass'], '--universe', r['universe'],
             '--medium', 'LB_btn', '--mediadb', media, '--solver', args.solver], tmp)
        vpath = os.path.join(tmp, f'{rid}.validation.json')
        if not os.path.exists(vpath):
            continue
        v = json.load(open(vpath))
        check(v['checks']['reference_growth_on_medium'] > 0.1, f'{rid}: reference does not grow')
        check(v['checks']['pool_growth_on_medium'] > 0.1, f'{rid}: pool does not grow')
        check(v['checks']['required_reference_sinks'] == EXPECTED_SINKS[rid],
              f"{rid}: required sinks {v['checks']['required_reference_sinks']} != {EXPECTED_SINKS[rid]}")
        for kind in ('pool', 'harmonised', 'biomass_module'):
            text = open(os.path.join(tmp, f'{rid}.{kind}.xml')).read()
            check('boundaryCondition="true"' not in text, f'{rid}.{kind}.xml has boundaryCondition species')
        check('geneProduct' not in open(os.path.join(tmp, f'{rid}.pool.xml')).read(), f'{rid}: pool carries genes')

    bench = os.path.join(carveme.project_dir, 'data', 'benchmark')
    run([os.path.join(BIN, 'build_strain_model.py'), '--sample', 'ecoli_bench', '--species', 'Escherichia coli',
         '--reference-id', 'iML1515', '--draft', os.path.join(bench, 'models', 'ecol.xml'),
         '--diamond', os.path.join(bench, 'fasta', 'Ecoli_K12_MG1655.tsv'),
         '--biomass-module', 'iML1515.biomass_module.xml', '--biomass-json', 'iML1515.biomass.json',
         '--pool', 'iML1515.pool.xml', '--pool-sources', 'iML1515.pool_sources.tsv',
         '--medium', 'LB_btn', '--mediadb', media, '--solver', args.solver], tmp)
    build_path = os.path.join(tmp, 'ecoli_bench.build.json')
    if os.path.exists(build_path):
        b = json.load(open(build_path))
        check(b['gapfill']['milp']['status'] == 'OPTIMAL', f"gap-fill MILP {b['gapfill']['milp']}")
        with open(os.path.join(tmp, 'ecoli_bench.gapfill.tsv')) as f:
            added = {r['reaction']: r['source'] for r in csv.DictReader(f, delimiter='\t')}
        check(all(s == 'reference' for s in added.values()), f'non-reference gap-fill: {added}')
        if set(added) != EXPECTED_ECOLI_GAPFILL:
            warnings.append(f'E. coli gap-fill differs from the SCIP reference set: '
                            f'+{sorted(set(added) - EXPECTED_ECOLI_GAPFILL)} -{sorted(EXPECTED_ECOLI_GAPFILL - set(added))}')
        run([os.path.join(BIN, 'model_qc.py'), '--model', 'ecoli_bench.xml.gz', '--build-json', 'ecoli_bench.build.json',
             '--gapfill-tsv', 'ecoli_bench.gapfill.tsv', '--sample', 'ecoli_bench'], tmp)
        qpath = os.path.join(tmp, 'ecoli_bench.qc.json')
        if os.path.exists(qpath):
            q = json.load(open(qpath))
            check(q['objective_ok'], 'objective is not the iML1515 biomass')
            check(q['growth_matches_build'], f"COBRApy {q['growth_default']} vs reframed {q['growth_build_reframed']}")
            check(q['growth_without_carbon'] < 1e-6, 'free lunch')
            check((q['egc_atp'] or 0) < 1e-6, 'energy-generating cycle')
            check(q.get('sbml_sbml_error', 0) == 0, f"SBML errors: {q.get('sbml_sbml_error')}")

    for w in warnings:
        print('WARNING:', w)
    for e in errors:
        print('ERROR:', e)
    print(f'{"FAILED" if errors else "PASSED"} ({len(errors)} errors, {len(warnings)} warnings); files in {tmp}')
    if not args.keep and not errors:
        subprocess.run(['rm', '-rf', tmp])
    sys.exit(1 if errors else 0)


if __name__ == '__main__':
    main()
