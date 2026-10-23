#!/usr/bin/env python
"""Collect per-sample results into one table (report only; nothing is filtered).

Inputs are the files staged by the SUMMARISE process:
  samplesheet.csv              sample, species, fasta (the run's input)
  qc/*.qc.json, build/*.build.json, carve/*.carve.json
  checkm2/<sample>/quality_report.tsv, quast/<sample>/report.tsv
  references/*.validation.json
Outputs:
  resisto_gem_summary.tsv      one row per sample in the samplesheet (missing steps left empty)
  reference_summary.tsv        one row per prepared reference
  software_versions.json       versions reported by the model-building steps (must be identical)
"""
import argparse
import csv
import glob
import json
import os


def load_json(pattern, key):
    out = {}
    for path in glob.glob(pattern):
        with open(path) as f:
            d = json.load(f)
        out[d[key]] = d
    return out


def read_tsv_row(path):
    with open(path) as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    return rows[0] if rows else {}


def quast_metrics(path):
    """QUAST report.tsv is transposed: first column = metric, second = value."""
    metrics = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 2:
                metrics[parts[0]] = parts[1]
    return metrics


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--samplesheet', required=True)
    args = p.parse_args()

    with open(args.samplesheet) as f:
        samples = list(csv.DictReader(f))
    qc = load_json('qc/*.qc.json', 'sample')
    build = load_json('build/*.build.json', 'sample')
    carve = {os.path.basename(p)[:-len('.carve.json')]: json.load(open(p)) for p in glob.glob('carve/*.carve.json')}
    refs = load_json('references/*.validation.json', 'reference_id')

    rows = []
    for s in samples:
        sid = s['sample']
        row = {'sample': sid, 'species': s['species']}
        cm = os.path.join('checkm2', sid, 'quality_report.tsv')
        if os.path.exists(cm):
            r = read_tsv_row(cm)
            row.update({'checkm2_completeness': r.get('Completeness'),
                        'checkm2_contamination': r.get('Contamination')})
        qr = os.path.join('quast', sid, 'report.tsv')
        if os.path.exists(qr):
            q = quast_metrics(qr)
            row.update({'quast_contigs': q.get('# contigs'), 'quast_total_length': q.get('Total length'),
                        'quast_n50': q.get('N50'), 'quast_gc': q.get('GC (%)')})
        c = carve.get(sid, {})
        row['carve_milp_status'] = c.get('milp', {}).get('status')
        b = build.get(sid, {})
        if b:
            row.update({'reference_id': b.get('reference_id'), 'model_id': b.get('model_id'),
                        'draft_reactions': b['draft']['reactions'],
                        'gapfill_milp_status': b['gapfill']['milp'].get('status'),
                        'gapfill_added': b['gapfill']['added'],
                        'gapfill_reference_only': b['gapfill']['from_reference'],
                        'gapfill_reference_and_universe': b['gapfill']['from_reference+universe'],
                        'gapfill_universe_only': b['gapfill']['from_universe'],
                        'reactions': b['final']['reactions'], 'metabolites': b['final']['metabolites'],
                        'genes': b['final']['genes'],
                        'growth_gapfill_medium': b['growth']['gapfill_medium']})
        q = qc.get(sid, {})
        if q:
            row.update({k: q.get(k) for k in ['growth_default', 'growth_matches_build', 'growth_without_carbon',
                                              'egc_atp', 'mass_imbalanced', 'charge_imbalanced',
                                              'gpr_coverage', 'blocked_reactions', 'flags']})
        row['status'] = 'model' if q else ('failed_after_carve' if c else 'no_model')
        rows.append(row)

    fields = []
    for r in rows:
        fields += [k for k in r if k not in fields]
    with open('resisto_gem_summary.tsv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)

    ref_rows = []
    for rid, v in sorted(refs.items()):
        ref_rows.append({'reference_id': rid, 'biomass': v['biomass'], 'universe': v['universe'],
                         'medium': v['medium'], 'atpm_bounds': v['atpm_bounds'],
                         'biomass_module': ';'.join(v['biomass_module']),
                         'required_reference_sinks': ';'.join(v['checks']['required_reference_sinks']),
                         'pool_reactions': v['pool']['reactions'],
                         'pool_reference_only': v['pool']['reference_only'],
                         'pool_reference_variants': v['pool']['reference_variants'],
                         'reference_growth': v['checks']['reference_growth_on_medium'],
                         'pool_growth': v['checks']['pool_growth_on_medium']})
    if ref_rows:
        with open('reference_summary.tsv', 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(ref_rows[0]), delimiter='\t', lineterminator='\n')
            w.writeheader()
            w.writerows(ref_rows)

    versions = {json.dumps(d.get('software', {}), sort_keys=True) for d in list(build.values()) + list(refs.values())}
    with open('software_versions.json', 'w') as f:
        json.dump([json.loads(v) for v in sorted(versions)], f, indent=2)
    if len(versions) > 1:
        raise SystemExit('Model-building steps reported different software versions; see software_versions.json')
    print(f'{len(rows)} samples: ' + ', '.join(f'{k}={sum(r["status"] == k for r in rows)}'
                                             for k in ('model', 'failed_after_carve', 'no_model')))


if __name__ == '__main__':
    main()
