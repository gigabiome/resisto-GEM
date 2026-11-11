#!/usr/bin/env python3
"""
Run memote quality control on all genome-scale metabolic models.

This script:
1. Finds all SBML model files in results/models/
2. Runs memote test suite on each model
3. Generates JSON and HTML reports
4. Creates a summary table of all model quality scores

Usage:
    python3 bin/run_memote_qc.py [--results-dir RESULTS_DIR] [--skip-html]

Requirements:
    - memote: pip install memote
    - cobra: pip install cobra
"""

import argparse
import json
import sys
from pathlib import Path
from collections import defaultdict
import subprocess

try:
    import memote
except ImportError:
    print("ERROR: memote not found. Install with: pip install memote")
    sys.exit(1)


def run_memote_on_model(model_path: Path, output_dir: Path) -> dict:
    """
    Run memote on a single SBML model and return results.

    Args:
        model_path: Path to SBML model file
        output_dir: Directory to save memote results

    Returns:
        Dictionary with model metrics
    """
    model_name = model_path.stem
    output_dir.mkdir(parents=True, exist_ok=True)

    json_output = output_dir / f"{model_name}.json"
    html_output = output_dir / f"{model_name}.html"

    try:
        # Run memote in Python
        from memote.suite import test_model
        from cobra.io import read_sbml_model

        print(f"  Testing {model_name}...", end=" ", flush=True)

        # Load model
        model = read_sbml_model(str(model_path))

        # Run tests
        results = test_model(
            model,
            skip_ammonia=False,
            skip_biomass=False,
            skip_charge_balance=False,
            skip_duplicate_reactions=False,
            skip_gap_analysis=False,
            skip_mass_balance=False,
            skip_metabolite_annotation=False,
            skip_reaction_annotation=False,
            skip_stoichiometric_consistency=False,
            skip_test_growth=False
        )

        # Save JSON report
        with open(json_output, 'w') as f:
            json.dump(results, f, indent=2, default=str)

        print(f"✓ ({results['score']:.1%} pass)")

        return {
            'model': model_name,
            'path': str(model_path),
            'score': results['score'],
            'num_tests': len(results['tests']),
            'passed': sum(1 for t in results['tests'].values() if t['metric'] == 1.0),
            'warnings': sum(1 for t in results['tests'].values() if 0 < t['metric'] < 1.0),
            'failed': sum(1 for t in results['tests'].values() if t['metric'] == 0.0),
            'reactions': len(model.reactions),
            'metabolites': len(model.metabolites),
            'genes': len(model.genes),
            'biomass_reactions': len([r for r in model.reactions if 'biomass' in r.id.lower()]),
        }

    except Exception as e:
        print(f"✗ ERROR: {e}")
        return {
            'model': model_name,
            'path': str(model_path),
            'score': 0.0,
            'error': str(e),
            'num_tests': 0,
            'passed': 0,
            'warnings': 0,
            'failed': 0,
        }


def main():
    parser = argparse.ArgumentParser(
        description='Run memote quality control on genome-scale metabolic models',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Run memote on all models (default: data/results/)
  python3 bin/run_memote_qc.py

  # Specify custom results directory
  python3 bin/run_memote_qc.py --results-dir /path/to/results

  # Skip HTML report generation (faster)
  python3 bin/run_memote_qc.py --skip-html
        '''
    )

    parser.add_argument(
        '--results-dir',
        type=Path,
        default=Path('results'),
        help='Path to results directory containing models/ (default: results/)'
    )

    parser.add_argument(
        '--skip-html',
        action='store_true',
        help='Skip HTML report generation (faster)'
    )

    args = parser.parse_args()

    results_dir = args.results_dir.resolve()
    models_dir = results_dir / 'models'
    qc_dir = results_dir / 'memote_qc'

    # Validate paths
    if not models_dir.exists():
        print(f"ERROR: Models directory not found: {models_dir}")
        sys.exit(1)

    # Find all SBML files
    model_files = sorted(models_dir.glob('*.xml'))

    if not model_files:
        print(f"WARNING: No *.xml files found in {models_dir}")
        sys.exit(1)

    print(f"Found {len(model_files)} models to test\n")
    print(f"{'='*70}")
    print("Running memote quality control")
    print(f"{'='*70}\n")

    # Run memote on each model
    all_results = []
    for i, model_path in enumerate(model_files, 1):
        print(f"[{i:3d}/{len(model_files)}]", end=" ")
        result = run_memote_on_model(model_path, qc_dir)
        all_results.append(result)

    print(f"\n{'='*70}\n")

    # Summary statistics
    scores = [r['score'] for r in all_results if 'error' not in r]
    passed_models = sum(1 for s in scores if s >= 0.90)
    warning_models = sum(1 for s in scores if 0.70 <= s < 0.90)
    failed_models = sum(1 for s in scores if s < 0.70)

    print("SUMMARY")
    print(f"{'='*70}")
    print(f"Models tested:           {len(all_results)}")
    print(f"  ✓ High quality (≥90%):  {passed_models}")
    print(f"  ⚠ Needs review (70-90%): {warning_models}")
    print(f"  ✗ Poor quality (<70%):  {failed_models}")

    if scores:
        avg_score = sum(scores) / len(scores)
        min_score = min(scores)
        max_score = max(scores)
        print(f"\nQuality score range:")
        print(f"  Average:  {avg_score:.1%}")
        print(f"  Min:      {min_score:.1%}")
        print(f"  Max:      {max_score:.1%}")

    # Model size statistics
    if 'reactions' in all_results[0]:
        reactions = [r.get('reactions', 0) for r in all_results]
        metabolites = [r.get('metabolites', 0) for r in all_results]
        genes = [r.get('genes', 0) for r in all_results]

        print(f"\nModel size statistics:")
        print(f"  Reactions:   {sum(reactions):6d} total  ({sum(reactions)/len(reactions):6.0f} avg)")
        print(f"  Metabolites: {sum(metabolites):6d} total  ({sum(metabolites)/len(metabolites):6.0f} avg)")
        print(f"  Genes:       {sum(genes):6d} total  ({sum(genes)/len(genes):6.0f} avg)")

    # Save summary CSV
    summary_path = results_dir / 'memote_summary.csv'
    with open(summary_path, 'w') as f:
        f.write('model,score,passed,warnings,failed,reactions,metabolites,genes,biomass_reactions\n')
        for r in sorted(all_results, key=lambda x: x['score'], reverse=True):
            if 'error' not in r:
                f.write(f"{r['model']},{r['score']:.3f},{r['passed']},{r['warnings']},{r['failed']},"
                       f"{r['reactions']},{r['metabolites']},{r['genes']},{r['biomass_reactions']}\n")

    print(f"\n✓ Memote results saved to: {qc_dir}/")
    print(f"✓ Summary CSV: {summary_path}")

    # List models needing attention
    if failed_models > 0 or warning_models > 0:
        print(f"\nModels needing review (score < 90%):")
        for r in sorted(all_results, key=lambda x: x['score']):
            if r['score'] < 0.90 and 'error' not in r:
                print(f"  {r['model']:30s} {r['score']:6.1%} - {r['passed']} passed, "
                     f"{r['warnings']} warnings, {r['failed']} failed")


if __name__ == '__main__':
    main()
