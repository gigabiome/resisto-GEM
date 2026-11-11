# Memote Quality Control for Genome-Scale Metabolic Models

This guide explains how to run memote quality assessment on your reconstructed GEMs.

## What is Memote?

Memote (Metabolic Model Testing) is a standard benchmarking tool that checks:
- **Annotation quality** — metabolites and reactions have proper IDs, names, formulas
- **Stoichiometric consistency** — balanced reactions
- **Mass/charge balance** — reactions conserve mass and charge
- **Metabolic feasibility** — models can produce biomass
- **Functional completeness** — pathways for central metabolism present
- **Gap analysis** — identifies blocked metabolic pathways

**Output:** A score from 0-100% indicating model quality.

## Prerequisites

Install memote on the cluster:
```bash
pip install memote
# or via conda
conda install -c bioconda memote
```

Requires:
- `cobra` (constraint-based modeling)
- `python 3.6+`

## Running Memote QC

### Option 1: Python script (more detailed, slower)

```bash
cd /path/to/resisto-GEM
python3 bin/run_memote_qc.py
```

**Output:**
- `results/memote_qc/*.json` — detailed test results for each model
- `results/memote_summary.csv` — summary table (one row per model)

### Option 2: CLI script (faster, more stable)

```bash
bash bin/run_memote_qc_cli.sh
```

Same outputs, but uses memote's command-line interface which is often more reliable.

### Custom results directory

```bash
python3 bin/run_memote_qc.py --results-dir /path/to/results
bash bin/run_memote_qc_cli.sh /path/to/results
```

## Understanding the Results

### Quality Score Interpretation

| Score | Status | Action |
|-------|--------|--------|
| ≥90% | **High quality** | Ready for downstream analysis |
| 70-90% | **Needs review** | Check which tests are failing; may still be usable |
| <70% | **Poor quality** | Requires model refinement or gap-filling adjustment |

### Summary CSV Columns

- `model`: Model ID (biosample)
- `score`: Overall quality score (0-1)
- `passed`: Number of tests passed
- `warnings`: Tests with partial credit (0 < score < 1)
- `failed`: Number of failed tests
- `reactions`: Number of reactions in model
- `metabolites`: Number of metabolites in model
- `genes`: Number of genes with gene-protein-reaction rules
- `biomass_reactions`: Number of biomass reaction variants

### Common Issues and Fixes

**Low metabolite annotation score**
- Many metabolites missing BiGG IDs or formulas
- Impact: Minor (doesn't prevent FBA)
- Fix: Can be addressed in post-processing

**Failed stoichiometry balance**
- Reactions don't balance chemically
- Impact: Serious — violates physical laws
- Fix: Likely due to gap-filling errors; review reference pool

**Cannot produce biomass**
- Model lacks cofactors or biosynthetic pathways
- Impact: Critical — model is non-viable
- Fix: Improve gap-filling, add missing pathways manually

**Blocked reactions**
- Reactions produce/consume metabolites that can't be produced/consumed
- Impact: Reduces model predictive power
- Fix: Identify missing pathway links; resolve with reference GEMs

## Interpreting Detailed Results

Each model gets a JSON report with per-test details:

```json
{
  "score": 0.85,
  "tests": {
    "test_annotation_metabolite_formula": {
      "metric": 0.92,
      "message": "92% of metabolites have a chemical formula"
    },
    "test_find_blocked_reactions": {
      "metric": 0.0,
      "message": "Model has 5 blocked reactions"
    }
  }
}
```

Look for `metric: 0.0` to see which tests are failing and why.

## Next Steps

### For High-Quality Models (≥90%)
- Proceed to metabolic analysis (FBA, phenotype simulation)
- Use for benchmarking or strain engineering studies
- Ready for publication/archive

### For Medium-Quality Models (70-90%)
1. Review failing tests in JSON reports
2. Check if issues are annotation-only (low impact) vs. functional (high impact)
3. Adjust gap-filling parameters if many models fail the same test
4. Consider batch refinement of the full dataset

### For Poor Models (<70%)
1. Investigate whether the problem is:
   - **Gene prediction errors** (check QUAST/CheckM2 scores)
   - **Reference pool mismatch** (species not well-covered by reference GEMs)
   - **Annotation issues** (missing genes or domains in gene prediction)
2. Options:
   - Re-run with improved reference GEM pools
   - Manually curate the model
   - Exclude from benchmarking if severely broken

## Species-Level Summary

After running memote on all models, you can generate a species-level summary:

```bash
python3 << 'EOF'
import pandas as pd

df = pd.read_csv('results/memote_summary.csv')
df['species'] = df['model'].str.split('.').str[0]  # Extract species from model ID

summary = df.groupby('species').agg({
    'score': ['mean', 'min', 'max', 'std'],
    'reactions': 'mean',
    'genes': 'mean',
    'biomass_reactions': 'sum'
}).round(3)

print(summary)
summary.to_csv('results/memote_by_species.csv')
EOF
```

## Troubleshooting

**"ERROR: memote not found"**
- Install: `pip install memote cobra`
- Or: `conda install -c bioconda memote`

**Script hangs on a model**
- Memote can be slow for very large models (>10k reactions)
- Set a timeout or run models in parallel with different workers

**"Cannot read model: Parsing failed"**
- SBML file is corrupted or malformed
- Check with: `cobra.io.read_sbml_model(path)` directly
- May need to re-export from CarveMe

**Low scores across all models**
- Check if the reference GEM is appropriate for your species
- Verify gap-filling worked correctly
- May indicate systematic issue in pipeline

## Running at Scale

For all 12,993 genomes, consider:

1. **Parallelization**: Run on multiple cluster nodes
   ```bash
   # Split models into chunks, each run in parallel job
   ls results/models/*.xml | split -n l/10 - chunk_
   # Submit 10 parallel jobs, each processing one chunk
   ```

2. **Time estimates**: ~5-10 seconds per model
   - 140 models = 10-15 minutes
   - 12,993 models = 18-36 hours with 1 worker
   - Use `--nodes=10` equivalent to parallelize

3. **Storage**: JSON reports are ~50KB each
   - 140 models × 50KB = 7 MB
   - 12,993 models × 50KB = 650 MB

## References

- **Memote docs**: https://memote.readthedocs.io/
- **Memote benchmark scores**: Standard thresholds for published benchmarks
- **CarveMe integration**: Memote complements CarveMe's reconstruction
