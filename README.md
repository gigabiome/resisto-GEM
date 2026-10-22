# resisto-GEM

Strain-specific genome-scale metabolic models (GEMs) for clinical isolates of seven species:
*E. coli*, *K. pneumoniae*, *Enterobacter* spp., *A. baumannii*, *P. aeruginosa*, *S. aureus* and *E. faecium*.

**Status: design phase.** The pipeline is not implemented yet. The current concept and its open decisions are in [docs/figures/gem_pipeline_workflow.pdf](docs/figures/gem_pipeline_workflow.pdf).

Planned flow:
1. Gene prediction with [AMRBench_QC](https://github.com/gigabiome/AMRBench_QC).
2. Draft reconstruction with CarveMe.
3. Biomass replacement and gap-filling against a curated reference GEM for each species.

It must scale from the 140-genome test set to about 12,000 genomes.

## Repository layout

| Path | Contents |
|---|---|
| `data/test_sample_sheet.csv` | Test set: 140 assemblies, 20 per species (`sample,species,fasta`) |
| `data/genomes/` | The assemblies themselves. **Not tracked in git**: copy them in separately. `manifest.csv` is tracked |
| `data/bv_brc_metadata.parquet` | BV-BRC metadata for the genomes |
| `reference_gems/` | Catalogue of curated reference GEMs (`references.tsv`, `README.md`). The model files and papers are **not tracked**; download them with the script below |
| `scripts/` | Download, patch and summary scripts for the reference GEMs |
| `summaries/` | Summary table of the reference GEMs (CSV and PDF) |
| `docs/` | Tool survey (`gem_reconstruction_tools.md`) and the workflow figure |

## Getting the untracked data

```bash
# Reference GEMs and papers (about 240 MB), then the patched Ecoli-GEM
bash scripts/download_reference_gems.sh
python scripts/patch_ecoli_gem.py
```

The genome assemblies must be placed under `data/genomes/<Species>/<sample>.fna`, matching the paths in `data/test_sample_sheet.csv`. Those paths are relative to the repository root, so run tools from there.

## Requirements (current scripts)

- Python with `cobra` (tested with 0.30.0), `pandas` and `matplotlib`
- `curl` and `unzip` for the download script
