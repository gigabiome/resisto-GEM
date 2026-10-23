# resisto-GEM

Strain-specific genome-scale metabolic models (GEMs) for clinical isolates of seven species:
*E. coli*, *K. pneumoniae*, *Enterobacter* spp., *A. baumannii*, *P. aeruginosa*, *S. aureus* and *E. faecium*.

A Nextflow pipeline built for SLURM clusters, scaling from the 140-genome test set to about 12,000 genomes:

1. **Gene prediction and assembly QC** with the [AMRBench_QC](https://github.com/gigabiome/AMRBench_QC) modules (QUAST, CheckM2, Bakta), included as a pinned git submodule.
2. **Draft reconstruction** with CarveMe 1.6.6, guided by the species reference GEM (`--reference`).
3. **Biomass replacement and gap-filling.** The draft gets the reference GEM's biomass, then is gap-filled from the CarveMe universe plus the reference's own reactions, preferring reference reactions.
4. **Independent QC** with COBRApy, and one summary table for the whole run.

**Status:** implemented and checked as follows:
- the scripts are tested on all six references and on CarveMe's E. coli benchmark (`tests/regression.py`);
- the Nextflow wiring is checked with lint, preview and a stub run of all 140 test genomes.

It has **not yet been run on real genomes**; the first run will be on the cluster. See [docs/method.md](docs/method.md) for exactly what it does and why.

## Quick start (cluster)

```bash
git clone --recurse-submodules https://github.com/gigabiome/resisto-GEM.git && cd resisto-GEM
bash scripts/download_reference_gems.sh                      # reference GEMs (checksummed)
bash containers/resisto/build.sh /path/to/CPLEX_Studio       # CarveMe + CPLEX image (docs/solver.md)
nextflow run main.nf -profile apptainer,test -c cluster.config \
    --resisto_container containers/resisto/resisto-gem.sif --outdir results
```

Details are in [docs/usage.md](docs/usage.md).

## Repository layout

| Path | Contents |
|---|---|
| `main.nf`, `workflows/`, `subworkflows/`, `modules/local/`, `conf/`, `nextflow.config` | The pipeline |
| `external/AMRBench_QC` | Git submodule, pinned at `52001eb`; its QUAST/CheckM2/Bakta modules are used unchanged |
| `bin/` | Model-building scripts run by the pipeline (reference preparation, carve wrapper, build, QC, summary) |
| `assets/references.csv` | Species → reference GEM, biomass reaction, CarveMe universe, SHA-256 |
| `assets/media.tsv` | Gap-filling medium `LB_btn` (CarveMe LB + biotin) |
| `envs/resisto.yml`, `containers/resisto/` | Pinned environment and Apptainer recipe (CPLEX added at build time) |
| `tests/regression.py` | Regression test of the model-building scripts |
| `data/test_sample_sheet.csv` | Test set: 140 assemblies, 20 per species. The assemblies in `data/genomes/` are **not tracked in git** |
| `reference_gems/` | Reference GEM catalogue. Models and papers are **not tracked**; run `scripts/download_reference_gems.sh` |
| `scripts/` | Download, patch and summary scripts for the reference GEMs |
| `summaries/`, `docs/` | Reference GEM summary; tool survey, method, usage, solver set-up, workflow figure |

## Requirements

- Nextflow ≥ 24.04 (tested with 26.04.7), Java 17+
- Apptainer/Singularity, Docker or conda
- A CPLEX licence for production runs (`docs/solver.md`)
