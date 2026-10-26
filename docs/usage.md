# resisto-GEM: usage

## Before the first run (on the cluster)

1. Clone with the AMRBench_QC submodule:

   ```bash
   git clone --recurse-submodules https://github.com/gigabiome/resisto-GEM.git
   ```

   In an existing clone, run `git submodule update --init`. Without the submodule, Nextflow stops with `Can't find a matching module file for include: ../external/AMRBench_QC/...`. Nextflow fetches it automatically when it pulls the pipeline itself (`nextflow run gigabiome/resisto-GEM`), because `manifest.recurseSubmodules = true`.

2. Download the reference GEMs (about 240 MB) into `reference_gems/`:

   ```bash
   bash scripts/download_reference_gems.sh
   ```

   The pipeline checks each reference against the SHA-256 in `assets/references.csv` and stops on a mismatch.

3. Set up the solver environment: Gurobi from the cluster module plus `gurobipy` in `envs/resisto.yml` (`docs/solver.md`).

4. Optionally, check the model-building scripts in that environment:

   ```bash
   python tests/regression.py --solver gurobi
   ```

5. Put the assemblies where the samplesheet points. Relative paths are resolved from the launch directory.

## Samplesheet

```csv
sample,species,fasta
562.30349,Escherichia coli,data/genomes/Escherichia_coli/562.30349.fna
```

- `sample` must be unique and contain no spaces.
- `species` must match a row of `assets/references.csv` exactly: `Escherichia coli`, `Klebsiella pneumoniae`, `Enterobacter spp.`, `Acinetobacter baumannii`, `Pseudomonas aeruginosa`, `Staphylococcus aureus` or `Enterococcus faecium`.
- `fasta` may be gzipped.

The run stops before any task if a sample is duplicated, a file is missing, or a species has no reference.

## Running

```bash
nextflow run main.nf -profile apptainer -c cluster.config \
    --input data/test_sample_sheet.csv \
    --resisto_container /path/to/resisto-gem.sif \
    --bakta_db /path/to/bakta/db --checkm2_db /path/to/CheckM2_database \
    --outdir results
```

- `cluster.config` is your SLURM configuration (executor, queues, scratch, container cache).
- `-profile test` sets `--input data/test_sample_sheet.csv`; launch from the repository root.
- Add `-resume` to continue a stopped run.
- Add `-stub` to check the run's set-up on the cluster without running the tools. The AMRBench_QC modules have no stub blocks, so they run for real even in stub mode.

## Parameters

| Parameter | Default | Meaning |
|---|---|---|
| `--input` | required | Samplesheet |
| `--outdir` | `results` | Output folder |
| `--reference_registry` | `assets/references.csv` | Species → reference GEM, biomass, universe, SHA-256 |
| `--reference_dir` | `reference_gems/` | Where the reference files are |
| `--mediadb` / `--medium` | `assets/media.tsv` / `LB_btn` | Gap-filling medium (CarveMe LB + biotin) |
| `--solver` | `gurobi` | `gurobi`, `cplex` or `scip` |
| `--resisto_modules` | none | Cluster module(s) loaded by the model-building jobs, e.g. `gurobi/12.0.3` |
| `--reference_bonus` | `1.0` | Score added to reference reactions in gap-filling (cost 1/(1+score)) |
| `--carve_reference_score` | `0.0` | CarveMe `--reference-score` (CarveMe default) |
| `--min_growth` | `0.1` | Growth required by gap-filling (CarveMe default) |
| `--gapfill_max_uptake` | `100` | Medium uptake bound while gap-filling (CarveMe CLI) |
| `--init_max_uptake` | `10` | Medium uptake bound stored in the models (CarveMe `-i`) |
| `--atpm_source` | `reference` | ATP maintenance bounds from the reference, or CarveMe's (`carveme`) |
| `--diamond_args` | `--more-sensitive --top 10 --quiet` | CarveMe's DIAMOND defaults; threads are added per task |
| `--qc_blocked` | `true` | Count blocked reactions in QC (FVA) |
| `--resisto_env` / `--resisto_container` / `--resisto_conda` | none / none / `envs/resisto.yml` | Environment for the model-building steps (CarveMe + solver); see `docs/solver.md` |
| `--bakta_db`, `--bakta_db_type`, `--checkm2_db`, `--db_dir`, `--bakta_args` | as AMRBench_QC | Databases are downloaded into `--db_dir` if not given (Bakta full: ~70 GB) |
| `--publish_intermediates` | `false` | Also publish CarveMe drafts (~9 MB) and DIAMOND tables |
| `--bakta_publish_all` | `false` | Publish all Bakta files, not only `faa/gff3/tsv/txt` |
| `--max_retries` | `2` | Retries per task, with memory and time growing each time |
| `--fail_fast` | `false` | Stop when a genome fails all retries; by default it is skipped and listed in the summary |

## Defaults and why

| Choice | Value | Reason |
|---|---|---|
| Gene prediction | Bakta via AMRBench_QC | Same predictions as the AMR benchmark |
| Reconstruction | CarveMe 1.6.6 | Chosen tool; BiGG namespace matches the references |
| References | see `docs/method.md` | iCN718 rather than iACB23LX (whose biomass is CarveMe's generic one); KpSC-pan for *Enterobacter* (no curated *Enterobacter* GEM) |
| Medium | LB + biotin, all species | One medium; biotin is required by iDR478's biomass and is present in real LB |
| QC | report only | Filtering is a downstream decision |
| Solver | Gurobi (cluster module) | A full commercial solver is needed: SCIP returns suboptimal carving solutions |

## Outputs

| Path | Contents |
|---|---|
| `resisto_gem_summary.tsv` | One row per sample: assembly QC, MILP status, gap-fill counts, model size, growth, QC flags, status |
| `reference_summary.tsv`, `references/` | Per reference: harmonised model, biomass module, pool, pool sources, validation |
| `models/<sample>.xml.gz` | Strain model (SBML L3 FBC, BiGG IDs) |
| `models/<sample>.gapfill.tsv` | Gap-filled reactions with source, score, cost and candidate GPR |
| `models/<sample>.build.json` | Build provenance: inputs, parameters, MILP status, growth, software versions |
| `qc/<sample>.qc.json` | COBRApy QC |
| `carveme/<sample>.carve.json` | Carving MILP status and the exact `carve` arguments |
| `bakta/`, `checkm2/`, `quast/` | AMRBench_QC outputs (reduced; see parameters) |
| `software_versions.json` | Versions reported by every model-building step (the run fails if they differ) |
| `pipeline_info/` | Solver check, Nextflow trace, report and timeline |

Load models with COBRApy (`cobra.io.read_sbml_model`) or reframed (`load_cbmodel(..., flavor='bigg')`). If you modify a model with reframed, write it with `bin/resisto_common.py:save_model`, not `save_cbmodel` directly; otherwise COBRApy will misread the result (`docs/method.md`, section 6).
