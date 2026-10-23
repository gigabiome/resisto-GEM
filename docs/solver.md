# MILP solver (CPLEX)

CarveMe's carving and gap-filling are mixed-integer problems with thousands of binary variables. The pipeline is configured for **IBM CPLEX** (`--solver cplex`, the default).

**CPLEX is licensed.** It is not in `envs/resisto.yml` and must never be committed or pushed to a public container registry. The pipeline's first step (`CHECK_SOLVER`) solves a 2,500-variable binary problem and stops the run if the solver is missing, or is a size-limited edition (CPLEX Community Edition: 1,000 variables). This check runs before any genome is processed.

## Option A: Apptainer image (recommended for SLURM)

```bash
bash containers/resisto/build.sh /path/to/CPLEX_Studio_install resisto-gem.sif
```

The script copies the CPLEX Python API for Python 3.10 from your CPLEX Studio folder into the build and installs it on top of `envs/resisto.yml`. If your CPLEX release lays out its Python API differently, put the directory that holds its `setup.py` at `containers/resisto/cplex_python/` and run `apptainer build resisto-gem.sif resisto.def` from `containers/resisto/`.

The CPLEX release must support Python 3.10. Keep the `.sif` file private.

Then run with:

```bash
nextflow run main.nf -profile apptainer --resisto_container /path/to/resisto-gem.sif ...
```

## Option B: conda environment

```bash
conda env create -f envs/resisto.yml -p /path/to/envs/resisto-gem
conda activate /path/to/envs/resisto-gem
cd /path/to/CPLEX_Studio/cplex/python/3.10/x86-64_linux && python -m pip install .
```

The directory above is the usual layout of CPLEX Studio 22.1 but has not been verified for your release. Use the directory that holds the `setup.py` for Python 3.10.

Then run with:

```bash
nextflow run main.nf -profile conda --resisto_conda /path/to/envs/resisto-gem ...
```

## Checking the installation by hand

```bash
python bin/check_solver.py --solver cplex
```

```bash
python tests/regression.py --solver cplex
```

The regression test was written with SCIP. With CPLEX, an alternative but equally optimal gap-fill set is reported as a warning, not an error.

## Other solvers

`--solver gurobi` works if `gurobipy` and a full licence reachable from the compute nodes are installed in the environment. `--solver scip` needs no licence, but:
- CarveMe sets a 600 s limit per MILP for SCIP;
- SCIP accepts suboptimal solutions (the carving MILP was SUBOPTIMAL in testing).

The solver used and every MILP status are recorded per model (`carveme/*.carve.json`, `models/*.build.json`).

## Threads

reframed exposes no thread setting, so CPLEX and Gurobi would otherwise use every core on a node. `bin/resisto_common.py:set_solver` caps them at the CPUs Nextflow allocates to the task.
