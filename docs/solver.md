# MILP solver (Gurobi)

CarveMe's carving and gap-filling are mixed-integer problems with thousands of binary variables. The pipeline uses **Gurobi** by default (`--solver gurobi`), because it is available on the cluster. CPLEX (`--solver cplex`, sections at the end) and SCIP (`--solver scip`) also work.

**The solver must be the full version.** The pipeline's first step (`CHECK_SOLVER`) solves a 2,500-variable binary problem and stops the run if it can't. That catches a missing solver, or the size-limited licence that comes with `pip install gurobipy`: tested, it fails with "Model too large for size-limited license". This check runs before any genome is processed.

## Gurobi on the cluster (recommended)

1. **Find the cluster's Gurobi module and version:**

   ```bash
   module avail gurobi
   module load gurobi/<version>
   gurobi_cl --version
   ```

   The module normally also sets the licence, e.g. `GRB_LICENSE_FILE` pointing at a licence server. If `gurobi_cl` reports a licence error, ask HPC support how jobs should obtain the licence.

2. **Create the pipeline's Python environment and add `gurobipy` at the same version as the module.** Its licence must cover the `gurobipy` version, so match the module's major.minor, e.g. module 12.0.3 → `gurobipy==12.0.3`.

   ```bash
   conda env create -f envs/resisto.yml -p $HOME/envs/resisto-gem
   conda activate $HOME/envs/resisto-gem
   pip install "gurobipy==<module version>"
   ```

3. **Check it** with the module loaded:

   ```bash
   module load gurobi/<version>
   python bin/check_solver.py --solver gurobi
   python tests/regression.py --solver gurobi
   ```

   The regression test was written with SCIP. Any difference in Gurobi's gap-fill set is an equally optimal alternative and is reported as a warning, not an error.

4. **Run the pipeline:**

   ```bash
   nextflow run main.nf -profile apptainer -c cluster.config \
       --resisto_env $HOME/envs/resisto-gem --resisto_modules gurobi/<version> ...
   ```

   - `--resisto_modules` makes every model-building job run `module load gurobi/<version>`, so the licence settings reach the compute nodes. Tested: the job scripts contain the `module load` line.
   - `--resisto_env` puts the environment's `bin/` first on `PATH` for those jobs only. The AMRBench_QC steps keep their Apptainer containers.
   - The licence server must be reachable from the compute nodes.

# CPLEX (alternative)

Use `--solver cplex` with one of the options below.

## CPLEX option A: Apptainer image (CPLEX ≤ 22.1.1 only)

```bash
bash containers/resisto/build.sh /path/to/CPLEX_Studio_install resisto-gem.sif
```

The script copies the CPLEX Python API for Python 3.10 from your CPLEX Studio folder into the build and installs it on top of `envs/resisto.yml`. If your CPLEX release lays out its Python API differently, put the directory that holds its `setup.py` at `containers/resisto/cplex_python/` and run `apptainer build resisto-gem.sif resisto.def` from `containers/resisto/`.

The CPLEX release must support Python 3.10. Keep the `.sif` file private.

**Limitations:** `build.sh` needs the `setup.py` that CPLEX Studio shipped up to 22.1.1; IBM removed it in 22.1.2, so use Option B for newer releases. Building an image also needs `apptainer build` rights on the cluster (often `--fakeroot`).

Then run with:

```bash
nextflow run main.nf -profile apptainer --solver cplex --resisto_container /path/to/resisto-gem.sif ...
```

## CPLEX option B: your own environment on the cluster (no root needed)

Everything goes in your home (or project) folder. The AMRBench_QC steps keep running in Apptainer.

1. **Check first whether your cluster already provides CPLEX** (e.g. `module avail cplex`, or ask HPC support). If it does, skip step 2 and use that installation folder in step 4.

2. **Get CPLEX Optimization Studio and install it in your home folder.**
   - Academic users get it free from IBM SkillsBuild (https://academic.ibm.com/a2mt/). Register with your institutional e-mail, then go to Data Science, ILOG CPLEX Optimization Studio, download method HTTP, and choose the **Linux x86-64** installer.
   - Copy the installer to the cluster and run it on the login node:

     ```bash
     sh cplex_studio2220.linux_x86_64.bin
     ```

   - When the installer asks where to install, give a folder you own, e.g. `$HOME/opt/ibm/ILOG/CPLEX_Studio2220`. The default under `/opt` needs root.
   - Adjust `2220` to your release (e.g. `2212` for 22.1.2).

3. **Create the pipeline's Python environment** (with conda, mamba or micromamba):

   ```bash
   conda env create -f envs/resisto.yml -p $HOME/envs/resisto-gem
   ```

4. **Add the CPLEX Python API, then switch it to your full licence.** On PyPI, `pip install cplex` gives the size-limited Community Edition (1,000 variables). `docplex config --upgrade` replaces its libraries with those of your installation.

   ```bash
   conda activate $HOME/envs/resisto-gem
   pip install "cplex==22.2.0.1" "docplex==2.32.264"
   docplex config --upgrade $HOME/opt/ibm/ILOG/CPLEX_Studio2220
   ```

   Match the `cplex` version to your CPLEX Studio release. 22.2.0 → `cplex==22.2.0.1`; 22.1.2 → `cplex==22.1.2.1`. Both have Python 3.10 wheels.

5. **Check it** (the solver check fails on the Community Edition, so it confirms the upgrade worked):

   ```bash
   python bin/check_solver.py --solver cplex
   python tests/regression.py --solver cplex
   ```

6. **Run the pipeline with that environment:**

   ```bash
   nextflow run main.nf -profile apptainer -c cluster.config --solver cplex --resisto_env $HOME/envs/resisto-gem ...
   ```

   `--resisto_env` puts the environment's `bin/` first on `PATH` for the model-building steps only and runs them without a container. Your home folder must be readable from the compute nodes, as it normally is.

   (With `-profile conda`, `--resisto_conda $HOME/envs/resisto-gem` works too, but then the AMRBench_QC steps also use conda instead of their containers.)

Sources for steps 2 and 4: the M4OPT CPLEX installation guide (https://m4opt.readthedocs.io/en/stable/install/cplex.html) and IBM's note that the pip package is the Community Edition (https://supportcontent.ibm.com/support/pages/node/7160884). Older CPLEX releases (≤ 22.1.1) also shipped a `setup.py` in the installation folder; IBM removed it in 22.1.2.

## Checking the installation by hand

```bash
python bin/check_solver.py --solver gurobi     # or cplex
```

```bash
python tests/regression.py --solver gurobi     # or cplex
```

## SCIP

`--solver scip` needs no licence, but:
- CarveMe sets a 600 s limit per MILP for SCIP;
- SCIP accepts suboptimal solutions (the carving MILP was SUBOPTIMAL in testing).

The solver used and every MILP status are recorded per model (`carveme/*.carve.json`, `models/*.build.json`).

## Threads

reframed exposes no thread setting, so Gurobi and CPLEX would otherwise use every core on a node. `bin/resisto_common.py:set_solver` caps them at the CPUs Nextflow allocates to the task. For Gurobi this is tested: a new model gets `Threads` equal to the task CPUs.
