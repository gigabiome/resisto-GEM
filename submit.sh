#!/bin/bash
#SBATCH --job-name=resisto-gem
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=48:00:00
#SBATCH --output=logs/resisto-gem_%j.out
#SBATCH --error=logs/resisto-gem_%j.err

set -euo pipefail
mkdir -p logs work

module load anaconda3
module load gurobi

source ~/.bashrc
export CONDA_ALWAYS_YES=true

# Activate the conda environment (create it first if needed with: conda env create -f envs/resisto.yml)
conda activate resisto-gem

# Add gurobipy matching cluster Gurobi version (one-time step)
# pip install gurobipy==12.0.3

# Run the pipeline with -profile conda to use the activated environment
nextflow run main.nf \
  --input data/test_sample_sheet.csv \
  --resisto_modules gurobi/12.0.3 \
  -profile conda \
  -work-dir work \
  -c cluster.config
