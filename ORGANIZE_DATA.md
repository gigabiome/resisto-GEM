# Organizing BV-BRC Data by Species

This guide explains how to organize the 12,993 BV-BRC genomes by species on the cluster.

## Setup on Cluster

1. **Create the data directory structure:**
   ```bash
   mkdir -p data/genomes
   ```

2. **Copy genome files:**
   ```bash
   # Copy all *.fna files from BV_BRC_data to data/genomes/
   cp ../AMRBench/BV_BRC_data/*.fna data/genomes/
   ```

3. **Copy metadata:**
   ```bash
   cp ../AMRBench/bv_brc_metadata.parquet data/
   ```

## Running the Organization Script

### Prerequisites
Ensure pandas and pyarrow are available (usually in conda or pip):
```bash
pip install pandas pyarrow
# or
conda install -c conda-forge pandas pyarrow
```

### Basic Usage
From the resisto-GEM directory on the cluster:

```bash
python3 bin/organize_by_species.py
```

This will:
- Read `data/bv_brc_metadata.parquet`
- Organize all `*.fna` files in `data/genomes/` into species-specific subdirectories
- Create `data/genomes/sample_sheet_by_species.csv` with metadata for all genomes

### Options

**Dry run** (preview what will happen without making changes):
```bash
python3 bin/organize_by_species.py --dry-run
```

**Custom data directory:**
```bash
python3 bin/organize_by_species.py --data-dir /path/to/data
```

**Custom metadata location:**
```bash
python3 bin/organize_by_species.py --metadata /path/to/metadata.parquet
```

## Expected Output

The script will create directories like:
```
data/genomes/
├── Escherichia_coli/              (108,330 genomes)
├── Klebsiella_pneumoniae/         (32,313 genomes)
├── Enterobacter_spp/              (18,627 genomes)
├── Acinetobacter_baumannii/       (10,690 genomes)
├── Pseudomonas_aeruginosa/        (4,527 genomes)
├── Staphylococcus_aureus/         (3,943 genomes)
├── Enterococcus_faecium/          (2,678 genomes)
├── Staphylococcus_epidermidis/    (1,498 genomes)
└── sample_sheet_by_species.csv
```

The sample sheet contains three columns:
- `sample`: biosample ID (e.g., "562.141493")
- `species`: species name (e.g., "Escherichia coli")
- `fasta`: relative path to the FASTA file

## Using the Organized Data

After organization, you can run the Nextflow pipeline using the sample sheet:

```bash
nextflow run main.nf \
  --input data/genomes/sample_sheet_by_species.csv \
  -c create.config \
  --resisto_env /path/to/conda/envs/resisto-gem \
  --resisto_modules gurobi/9.5.2
```

Or process specific species:

```bash
# Extract just E. coli genomes
awk -F',' '$2 == "Escherichia coli" {print $0}' data/genomes/sample_sheet_by_species.csv > ecoli_samples.csv

# Run pipeline on E. coli
nextflow run main.nf --input ecoli_samples.csv -c create.config ...
```

## Troubleshooting

**"ERROR: pandas not found"**
- Install pandas: `pip install pandas pyarrow`
- Or activate conda environment: `conda activate resisto-gem` then run script

**"ERROR: Metadata file not found"**
- Ensure `data/bv_brc_metadata.parquet` exists
- Check path with: `ls -lh data/bv_brc_metadata.parquet`

**"WARNING: No *.fna files found"**
- Ensure genomes are copied to `data/genomes/`
- Check with: `ls data/genomes/*.fna | wc -l`

**Script is very slow**
- This is normal for 12,993 files. The first run may take 10-30 minutes
- Use `--dry-run` first to verify without moving files
