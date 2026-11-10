#!/usr/bin/env python3
"""
Organize bacterial genomes by species using BV-BRC metadata.

This script reads metadata from a parquet file and organizes genome FASTA files
into subdirectories named after their species.

Usage:
    python3 organize_by_species.py [--data-dir DATA_DIR] [--metadata FILE] [--dry-run]

Expected directory structure:
    data/
    ├── genomes/           (contains all *.fna files to organize)
    └── bv_brc_metadata.parquet
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict
import shutil

try:
    import pandas as pd
except ImportError:
    print("ERROR: pandas not found. Install with: pip install pandas pyarrow")
    sys.exit(1)


def sanitize_species_name(species: str) -> str:
    """Convert species name to valid directory name."""
    return species.replace(' ', '_').replace('.', '').replace('(', '').replace(')', '')


def organize_genomes(data_dir: Path, metadata_file: Path, dry_run: bool = False) -> None:
    """
    Organize genomes into species directories.

    Args:
        data_dir: Path to data directory containing genomes/ subdirectory
        metadata_file: Path to bv_brc_metadata.parquet
        dry_run: If True, print operations without executing
    """

    genomes_dir = data_dir / 'genomes'

    # Validate paths
    if not genomes_dir.exists():
        print(f"ERROR: Genomes directory not found: {genomes_dir}")
        sys.exit(1)

    if not metadata_file.exists():
        print(f"ERROR: Metadata file not found: {metadata_file}")
        sys.exit(1)

    # Load metadata
    print(f"Loading metadata from {metadata_file}")
    try:
        df_metadata = pd.read_parquet(metadata_file)
    except Exception as e:
        print(f"ERROR: Failed to read metadata: {e}")
        sys.exit(1)

    # Get species column
    species_col = None
    for col in ['species_collapsed', 'Species', 'species']:
        if col in df_metadata.columns:
            species_col = col
            break

    if not species_col:
        print(f"ERROR: No species column found. Available: {df_metadata.columns.tolist()}")
        sys.exit(1)

    print(f"Using species column: {species_col}")

    # Create biosample -> species mapping (unique)
    df_unique = df_metadata[['biosample', species_col]].drop_duplicates(subset=['biosample'])
    species_map = dict(zip(df_unique['biosample'].astype(str), df_unique[species_col].astype(str)))

    print(f"Loaded {len(species_map)} unique biosamples from metadata\n")

    # Find all genome files
    genome_files = sorted(genomes_dir.glob('*.fna'))

    if not genome_files:
        print(f"WARNING: No *.fna files found in {genomes_dir}")
        sys.exit(1)

    print(f"Found {len(genome_files)} genome files to organize\n")

    # Group by species
    species_groups = defaultdict(list)
    unmapped = []

    for genome_file in genome_files:
        biosample = genome_file.stem

        if biosample in species_map:
            species = species_map[biosample]
            species_groups[species].append(genome_file)
        else:
            unmapped.append(biosample)

    # Print summary
    print(f"Species distribution:")
    for species in sorted(species_groups.keys()):
        count = len(species_groups[species])
        pct = 100.0 * count / len(genome_files)
        print(f"  {species:45s}: {count:6d} ({pct:5.1f}%)")

    if unmapped:
        print(f"  {'UNMAPPED':45s}: {len(unmapped):6d} ({100.0*len(unmapped)/len(genome_files):5.1f}%)")

    print()

    # Organize files
    if dry_run:
        print("[DRY RUN MODE - No changes will be made]\n")

    total_moved = 0
    total_existed = 0
    total_unmapped = len(unmapped)

    for species in sorted(species_groups.keys()):
        species_dir_name = sanitize_species_name(species)
        species_dest = genomes_dir / species_dir_name

        files_in_species = species_groups[species]

        if not dry_run:
            species_dest.mkdir(parents=True, exist_ok=True)

        for genome_file in files_in_species:
            dest_file = species_dest / genome_file.name

            if dest_file.exists():
                total_existed += 1
            else:
                if not dry_run:
                    try:
                        shutil.move(str(genome_file), str(dest_file))
                        total_moved += 1
                    except Exception as e:
                        print(f"  ERROR moving {genome_file.name}: {e}")
                else:
                    print(f"  Would move: {genome_file.name} -> {species_dir_name}/")
                    total_moved += 1

    # Print results
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    print(f"Total genomes processed:    {len(genome_files)}")
    print(f"  ✓ Moved to species dirs:  {total_moved}")
    print(f"  → Already in place:       {total_existed}")
    print(f"  ✗ Unmapped:               {total_unmapped}")
    print(f"Total species groups:       {len(species_groups)}")

    if dry_run:
        print("\n[DRY RUN - No changes were made]")
    else:
        print("\n✓ Organization complete!")

    # Create sample sheet
    sample_sheet_path = genomes_dir / 'sample_sheet_by_species.csv'

    sample_rows = []
    for species in sorted(species_groups.keys()):
        species_dir_name = sanitize_species_name(species)
        for genome_file in species_groups[species]:
            sample_rows.append({
                'sample': genome_file.stem,
                'species': species,
                'fasta': f"data/genomes/{species_dir_name}/{genome_file.name}"
            })

    if sample_rows:
        df_samples = pd.DataFrame(sample_rows)
        if not dry_run:
            df_samples.to_csv(sample_sheet_path, index=False)
            print(f"✓ Sample sheet created: {sample_sheet_path}")
        else:
            print(f"[DRY RUN] Would create sample sheet: {sample_sheet_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Organize bacterial genomes by species using metadata',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Organize genomes (assume data/ directory in current path)
  python3 organize_by_species.py

  # Specify custom data directory
  python3 organize_by_species.py --data-dir /path/to/data

  # Dry run to see what would happen
  python3 organize_by_species.py --dry-run
        '''
    )

    parser.add_argument(
        '--data-dir',
        type=Path,
        default=Path('data'),
        help='Path to data directory containing genomes/ and metadata (default: data/)'
    )

    parser.add_argument(
        '--metadata',
        type=Path,
        default=None,
        help='Path to metadata parquet file (default: data/bv_brc_metadata.parquet)'
    )

    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Print operations without executing'
    )

    args = parser.parse_args()

    data_dir = args.data_dir.resolve()

    if args.metadata is None:
        metadata_file = data_dir / 'bv_brc_metadata.parquet'
    else:
        metadata_file = args.metadata.resolve()

    print(f"Data directory:  {data_dir}")
    print(f"Metadata file:   {metadata_file}\n")

    organize_genomes(data_dir, metadata_file, dry_run=args.dry_run)


if __name__ == '__main__':
    main()
