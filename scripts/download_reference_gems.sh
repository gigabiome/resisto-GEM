#!/usr/bin/env bash
# Download curated reference GEMs + their source papers for each species in data/genomes.
# Re-runnable: existing files are skipped. Usage: bash scripts/download_reference_gems.sh
set -uo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$ROOT/reference_gems"
UA="Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0 Safari/537.36"
BIGG="https://bigg.ucsd.edu/static/models"
PMC_S3="https://pmc-oa-opendata.s3.amazonaws.com"

fetch() {  # fetch <url> <dest> [expected `file` type substring]
  local url="$1" dest="$2" want="${3:-}"
  mkdir -p "$(dirname "$dest")"
  if [[ -s "$dest" ]]; then echo "  skip  $(basename "$dest")"; return 0; fi
  if ! curl -sSL --fail -m 600 -A "$UA" "$url" -o "$dest.part"; then
    echo "  FAIL  $(basename "$dest")  <- $url"; rm -f "$dest.part"; return 1
  fi
  if [[ -n "$want" ]] && ! file -b "$dest.part" | grep -qi "$want"; then
    echo "  FAIL  $(basename "$dest")  (got $(file -b "$dest.part" | cut -c1-30), expected $want)"; rm -f "$dest.part"; return 1
  fi
  mv "$dest.part" "$dest"; echo "  ok    $(basename "$dest")"
}

biomodels() {  # biomodels <species_dir> <subdir> <MODEL_ID> <file.omex>  (OMEX = zip, extracted alongside)
  local dest="$OUT/$1/models/$2/${4%.omex}"
  fetch "https://www.ebi.ac.uk/biomodels/model/download/$3?filename=$4" "$dest/$4" Zip \
    && unzip -oq "$dest/$4" -d "$dest"
}

bigg() {  # bigg <species_dir> <model_id>
  fetch "$BIGG/$2.xml.gz" "$OUT/$1/models/$2.xml.gz" gzip && gunzip -kf "$OUT/$1/models/$2.xml.gz"
  fetch "$BIGG/$2.json"   "$OUT/$1/models/$2.json"
}

pmc_pdf() {  # pmc_pdf <species_dir> <PMCID> <filename_stem>
  fetch "$PMC_S3/$2.1/$2.1.pdf" "$OUT/$1/papers/$3.pdf" PDF
}

echo "== Escherichia coli"
bigg Escherichia_coli iML1515
# Ecoli-GEM: SysBioChalmers community-curated continuation of iML1515 (v0.2.0, pinned commit)
ECG="https://raw.githubusercontent.com/SysBioChalmers/Ecoli-GEM/fea539283a3d4bfac5c115f7de93261ae614ef8b"
for f in model/Ecoli-GEM.xml model/Ecoli-GEM.json model/Ecoli-GEM.yml model/Ecoli-GEM.xlsx memote_report.html; do
  fetch "$ECG/$f" "$OUT/Escherichia_coli/models/Ecoli-GEM_v0.2.0/$(basename "$f")"
done
# iML1515 paper is an author manuscript (no OA PDF); PMC only exposes the text.
fetch "$PMC_S3/PMC6521705.1/PMC6521705.1.txt" "$OUT/Escherichia_coli/papers/Monk2017_iML1515_NatBiotechnol.txt"
pmc_pdf Escherichia_coli PMC9393557 Monk2022_Escherichia_222strains_PhilTransB

echo "== Staphylococcus aureus"
bigg Staphylococcus_aureus iYS854
pmc_pdf Staphylococcus_aureus PMC6326480 Seif2019_iYS854_PLoSCompBiol
pmc_pdf Staphylococcus_aureus PMC8241996 Renz2021_Saureus_114GEMs_npjSBA

echo "== Acinetobacter baumannii"
bigg Acinetobacter_baumannii iCN718
pmc_pdf Acinetobacter_baumannii PMC5902709 Norsigian2018_iCN718_FrontGenet
pmc_pdf Acinetobacter_baumannii PMC7696501 Zhao2020_iAB5075_Microorganisms
# Leonidou et al. 2024: new iACB23LX (ATCC 17978) + their curated versions of 7 earlier A. baumannii GEMs
biomodels Acinetobacter_baumannii Leonidou2024 MODEL2309120001 iACB23LX.omex
biomodels Acinetobacter_baumannii Leonidou2024 MODEL2406250005 iLP844.omex
biomodels Acinetobacter_baumannii Leonidou2024 MODEL2406250006 iJS784.omex
biomodels Acinetobacter_baumannii Leonidou2024 MODEL2406250007 iCN718.omex
biomodels Acinetobacter_baumannii Leonidou2024 MODEL2406250008 iATCC1906.omex
biomodels Acinetobacter_baumannii Leonidou2024 MODEL2406250009 iAB5075.omex
biomodels Acinetobacter_baumannii Leonidou2024 MODEL2406250010 AbyMBEL891.omex
biomodels Acinetobacter_baumannii Leonidou2024 MODEL2406250011 AGORA.omex
pmc_pdf Acinetobacter_baumannii PMC11463759 Leonidou2024_iACB23LX_PLoSPathog

echo "== Pseudomonas aeruginosa"
fetch "https://raw.githubusercontent.com/dawsonpayne/iPau21/HEAD/iPau21.xml" \
      "$OUT/Pseudomonas_aeruginosa/models/iPau21.xml" XML
pmc_pdf Pseudomonas_aeruginosa PMC8501023 Payne2021_iPau21_npjSBA
pmc_pdf Pseudomonas_aeruginosa PMC5344303 Bartell2017_iPAE1146_NatCommun
# Dahal et al. 2023: iSD1509 (PA14, BiGG namespace)
biomodels Pseudomonas_aeruginosa . MODEL2205090001 iSD1509.omex
pmc_pdf Pseudomonas_aeruginosa PMC9918512 Dahal2023_iSD1509_CommunBiol

echo "== Klebsiella pneumoniae"
KP="https://raw.githubusercontent.com/kelwyres/KpSC-pan-metabolic-model/HEAD"
for f in KpSC_pan_metabolic_model_v2.xml KpSC_pan_metabolic_model_v2.json \
         KpSC_pan_metabolic_model_v2_nucl.fna KpSC_pan_metabolic_model_v2_prots.faa \
         KpSC_pan_metabolic_model_v2_reaction_metadata.tsv; do
  fetch "$KP/$f" "$OUT/Klebsiella_pneumoniae/models/$f"
done
bigg Klebsiella_pneumoniae iYL1228
pmc_pdf Klebsiella_pneumoniae PMC10926698 Cooper2024_KpSCpan_v2_MicrobGenom
pmc_pdf Klebsiella_pneumoniae PMC9104693  Hawkey2022_Klebsiella_37GEMs_GenomeRes
pmc_pdf Klebsiella_pneumoniae PMC10564454 Vezina2023_Bactabolize_eLife
# Jenior et al. 2022: transcriptome-constrained clinical vs laboratory context-specific models (from BiGG iYL1228)
JEN="https://raw.githubusercontent.com/mjenior/Klebsiella_2021/a198c7a00bdb426824301961fe11e682230d3109/data"
for ctx in clinical_maxfit laboratory_maxfit clinical_maxfit_reps laboratory_maxfit_reps; do
  fetch "$JEN/$ctx/model.sbml"      "$OUT/Klebsiella_pneumoniae/models/Jenior2022_context/$ctx/model.xml" XML
  fetch "$JEN/$ctx/parameters.txt"  "$OUT/Klebsiella_pneumoniae/models/Jenior2022_context/$ctx/parameters.txt"
done
pmc_pdf Klebsiella_pneumoniae PMC9616910 Jenior2022_Kpn_valine_npjSBA

echo "== Enterococcus faecium"
fetch "https://public-pages-files-2025.frontiersin.org/articles/1941652/file/Data_Sheet_2.zip/1941652_data-sheet_2/1" \
      "$OUT/Enterococcus_faecium/models/iDR478_DataSheet2.zip" Zip \
  && unzip -oq "$OUT/Enterococcus_faecium/models/iDR478_DataSheet2.zip" -d "$OUT/Enterococcus_faecium/models/iDR478" -x "__MACOSX/*"
fetch "https://public-pages-files-2025.frontiersin.org/articles/1941652/file/Table_1.xlsx/1941652_table_1/1" \
      "$OUT/Enterococcus_faecium/models/iDR478_Table1.xlsx"
fetch "https://public-pages-files-2025.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2026.1941652/pdf" \
      "$OUT/Enterococcus_faecium/papers/Rasmi2026_iDR478_FrontMicrobiol.pdf" PDF

echo "== Enterobacter spp. (no curated species GEM; closest-relative templates + AGORA2 paper)"
pmc_pdf Enterobacter_spp PMC10497413 Heinken2023_AGORA2_NatBiotechnol

echo "Done -> $OUT"
