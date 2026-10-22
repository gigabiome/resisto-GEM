# Reference GEMs

Curated genome-scale metabolic models (GEMs) and their source papers. These are the references for building strain-specific models from the 140 assemblies in `../data/genomes` (20 per species). The full catalogue, with paths and DOIs, is in `references.tsv`. To re-download, run `bash ../scripts/download_reference_gems.sh`. Files that are already present are skipped.

| Species | Primary reference | Strain | Namespace | Rxn / Met / Genes | Grows out of the box? |
|---|---|---|---|---|---|
| *E. coli* | iML1515 (BiGG) | K-12 MG1655 | BiGG | 2712 / 1877 / 1516 | yes (0.88 h⁻¹) |
| *E. coli* (alt.) | Ecoli-GEM v0.2.0 (SysBioChalmers) | K-12 MG1655 | BiGG | 2710 / 1877 / 1516 | original: **no**; `Ecoli-GEM_HSDyfix.xml`: yes (0.88 h⁻¹) |
| *S. aureus* | iYS854 (BiGG) | USA300 | BiGG | 1455 / 1335 / 866 | yes |
| *A. baumannii* | **iACB23LX** (Leonidou 2024, BioModels) | ATCC 17978 | BiGG | 2321 / 1660 / 1164 | yes (0.55 h⁻¹) |
| *A. baumannii* (alt.) | iCN718 (BiGG), plus 7 re-curated models from Leonidou 2024 | AYE and others | BiGG / mixed | 1015 / 888 / 709 | yes |
| *P. aeruginosa* | **iSD1509** (Dahal 2023, BioModels) | PA14 | BiGG | 2021 / 1642 / 1510 | yes (1.59 h⁻¹, LB medium by default) |
| *P. aeruginosa* (alt.) | iPau21 (GitHub) | PA14 | **ModelSEED** | 1588 / 1325 / 1169 | **no**: all exchanges are closed, so a medium must be set first |
| *K. pneumoniae* | KpSC-pan v2 (GitHub) and iYL1228 | pan-model of 507 strains / MGH 78578 | BiGG | 3547 / 2577 / 2400 | yes |
| *E. faecium* | iDR478 (Frontiers supplement) | TX0016 | BiGG | 938 / 638 / 480 | yes |
| *Enterobacter* spp. | none exists (see below) | — | — | — | — |

All models were loaded with COBRApy 0.30 in the `gem` conda env and optimised.

## Notes for the pipeline

- **Enterobacter has no curated GEM.** Literature searches up to Oct 2026 found none. The best options are:
  - use the BiGG-namespace Enterobacterales models as gap-fill templates: iML1515 + KpSC-pan v2, both closely related to *Enterobacter*;
  - or pull *E. cloacae* models from AGORA2 (vmh.life). These are semi-curated and use the VMH namespace, which would need ID mapping. The paper is included here, but the models were not downloaded, because vmh.life does not serve direct file links.
- **Klebsiella (Jenior 2022)** is in `models/Jenior2022_context/`. These are four transcriptome-constrained, context-specific models built from BiGG iYL1228: clinical vs laboratory isolates, each in a single-sample and a replicate (`_reps`) version. They are pruned to about 300–390 reactions, so they are **not suitable as gap-filling references**. Use them to compare clinical and laboratory metabolic states. Jenior et al. used iYL1228 unmodified, and that model is already in `models/`.
- **Klebsiella:** KpSC-pan v2 ships with protein and nucleotide sequences for its GPRs, so it plugs straight into [Bactabolize](https://github.com/kelwyres/Bactabolize). Bactabolize is a reference-based strain-specific GEM builder made for exactly this task. The same approach (BLAST/DIAMOND the reference GPR genes against each assembly, then gap-fill) carries over to the other species.
- **A. baumannii (Leonidou 2024)** is in `models/Leonidou2024/`. It contains the new iACB23LX (MEMOTE score 89%), plus the authors' re-curated versions of iCN718, iAB5075, iATCC1906 (*sic*: this is iATCC19606), iLP844, iJS784, AbyMBEL891 and AGORA AB0057. Each one is an unpacked OMEX archive with a FROG report. Two caveats:
  - Every re-curated model except iACB23LX ships with a computed *minimal medium* that caps growth at 0.1 h⁻¹. They grow normally once exchanges are opened.
  - Two models are not usable: **iJS784** cannot grow even with all exchanges open (the paper reports this too), and **AbyMBEL891** has no GPRs and uses legacy `R###` IDs.
  - iLP844 and iJS784 use ModelSEED IDs, and iATCC1906 uses KEGG IDs.
- **P. aeruginosa:** iSD1509 is now the primary reference. It is newer and larger, it grows as supplied, and it uses BiGG IDs like the other references, so no namespace mapping is needed. iPau21 uses ModelSEED IDs (`rxn…`, `cpd…`). Every other reference uses BiGG IDs, so map namespaces (e.g. via MetaNetX) before doing any cross-species comparison.
- **Ecoli-GEM** ([SysBioChalmers/Ecoli-GEM](https://github.com/SysBioChalmers/Ecoli-GEM), v0.2.0, pinned to commit `fea5392`) is a community-curated continuation of iML1515. It drops GLYCK and XYLI2, changes the direction of GHMT2r, HSDy and SUCCt1pp, corrects the THD2pp stoichiometry and changes 43 GPRs, following Mao 2022, Bernstein 2023 and Carrao 2024. As released, it **cannot grow**: HSDy is irreversible in the oxidative direction, which blocks Thr/Met/Ile synthesis (upstream issue #34; fix in PR #36, still open). `../scripts/patch_ecoli_gem.py` applies the upstream fix (HSDy runs aspsa → hom__L only) and writes `Ecoli-GEM_HSDyfix.xml/json`; use that file. The repo's MEMOTE report is saved alongside.
- **E. coli:** the iML1515 paper is an author manuscript and has no open-access PDF. Only the PMC text is saved here. The 222-strain *Escherichia* paper (Monk 2022) provides strain tables, not model files.
- **Download manually** (paywalled or bot-protected): Liao et al. 2011, the iYL1228 paper (doi:10.1128/JB.01218-10), and the iML1515 PDF (doi:10.1038/nbt.3956).
