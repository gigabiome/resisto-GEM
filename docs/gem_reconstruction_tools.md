# Automated GEM reconstruction tools: reference list

Compiled 2026-10-22 as input for choosing the strain-specific GEM pipeline. The pipeline must cover 140 genomes now and scale to about 12,000 later.

**How this was verified.** Every DOI below was resolved against Crossref on 2026-10-22, so titles, first authors, years and venues are confirmed. Repository status (latest release, last push, licence) comes from the GitHub API on the same date. Runtime figures are quoted from the cited papers and attributed to them. They are **not** our own measurements, and they come from other hardware and other genomes.

## A. De novo reconstruction from a genome, for any species

| Tool | Approach | Namespace | Language / solver | Latest release (last push) | Licence | Primary reference |
|---|---|---|---|---|---|---|
| **CarveMe** | Top-down: carves a manually curated universal model, using DIAMOND homology scores | BiGG | Python. MILP solver: CPLEX or Gurobi (full academic licence), or SCIP (free; docs say ≥10 min per run) | 1.6.5, 2025-08-29 (tag 1.6.6; push 2025-09-12) | Apache-2.0 | Machado D *et al.* (2018) *Nucleic Acids Res* 46:7542. doi:10.1093/nar/gky537 |
| **gapseq** | Bottom-up: pathway prediction from a curated reference sequence DB, then informed gap-filling | ModelSEED | R + shell. Solver: GLPK or CPLEX | v2.1.0, 2026-05-30 (push 2026-10-20). Reference DB v1.5 (UniProt 2026_01) | GPL-3.0 | Zimmermann J *et al.* (2021) *Genome Biol* 22:81. doi:10.1186/s13059-021-02295-1 |
| **ModelSEED / ModelSEEDpy** (also run in KBase) | Template-based (RAST annotation → template → gap-fill) | ModelSEED | Python | no tagged releases (push 2026-02-11) | MIT | Henry CS *et al.* (2010) *Nat Biotechnol* 28:977. doi:10.1038/nbt.1672 · Seaver SMD *et al.* (2021) *Nucleic Acids Res* (ModelSEED biochemistry DB). doi:10.1093/nar/gkaa746 · Faria JP *et al.* (2023) ModelSEED v2, **bioRxiv preprint**. doi:10.1101/2023.10.04.556561 · KBase: Arkin AP *et al.* (2018) *Nat Biotechnol*. doi:10.1038/nbt.4163 |
| **Reconstructor** | KEGG-annotation-based draft with parsimonious gap-filling; outputs COBRApy-native models | ModelSEED | Python | v1.2.0, 2025-10-25 | MIT | Jenior ML *et al.* (2023) *Bioinformatics* 39:btad367. doi:10.1093/bioinformatics/btad367 |
| **RAVEN 2** | Homology to template models and/or de novo from KEGG/MetaCyc | KEGG / MetaCyc / template | MATLAB | v2.11.4, 2026-10-07 | GPL-3.0 | Wang H *et al.* (2018) *PLoS Comput Biol* 14:e1006541. doi:10.1371/journal.pcbi.1006541 |
| **merlin** | Interactive (GUI) annotation and reconstruction workbench | not verified | not verified | not checked | not checked | Capela J *et al.* (2022) *Nucleic Acids Res* 50:6052. doi:10.1093/nar/gkac459 |
| **Pathway Tools (PathoLogic + MetaFlux)** | Pathway/genome database → FBA model | MetaCyc | not verified (licence required; terms not checked) | not checked | not checked | Karp PD *et al.* (2021) *Brief Bioinform* 22:109. doi:10.1093/bib/bbz104 |
| **AuReMe / Pantograph** | Traceable "à-la-carte" reconstruction, combining template and annotation sources | MetaCyc / mixed | Python | no releases (push 2024-04-24) | GPL-3.0 | Aite M *et al.* (2018) *PLoS Comput Biol* 14:e1006146. doi:10.1371/journal.pcbi.1006146 |
| **AGORA / AGORA2 (DEMETER pipeline)** | Semi-automated refinement pipeline used to build 7,302 gut-microbe models | VMH | MATLAB (COBRA Toolbox) | COBRA Toolbox v3.7.1, 2026-07-13 | — | Magnúsdóttir S *et al.* (2017) *Nat Biotechnol* 35:81. doi:10.1038/nbt.3703 · Heinken A *et al.* (2023) *Nat Biotechnol* 41:1320. doi:10.1038/s41587-022-01628-0 |

## B. Reference-based, strain-specific (curated species model → strain models)

These match our setup, since a curated reference GEM exists for 6 of our 7 species.

| Tool | Approach | Namespace | Latest release | Licence | Primary reference |
|---|---|---|---|---|---|
| **Bactabolize** | Removes reference genes/reactions absent from each assembly (BLAST against the reference's GPR sequences), then gap-fills; predicts growth phenotypes | inherits the reference (KpSC-pan → BiGG) | 1.0.5, 2025-05-23 (tag 1.0.6) | GPL-3.0 | Vezina B *et al.* (2023) *eLife* 12:RP87406. doi:10.7554/eLife.87406 |
| **Multi-strain workflow (Norsigian *et al.*)** | Protocol: bidirectional BLAST against a curated reference, then remove absent genes/reactions, then manual QC. A protocol, not a packaged tool | inherits the reference | — | — | Norsigian CJ *et al.* (2020) *Nat Protoc* 15:1. doi:10.1038/s41596-019-0254-3 |
| **gempipe** | Drafts a pan-GEM (using CarveMe or a curated reference), then derives strain-specific models and supports curation | BiGG (via CarveMe) | v1.38.5, 2025-12-15 | GPL-3.0 | Lazzari G *et al.* (2026) *mSystems*. doi:10.1128/msystems.01007-25 (preprint doi:10.1101/2025.07.03.662949) |
| **CarveMe, reference mode** | `carve --reference` propagates a curated model to new strains instead of the universal model | inherits the reference (must be BiGG) | as CarveMe | Apache-2.0 | Machado 2018 (above). Benchmarked against Bactabolize in Vezina 2023 |

## C. Layers on top of the builders (pan / consensus / workflow)

| Tool | What it adds | Reference |
|---|---|---|
| **pan-Draft** (part of gapseq: `gapseq pan`) | Species-representative pan-models built from many (incomplete) genomes | De Bernardini N *et al.* (2024) *Genome Biol* 25:280. doi:10.1186/s13059-024-03425-1 |
| **GEMsembler** | Consensus models assembled from several tools' outputs (CarveMe, gapseq, ModelSEED, AGORA), plus structural comparison | Matveishina E *et al.* (2025) *mSystems*. doi:10.1128/msystems.00574-25 (repo tag v0.1.0) |
| **metaGEM** | Snakemake workflow wrapping CarveMe for MAG-scale reconstruction on HPC | Zorrilla F *et al.* (2021) *Nucleic Acids Res* 49:e126. doi:10.1093/nar/gkab815 (last release v1.0.5, 2023) |
| **CarveMe-GutMicrobes** | CarveMe variant with a gut-microbe-specific universal model | Basile A *et al.* (2026) **bioRxiv preprint**. doi:10.64898/2026.06.26.734454 |

## D. Independent assessments and critiques, worth reading before choosing

- Mendoza SN *et al.* (2019) A systematic assessment of current genome-scale metabolic reconstruction tools. *Genome Biol* 20:158. doi:10.1186/s13059-019-1769-1
- Price MN (2023; online Dec 2022) Erroneous predictions of auxotrophies by CarveMe. *Nat Ecol Evol* 7:194. doi:10.1038/s41559-022-01936-3. Authors' reply: Machado D & Patil KR (2023). doi:10.1038/s41559-022-01939-0
- The gapseq paper (Zimmermann 2021) benchmarks gapseq against CarveMe and ModelSEED. It is **self-reported**: those authors built gapseq.
- The Bactabolize paper (Vezina 2023) benchmarks Bactabolize, CarveMe (universal and reference mode), gapseq and KBase on *Klebsiella*. That is the closest match to our use case, but it is also partly self-reported.

## Reported runtimes (relevant to scaling to 12,000 genomes)

These come from Vezina 2023 (10 KpSC genomes, Intel Xeon Gold 6150 HPC node), timing draft-model construction only. The last column is a naive linear projection to 12,000 genomes, for orientation only.

| Tool | Time per genome (reported) | ≈ CPU-time for 12,000 genomes |
|---|---|---|
| CarveMe, KpSC-pan reference mode | 20.0 s | ~67 h |
| CarveMe, universal | 30.3 s | ~100 h |
| Bactabolize, KpSC-pan | 98.1 s | ~330 h |
| KBase (batch, includes upload and queue) | 183.5 s | ~610 h (web platform, not a local HPC run) |
| gapseq (draft only, **excludes** gap-filling) | 5.46 h | ~65,500 h |

The pan-Draft paper reports that its own step takes "a few minutes", run on top of existing gapseq drafts. No comparable timing was found for ModelSEEDpy, Reconstructor, RAVEN, gempipe or merlin.

## Properties relevant to our design (facts, not a recommendation)

- **Namespace.** The curated references we downloaded are all BiGG (iML1515, Ecoli-GEM, iYS854, iACB23LX, iCN718, iSD1509, KpSC-pan, iYL1228, iDR478), except iPau21, which uses ModelSEED. BiGG-native tools (CarveMe, Bactabolize, gempipe) can use these directly as references. gapseq, ModelSEED and Reconstructor output ModelSEED IDs, so comparing their models with our references needs ID mapping (e.g. MetaNetX).
- **Solver licensing.** CarveMe at full speed needs CPLEX or Gurobi. Both are free for academics but cannot be bundled into a container. The free SCIP option is documented as much slower. gapseq runs with GLPK, which is free.
- **Species without a curated reference.** *Enterobacter* in our set has none. Whether the 12,000-genome set contains further species without references determines whether a purely reference-based tool is sufficient.
