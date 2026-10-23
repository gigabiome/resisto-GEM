# resisto-GEM: method

This document describes what the pipeline does to each genome and why. Each rule is based on behaviour checked in the tool sources (CarveMe 1.6.6, reframed 1.6.0, COBRApy 0.30.0) or measured on the reference models. Section 6 lists those findings.

## 1. Inputs

| Input | Source |
|---|---|
| Assemblies | `sample,species,fasta` samplesheet, e.g. `data/test_sample_sheet.csv` |
| Reference registry | `assets/references.csv`: species → reference GEM, biomass reaction, CarveMe universe, SHA-256 of the reference file |
| Medium | `assets/media.tsv`, medium `LB_btn`: CarveMe 1.6.6's own LB (65 compounds, generated from its `media_db.tsv`) plus biotin |

Species and their references:

| Species | Reference | Biomass reaction | Universe |
|---|---|---|---|
| *E. coli* | iML1515 | BIOMASS_Ec_iML1515_core_75p37M | gramneg |
| *K. pneumoniae* | KpSC-pan v2 | BIOMASS_Core_Feb2022 | gramneg |
| *Enterobacter* spp. | KpSC-pan v2 (no curated *Enterobacter* GEM exists) | BIOMASS_Core_Feb2022 | gramneg |
| *A. baumannii* | iCN718 | BIOMASS__3 | gramneg |
| *P. aeruginosa* | iSD1509 | BIOMASS_PA14_v27M | gramneg |
| *S. aureus* | iYS854 | BIOMASS_iYS_wild_type | grampos |
| *E. faecium* | iDR478 | biomass3 | grampos |

Biotin is in the medium because iDR478 has biotin as a biomass component and cannot grow on CarveMe's LB, which has none. Biotin was the only missing compound. Real LB contains biotin from yeast extract. The medium is the same for all species.

## 2. Per genome: assembly QC and gene prediction

These steps use the AMRBench_QC modules unchanged, from the git submodule `external/AMRBench_QC` pinned at commit `52001eb`:

- QUAST 5.3.0
- CheckM2 1.1.0
- Bakta 1.12.0, run with `--keep-contig-headers` as in AMRBench_QC, with the full database by default

Bakta's predicted proteins (`.faa`) are the input to CarveMe. QC metrics are reported only; no genome is filtered.

## 3. Once per reference (`bin/prepare_reference.py`)

1. **Checksum.** The reference file must match the SHA-256 in the registry.
2. **Formulas and charges.** These are taken from the SBML `fbc:` attributes. reframed reads them from the SBML notes instead. For iSD1509 this changes 211 formulas, 158 of them `;`-joined lists that reframed would silently drop on writing, and 57 charges. The other references are unaffected.
3. **Compartments** are renamed `x` → `C_x` (CarveMe convention). No reaction or metabolite IDs are changed.
4. **Biomass module.** This is the selected biomass reaction plus every "biomass-like" reference reaction that feeds one of its non-currency reactants or drains one of its non-currency products, found recursively. "Biomass-like" means an objective coefficient, `BIOMASS` in the ID, or a non-exchange boundary reaction.
   - Example: iYS854's biomass is built from 7 lumped precursors (`BIOMASS_PROTEIN_1`, …).
   - Example: iCN718's biomass produces `biomass_c`, which only `DM_biomass_c` removes.
   - Alternative biomass reactions are left out, such as iML1515's WT biomass and iSD1509's anaerobic one.
5. **Gap-fill pool.** This is the CarveMe universe without its generic `Growth` reaction, plus reference reactions that are not in the universe. Rules:
   - Where a reaction ID is in both with a different stoichiometry, the universe version is kept. The reference version is added as `R_<id>__<reference>` only if it uses metabolites the universe does not have. Example: iSD1509's `R_LPADSS` uses *P. aeruginosa* lipid X, and without the variant lipid A cannot be made. If only coefficients differ, the reference version is not added, which avoids near-duplicate reactions that could form energy-generating cycles.
   - Sinks are not added. These are non-exchange boundary reactions that can supply a metabolite, including mislabelled exchanges such as iDR478's `R_Ex_lac_L__e` (free lactate) and `R_MEOHE` (free methanol). CarveMe's medium handling only controls `R_EX_*` reactions, so such reactions would be uncontrolled uptakes.
   - Pool reactions carry no gene rules. Reference GPRs name the reference strain's genes, not the new strain's.
6. **Required sinks.** If pool + module cannot grow on the medium, the smallest set of reference sinks that restores growth is found greedily (alphabetical order) and added to the biomass module:
   - iYS854 needs `SK_2fe1s_c` and `SK_apoACP_c`;
   - iDR478 needs `SK_apoACP_c`;
   - the other references need none.

   These supply protein-bound cofactors that the references do not synthesise.
7. **Checks** (the run stops if either fails):
   - the reference grows on the medium with its own biomass;
   - pool + module grows on the medium. A CarveMe draft is a subset of the universe, so this guarantees that gap-filling is feasible for every genome.
8. **Balance diagnostics.** Each reference reaction is reported as balanced, imbalanced or unknown, both in the reference and in the pool. Shared BiGG metabolite IDs can carry different protonation states in the universe, so some reference reactions become imbalanced on universe metabolites: 64 for iML1515, 70 for iSD1509, 20 for iYS854 and 8 for iDR478. These are reported, not altered.

## 4. Per genome: reconstruction

1. **Carving** (`bin/carve_wrapper.py`) runs the unmodified CarveMe command:

   ```
   carve proteins.faa -u <universe> --reference <harmonised reference> --reference-score 0.0 \
         --diamond-args="--more-sensitive --top 10 --quiet --threads N" --solver cplex
   ```

   There is no `-g`: gap-filling is done in step 2.
   - `--reference` only changes the carving scores of reactions that are in both the universe and the reference. With the default `--reference-score 0` it doubles their annotation score and exempts unannotated reference reactions from the −1 penalty.
   - `--diamond-args` replaces CarveMe's DIAMOND defaults, so the same defaults are passed explicitly, plus a thread limit.
   - The wrapper caps CPLEX/Gurobi threads at the CPUs SLURM allocates, records the MILP status, and exits with an error if no model was written. `carve` itself returns 0 even on failure.
2. **Biomass and gap-filling** (`bin/build_strain_model.py`):
   1. CarveMe's `Growth` reaction is removed and the reference biomass module added. The reference biomass becomes the only objective.
   2. `R_ATPM` bounds come from the reference (NGAM): 6.86 for iML1515, 8.39 fixed for iSD1509, 6.8 fixed for KpSC-pan, 0.81 for iDR478, 0 for iYS854 and iCN718.
   3. Gap-filling uses CarveMe's own `gapFill()` MILP, as in `carve -g`, with the pool in place of the universe. The cost of adding reaction *r* is 1/(1 + *s_r*), where:
      - *s_r* is CarveMe's homology score for this genome, from its DIAMOND table and BiGG GPRs via `reaction_scoring()`, with BiGG GPRs filtered to pool reactions;
      - plus 1.0 if *r* is a reference reaction (`--reference_bonus`);
      - unscored reactions cost 1.

      The medium is `LB_btn` with uptake ≤ 100 (CarveMe CLI value), and growth must be ≥ 0.1.
   4. The medium is then stored as the model's default bounds, with uptake ≤ 10 (CarveMe `-i`).
   5. Gap-filled reactions get no GPR, as in CarveMe. The homology-based candidate GPR is listed in `<sample>.gapfill.tsv`.
   6. The model is written as gzip-compressed SBML L3 FBC with BiGG IDs and no `boundaryCondition` species (see section 6).
3. **QC** (`bin/model_qc.py`, COBRApy, an independent parser and LP solver). All checks are reported; none filters:
   - SBML validation;
   - the objective is the reference biomass;
   - growth matches the build step;
   - growth with every carbon-containing uptake closed (free lunch);
   - maximum ATP hydrolysis with all uptakes closed (energy-generating cycle);
   - mass and charge balance, judged only where every formula/charge is set (`fbc:` attributes read with libSBML);
   - blocked reactions (FVA);
   - GPR coverage.

## 5. What the summary contains

`resisto_gem_summary.tsv` has one row per samplesheet entry: CheckM2 and QUAST metrics, the carving and gap-filling MILP status, gap-fill counts by source (reference-only / shared / universe-only), model size, growth and QC flags.

Genomes that failed after retries are kept, with status `no_model` or `failed_after_carve`. `reference_summary.tsv` and `references/*.validation.json` document each prepared reference.

## 6. Findings that shaped the method

| Finding | Consequence |
|---|---|
| `carve --reference` does not change gap-filling or biomass (CarveMe 1.6.6 source) | Biomass replacement and reference gap-filling are separate steps |
| iACB23LX's biomass equals CarveMe's generic Gram-negative `Growth` | iCN718 is used for *A. baumannii* |
| reframed's `bigg` SBML flavour writes extracellular species with `boundaryCondition="true"`; COBRApy then creates `EX_<met>` reactions that collide with the model's own exchanges and open uptake (CarveMe *E. coli* benchmark model: growth 50.2 instead of 3.48) | All models are written without the flag. **Re-saving these models with reframed re-introduces it**; load them with COBRApy only from the files this pipeline writes |
| COBRApy reads a missing `fbc:charge` as 0; CarveMe's universe lacks charges for ~1,000 metabolites | Charge balance is judged only where a charge is set |
| SCIP returns CarveMe's carving MILP as SUBOPTIMAL (E. coli test: SUBOPTIMAL after 197 s); CarveMe accepts this silently | CPLEX is used; MILP status is recorded for every model |
| CarveMe's bundled DIAMOND database was built by DIAMOND build 163 (v2.1.9) | `diamond=2.1.9` is pinned |
| The test E. coli draft gained 9 iML1515-only reactions (SUF Fe-S assembly, KDO2-lipid IV export, undecaprenyl-diphosphate synthase, Ni²⁺ exchange), all with homology support | Confirms the pool adds reactions CarveMe's universe lacks |

## 7. Known limitations

- **No reference GPRs.** Reference reactions without a BiGG GPR in CarveMe's database get no homology score. Example: species-specific reactions of iYS854/iDR478 count as unscored (cost 0.5 with the reference bonus).
- **Alternative optima.** Gap-filling solutions can differ between solvers when costs tie.
- **Growth rates aren't predictions.** Growth rates on `LB_btn` are model properties, not predictions; CarveMe's LB is a defined approximation (e.g. it contains glucose).
- **Enterobacter.** It uses the *Klebsiella* pan-model; its biomass composition is that of KpSC.
