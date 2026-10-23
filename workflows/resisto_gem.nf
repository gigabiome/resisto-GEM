/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { INPUT_CHECK        } from '../subworkflows/local/input_check'

// Gene prediction and assembly QC: AMRBench_QC modules, used unchanged from the pinned submodule
include { GUNZIP             } from '../external/AMRBench_QC/modules/local/gunzip/main'
include { QUAST              } from '../external/AMRBench_QC/modules/local/quast/main'
include { CHECKM2_DBDOWNLOAD } from '../external/AMRBench_QC/modules/local/checkm2_dbdownload/main'
include { CHECKM2            } from '../external/AMRBench_QC/modules/local/checkm2/main'
include { BAKTA_DBDOWNLOAD   } from '../external/AMRBench_QC/modules/local/bakta_dbdownload/main'
include { BAKTA              } from '../external/AMRBench_QC/modules/local/bakta/main'

// Model reconstruction
include { CHECK_SOLVER       } from '../modules/local/check_solver/main'
include { PREPARE_REFERENCE  } from '../modules/local/prepare_reference/main'
include { CARVEME_CARVE      } from '../modules/local/carveme_carve/main'
include { BUILD_STRAIN_MODEL } from '../modules/local/build_strain_model/main'
include { MODEL_QC           } from '../modules/local/model_qc/main'
include { SUMMARISE          } from '../modules/local/summarise/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RESISTO_GEM {
    main:
    if (!params.input) {
        error("Please provide a samplesheet, e.g. --input data/test_sample_sheet.csv")
    }
    if (!(params.bakta_db_type in ['full', 'light'])) {
        error("--bakta_db_type must be 'full' or 'light', not '${params.bakta_db_type}'")
    }
    if (!(params.atpm_source in ['reference', 'carveme'])) {
        error("--atpm_source must be 'reference' or 'carveme'")
    }

    INPUT_CHECK(params.input, params.reference_registry, params.reference_dir)
    mediadb = channel.value(file(params.mediadb, checkIfExists: true))

    //
    // Assemblies: decompress gzipped FASTA (AMRBench_QC GUNZIP), then QC and gene prediction
    //
    supplied = INPUT_CHECK.out.assemblies.branch { _meta, fasta ->
        compressed: fasta.name.endsWith('.gz')
        plain: true
    }
    GUNZIP(supplied.compressed)
    assemblies = supplied.plain.mix(GUNZIP.out.fasta)

    QUAST(assemblies)

    if (params.checkm2_db) {
        checkm2_db = channel.value(file(params.checkm2_db, checkIfExists: true))
    } else {
        CHECKM2_DBDOWNLOAD()
        checkm2_db = CHECKM2_DBDOWNLOAD.out.db.collect()
    }
    CHECKM2(assemblies, checkm2_db)

    if (params.bakta_db) {
        bakta_db = channel.value(file(params.bakta_db, checkIfExists: true))
    } else {
        BAKTA_DBDOWNLOAD(params.bakta_db_type)
        bakta_db = BAKTA_DBDOWNLOAD.out.db.collect()
    }
    BAKTA(assemblies, bakta_db)

    //
    // Solver check, then one preparation per reference model used in this run
    //
    CHECK_SOLVER(params.solver)
    solver_ok = CHECK_SOLVER.out.report.collect()

    PREPARE_REFERENCE(INPUT_CHECK.out.references, mediadb, solver_ok)

    //
    // Per genome: CarveMe draft with the species reference, then biomass + reference gap-filling
    //
    harmonised = PREPARE_REFERENCE.out.harmonised.map { ref, xml -> [ref.id, ref, xml] }
    carve_input = BAKTA.out.faa
        .map { meta, faa -> [meta.reference_id, meta, faa] }
        .combine(harmonised, by: 0)
        .map { _key, meta, faa, ref, xml -> [meta, faa, ref, xml] }
    CARVEME_CARVE(carve_input, solver_ok)

    bundle = PREPARE_REFERENCE.out.bundle.map { ref, module, bjson, pool, sources -> [ref.id, ref, module, bjson, pool, sources] }
    build_input = CARVEME_CARVE.out.draft
        .map { meta, draft, diamond -> [meta.reference_id, meta, draft, diamond] }
        .combine(bundle, by: 0)
        .map { _key, meta, draft, diamond, ref, module, bjson, pool, sources -> [meta, draft, diamond, ref, module, bjson, pool, sources] }
    BUILD_STRAIN_MODEL(build_input, mediadb)

    MODEL_QC(BUILD_STRAIN_MODEL.out.model)

    //
    // One table for the whole run (samples without a model are kept, with empty model columns)
    //
    SUMMARISE(
        file(params.input),
        MODEL_QC.out.qc.map { _meta, f -> f }.collect().ifEmpty([]),
        BUILD_STRAIN_MODEL.out.model.map { _meta, _m, _g, b -> b }.collect().ifEmpty([]),
        CARVEME_CARVE.out.status.map { _meta, f -> f }.collect().ifEmpty([]),
        CHECKM2.out.results.map { _meta, d -> d }.collect().ifEmpty([]),
        QUAST.out.report.map { _meta, d -> d }.collect().ifEmpty([]),
        PREPARE_REFERENCE.out.validation.collect(),
    )

    emit:
    summary = SUMMARISE.out.summary
}
