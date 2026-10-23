process CARVEME_CARVE {
    tag "${meta.id}"
    label 'process_medium'
    label 'resisto'

    input:
    tuple val(meta), path(proteins, stageAs: 'proteins.faa'), val(ref), path(reference)
    path solver_check

    output:
    tuple val(meta), path("${meta.id}.draft.xml"), path("${meta.id}.diamond.tsv"), emit: draft
    tuple val(meta), path("${meta.id}.carve.json"), emit: status

    script:
    // CarveMe writes the DIAMOND table next to its input (proteins.faa -> proteins.tsv) and runs
    // DIAMOND with all cores unless --threads is given; --diamond-args replaces CarveMe's defaults,
    // so params.diamond_args must contain them (default: CarveMe's own "--more-sensitive --top 10 --quiet").
    """
    carve_wrapper.py \\
        --solver ${params.solver} \\
        --threads ${task.cpus} \\
        --status-json ${meta.id}.carve.json \\
        -- \\
        proteins.faa \\
        -u ${ref.universe} \\
        --reference ${reference} \\
        --reference-score ${params.carve_reference_score} \\
        --diamond-args="${params.diamond_args} --threads ${task.cpus}" \\
        -o ${meta.id}.draft.xml

    mv proteins.tsv ${meta.id}.diamond.tsv
    """

    stub:
    """
    touch ${meta.id}.draft.xml ${meta.id}.diamond.tsv
    echo '{"stub": true}' > ${meta.id}.carve.json
    """
}
