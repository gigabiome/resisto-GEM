process PREPARE_REFERENCE {
    tag "${ref.id}"
    label 'process_single'
    label 'resisto'

    input:
    tuple val(ref), path(reference)
    path mediadb
    path solver_check // only a dependency: no reference is prepared before the solver check passes

    output:
    tuple val(ref), path("${ref.id}.harmonised.xml"), emit: harmonised
    tuple val(ref), path("${ref.id}.biomass_module.xml"), path("${ref.id}.biomass.json"), path("${ref.id}.pool.xml"), path("${ref.id}.pool_sources.tsv"), emit: bundle
    path "${ref.id}.validation.json", emit: validation

    script:
    """
    prepare_reference.py \\
        --reference ${reference} \\
        --reference-id ${ref.id} \\
        --sha256 ${ref.sha256} \\
        --biomass ${ref.biomass} \\
        --universe ${ref.universe} \\
        --medium ${params.medium} \\
        --mediadb ${mediadb} \\
        --max-uptake ${params.gapfill_max_uptake} \\
        --atpm-source ${params.atpm_source} \\
        --min-growth ${params.min_growth} \\
        --solver ${params.solver} \\
        --threads ${task.cpus}
    """

    stub:
    """
    touch ${ref.id}.harmonised.xml ${ref.id}.biomass_module.xml ${ref.id}.pool.xml ${ref.id}.pool_sources.tsv
    echo '{}' > ${ref.id}.biomass.json
    echo '{"reference_id": "${ref.id}", "stub": true}' > ${ref.id}.validation.json
    """
}
