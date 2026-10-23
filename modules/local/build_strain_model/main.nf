process BUILD_STRAIN_MODEL {
    tag "${meta.id}"
    label 'process_low'
    label 'resisto'

    input:
    tuple val(meta), path(draft), path(diamond), val(ref), path(biomass_module), path(biomass_json), path(pool), path(pool_sources)
    path mediadb

    output:
    tuple val(meta), path("${meta.id}.xml.gz"), path("${meta.id}.gapfill.tsv"), path("${meta.id}.build.json"), emit: model

    script:
    """
    build_strain_model.py \\
        --sample ${meta.id} \\
        --species "${meta.species}" \\
        --reference-id ${ref.id} \\
        --draft ${draft} \\
        --diamond ${diamond} \\
        --biomass-module ${biomass_module} \\
        --biomass-json ${biomass_json} \\
        --pool ${pool} \\
        --pool-sources ${pool_sources} \\
        --medium ${params.medium} \\
        --mediadb ${mediadb} \\
        --min-growth ${params.min_growth} \\
        --gapfill-max-uptake ${params.gapfill_max_uptake} \\
        --init-max-uptake ${params.init_max_uptake} \\
        --reference-bonus ${params.reference_bonus} \\
        --solver ${params.solver} \\
        --threads ${task.cpus}
    """

    stub:
    """
    touch ${meta.id}.xml.gz ${meta.id}.gapfill.tsv
    echo '{"sample": "${meta.id}", "reference_id": "${ref.id}", "stub": true}' > ${meta.id}.build.json
    """
}
