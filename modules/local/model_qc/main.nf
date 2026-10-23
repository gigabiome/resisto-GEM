process MODEL_QC {
    tag "${meta.id}"
    label 'process_low'
    label 'resisto'

    input:
    tuple val(meta), path(model), path(gapfill), path(build_json)

    output:
    tuple val(meta), path("${meta.id}.qc.json"), emit: qc

    script:
    def blocked = params.qc_blocked ? "--blocked --processes ${task.cpus}" : ''
    """
    model_qc.py \\
        --model ${model} \\
        --build-json ${build_json} \\
        --gapfill-tsv ${gapfill} \\
        --sample ${meta.id} \\
        ${blocked}
    """

    stub:
    """
    echo '{"sample": "${meta.id}", "stub": true}' > ${meta.id}.qc.json
    """
}
