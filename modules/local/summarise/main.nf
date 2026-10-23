process SUMMARISE {
    label 'process_single'
    label 'resisto'

    input:
    path samplesheet, stageAs: 'samplesheet.csv'
    path qc, stageAs: 'qc/*'
    path build, stageAs: 'build/*'
    path carve, stageAs: 'carve/*'
    path checkm2, stageAs: 'checkm2/*'
    path quast, stageAs: 'quast/*'
    path references, stageAs: 'references/*'

    output:
    path 'resisto_gem_summary.tsv', emit: summary
    path 'reference_summary.tsv', emit: references
    path 'software_versions.json', emit: versions

    script:
    """
    summarise.py --samplesheet samplesheet.csv
    """

    stub:
    """
    touch resisto_gem_summary.tsv reference_summary.tsv
    echo '[]' > software_versions.json
    """
}
