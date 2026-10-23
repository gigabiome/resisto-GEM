//
// Read and validate the samplesheet (sample,species,fasta) and the reference registry.
// Stops the run on: missing/duplicate sample names, spaces in names, missing files,
// species without a registry entry, or a registry reference file that does not exist.
//

workflow INPUT_CHECK {
    take:
    samplesheet // path: csv with columns sample,species,fasta
    registry    // path: assets/references.csv
    reference_dir

    main:
    def refs = [:]
    file(registry, checkIfExists: true).splitCsv(header: true, strip: true).each { row ->
        if (refs.containsKey(row.species)) {
            error("Reference registry: species '${row.species}' appears more than once")
        }
        def ref_file = file("${reference_dir}/${row.reference_file}")
        if (!ref_file.exists()) {
            error("Reference registry: ${ref_file} not found (species '${row.species}'). Run scripts/download_reference_gems.sh or set --reference_dir.")
        }
        refs[row.species] = [
            id      : row.reference_id,
            biomass : row.biomass,
            universe: row.universe,
            sha256  : row.sha256,
            file    : ref_file,
        ]
    }

    def seen = [] as Set
    def rows = file(samplesheet, checkIfExists: true).splitCsv(header: true, strip: true)
    rows.each { row ->
        if (!row.sample) {
            error("Samplesheet: empty 'sample' in row ${row}")
        }
        if (row.sample.contains(' ')) {
            error("Samplesheet: sample name '${row.sample}' contains spaces")
        }
        if (!seen.add(row.sample)) {
            error("Samplesheet: sample '${row.sample}' appears more than once")
        }
        if (!refs.containsKey(row.species)) {
            error("Samplesheet: species '${row.species}' (sample ${row.sample}) has no entry in ${registry}")
        }
        if (!row.fasta) {
            error("Samplesheet: sample '${row.sample}' has no 'fasta'")
        }
        if (!file(row.fasta).exists()) {
            error("Samplesheet: file '${row.fasta}' for sample '${row.sample}' does not exist (relative paths are resolved from the launch directory)")
        }
    }

    assemblies = channel.fromList(rows).map { row ->
        def meta = [id: row.sample, species: row.species, reference_id: refs[row.species].id]
        [meta, file(row.fasta, checkIfExists: true)]
    }

    // One entry per reference model actually used by this samplesheet
    def used = rows.collect { row -> refs[row.species] }.unique { r -> r.id }
    references = channel.fromList(used).map { r ->
        [[id: r.id, biomass: r.biomass, universe: r.universe, sha256: r.sha256], r.file]
    }

    emit:
    assemblies // channel: [ val(meta), fasta ]
    references // channel: [ val(ref), reference_sbml ]
}
