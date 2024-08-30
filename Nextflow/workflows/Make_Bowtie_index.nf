
process make_bowtie_index {
    label 'process_high'
    tag "${fasta.baseName}"
    publishDir path: { params.save_reference ? "${params.outdir}/reference_genome" : params.outdir },
            saveAs: { params.save_reference ? it : null }, mode: params.publish_dir_mode

    input:
    path(fasta)

    output:
    path("${fasta.baseName}.index*")

    script:
    """
    bowtie-build --threads $task.cpus ${fasta} ${fasta.baseName}.index
    """
}