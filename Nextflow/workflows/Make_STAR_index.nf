/*
 * STEP 2  - Build STAR/bowtie1 index
 */


process make_STAR_index {
    label 'process_high'
    tag "${fasta.baseName}"
    publishDir path: { params.save_reference ? "${params.outdir}/reference_genome" : params.outdir },
            saveAs: { params.save_reference ? it : null }, mode: params.publish_dir_mode

    input:
    path(fasta)
    path(gtf)

    output:
    file "star"

    script:
    def avail_mem = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
    """
    mkdir star
    STAR \\
        --runMode genomeGenerate \\
        --runThreadN $task.cpus \\
        --sjdbGTFfile $gtf \\
        --genomeDir star/ \\
        --genomeFastaFiles $fasta \\
        $avail_mem
    """
}