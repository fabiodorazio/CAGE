/* STEP 1 - FastQC
 */
process FASTQC {
    tag "$name"
    label 'process_low'
    // conda "${projectDir}/environment.yml"
    publishDir "${params.outdir}/fastqc", mode: params.publish_dir_mode,
        saveAs: { filename ->
                    filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"
                }

    input:
    //tuple val(name), path(reads)
    path(reads)
    
    output:
    path "*_fastqc.{zip, html}"

    script:
    """
    fastqc --quiet --threads $task.cpus $reads
    """
}
