/**
 * STEP 7 - STAR/bowtie alignment
 */

    process star {
        label 'process_high'
        tag "$name"
        publishDir "${params.outdir}/STAR", mode: params.publish_dir_mode,
                saveAs: {filename ->
                    if (filename.indexOf(".bam") == -1) "logs/$filename"
                    else  filename }

        input:
        path(reads)
        path(index)
        path(gtf)

        output:
        set val(name), file("*.bam") into star_aligned
        file "*.out" into star_alignment_logs
        file "*SJ.out.tab"

        script:
        """
        STAR --genomeDir $index \\
            --sjdbGTFfile $gtf \\
            --readFilesIn $reads \\
            --runThreadN $task.cpus \\
            --outSAMtype BAM SortedByCoordinate \\
            --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 \\
            --seedSearchStartLmax 20 \\
            --outFilterMismatchNmax 1 \\
            --readFilesCommand zcat \\
            --runDirPerm All_RWX \\
            --outFileNamePrefix $name \\
            --outFilterMultimapNmax 1
        """
    }
        
    process bowtie {
        label 'process_high'
        tag "$name"
        publishDir "${params.outdir}/bowtie", mode: params.publish_dir_mode,
                saveAs: {filename ->
                    if (filename.indexOf(".bam") == -1) "logs/$filename"
                    else  filename }

        input:
        path(reads) 
        path(index_array)

        output:
        tuple val(name), path("*.bam") into bam_stats, bam_aligned
        path("*.out") into bowtie_alignment_logs

        script:
        index = index_array[0].baseName - ~/.\d$/
        """
        bowtie --sam \\
            -m 1 \\
            --best \\
            --strata \\
            -k 1 \\
            --tryhard \\
            --threads $task.cpus \\
            --phred33-quals \\
            --chunkmbs 64 \\
            --seedmms 2 \\
            --seedlen 20 \\
            --maqerr 70  \\
            ${index}  \\
            -q ${reads} \\
            --un ${reads.baseName}.unAl > ${name}.sam 2> ${name}.out
            samtools sort -@ $task.cpus -o ${name}.bam ${name}.sam
        """
    }
    
    }else{
        bowtie_alignment_logs= Channel.empty()
    }
} else {
    further_processed_reads_sortmerna.into{bam_stats; bam_aligned}
    star_alignment_logs = Channel.empty()
    bowtie_alignment_logs = Channel.empty()

}