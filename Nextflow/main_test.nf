nextflow.enable.dsl=2
include { FASTQC } from './workflows/FASTQC'




/*
READ INPUT FILES FROM SAMPLESHEET
*/
workflow Input_Files_WF {
    // uses a samplesheet if provided. Otherwise the input is given by the path in the config file
    if (params.samplesheet) {
        ch_read_files_fastqc = Channel
            .fromPath(params.samplesheet)
            .splitCsv(header: true)
            .map{ row -> [row.id, [file(row.fastq1), file(row.fastq2)]]}
            .ifEmpty { exit 1, "params.samplesheet was empty - no input files supplied" }.view()
        } else {
        ch_read_files_fastqc = Channel
            .fromPath(params.input)
            .ifEmpty { exit 1, "Cannot find any reads matching: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\n" }
}
}

/*
/* STEP 1 - FastQC
 */
workflow FastQc_WF {
    if(!params.skip_initial_fastqc){

        fastqc_results = FASTQC(ch_read_files_fastqc)
    }
    else {
        fastqc_results = Channel.empty()
    }
}

workflow {
    Input_Files_WF()
    FastQc_WF
}
/*
STEP 2 - Make index files
*/


workflow Make_STAR_index_WF {
    if(params.aligner == 'star' && !params.star_index && params.fasta && !params.skip_alignment){

        star_index = make_STAR_index(fasta_star_index, gtf_make_STAR_index.collect())
    }
}
workflow Make_Bowtie_index_WF {
    if(params.aligner == 'bowtie1' && !params.bowtie_index && params.fasta && !params.skip_alignment){

        bowtie_index = make_bowtie_index(fasta_bowtie_index, gtf_make_Bowtie_index.collect())
    }
}

/*
STEP 3 - Mapping
*/

workflow ALIGNMENT {
    if(!params.skip_alignment){
        // STAR
        if (params.aligner == 'star') {
            star_aligned = star(ch_read_files_fastqc, star_index.collect(), gtf_star.collect())

        }

        // Bowtie
        if (params.aligner == 'bowtie1'){
            bowtie_aligned = bowtie(ch_read_files_fastqc, bowtie_index.collect())
        }
    
    }
    
    
}

