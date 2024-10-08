/*
========================================================================================
                                    CAGESEQ WORKFLOW
========================================================================================
Adapted from:
https://github.com/nf-core/cageseq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2
include { SAMPLESHEET_CHECK } from './workflows/check_samplesheet/check_samplesheet.nf' 
include { FASTQC } from './workflows/FASTQC/FASTQC.nf'
include { make_bowtie_index } from './workflows/Make_Bowtie_index/'
include { make_STAR_index } from './workflows/Make_STAR_index/'
include { star } from './workflows/Alignment/'
include { bowtie } from './workflows/Alignment/'




def helpMessage() {

    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/cageseq --input '*_R1.fastq.gz' -profile docker

    Mandatory arguments:
        --input [file]                    Path to input data (must be surrounded with quotes)
        -profile [str]                    Configuration profile to use. Can use multiple (comma separated)
                                            Available: docker, singularity, conda, test, awsbatch, <institute> and more

    Trimming:
        --save_trimmed [bool]             Set to true to Save trimmed FastQ files
        --trim_ecop [bool]                Set to false to not trim the EcoP site
        --trim_linker [bool]              Set to false to not trim the linker
        --trim_5g [bool]                  Set to false to not trim the additonal G at the 5' end
        --trim_artifacts [bool]           Set to false to not trim artifacts
        --artifacts_5end [file]           Path to 5 end artifact file, if not given the pipeline will use a default file with all possible artifacts
        --artifacts_3end [file]           Path to 3 end artifact file, if not given the pipeline will use a default file with all possible artifacts

    References                            If not specified in the configuration file or you wish to overwrite any of the references
        --genome [str]                    Name of iGenomes reference
        --gtf [file]                      Path to gtf file, used to generate the STAR index, for STAR alignment and for the clustering QC

    Ribosomal RNA removal:
        --remove_ribo_rna [bool]          Removes ribosomal RNA using SortMeRNA
        --save_non_ribo_reads [bool]      Save FastQ file intermediates after removing rRNA
        --ribo_database_manifest [string] Path to file that contains file paths for rRNA databases, optional

    Alignment:
        --aligner [str]                   Specifies the aligner to use (available are: 'star', 'bowtie1')
        --star_index [file]               Path to STAR index, set to false if igenomes should be used
        --bowtie_index [file]             Path to bowtie index, set to false if igenomes should be used

    Clustering:
        --min_cluster [int]               Minimum amount of reads to build a cluster with paraclu. Default ${params.min_cluster}
        --tpm_cluster_threshold [int]     --tpm_cluster_threshold [int] Threshold for expression count of ctss considered in paraclu clustering. Default: ${params.tpm_cluster_threshold}

    Output:
        --bigwig [bool]                    Set this option to get ctss files in bigwig-format, in addition to the default in bed-format

    Skipping options:
        --skip_initial_fastqc [bool]      Skip FastQC run on input reads
        --skip_trimming [bool]            Skip all trimming steps
        --skip_trimming_fastqc [bool]     Skip FastQC run on trimmed reads
        --skip_alignment [bool]           Skip alignment step
        --skip_samtools_stats [bool]      Skip samtools stats QC step of aligned reads
        --skip_ctss_generation [bool]     Skip steps generating CTSS files including clustering, bed/bigwig and count table output generation
        --skip_ctss_qc [bool]             Skip running RSeQC's read distribution QC step on the clustered CTSS

    Other options:
        --outdir [file]                   The output directory where the results will be saved
        --publish_dir_mode [str]          Mode for publishing results in the output directory. Available: symlink, rellink, link, copy, copyNoFollow, move (Default: copy)
        --email [email]                   Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
        --email_on_fail [email]           Same as --email, except only send mail if the workflow is not successful
        --max_multiqc_email_size [str]    Threshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB)
        -name [str]                       Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic

    AWSBatch options:
        --awsqueue [str]                  The AWSBatch JobQueue that needs to be set when running on AWSBatch
        --awsregion [str]                 The AWS Region for your AWS Batch job to run on
        --awscli [str]                    Path to the AWS CLI tool
    """.stripIndent()
}

// Create configuration variables
params.star_index = params.genome ? params.genomes[ params.genome ].star ?: false : false // a param star_index is created if params.genome exists
params.bowtie_index = params.genome ? params.genomes[ params.genome ].bowtie1 ?: false : false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false


// Validate inputs

if (params.aligner != 'star' && params.aligner != 'bowtie1'){ // params taken from config
    exit 1, "Invalid aligner option: ${params.aligner}. Valid options: 'star', 'bowtie1'"
}

if(params.star_index && params.aligner == 'star'){ // index param taken from user command: see help message
    star_index = Channel
        .fromPath(params.star_index, checkIfExists: true)
        .ifEmpty(exit 1, "STAR index not found: ${params.star_index}")
}

else if(params.bowtie_index && params.aligner == 'bowtie1'){ // index param taken from user command: see help message
    bowtie_index = Channel
        .fromPath(params.bowtie_index, checkIfExists: true)
        .ifEmpty(exit 1, "Bowtie index not found: ${params.star_index}")
}

else if ( params.fasta ){
    fasta_star_index = Channel
        .fromPath(params.fasta, checkIfExists: true)

    fasta_bowtie_index = = Channel
        .fromPath(params.fasta, checkIfExists: true)
}

if( params.gtf ){
    Channel
    .fromPath(params.gtf, checkIfExists: true)
    .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
    .set{gtf_make_STAR_index, gtf_star}
} else {
    exit 1, "No GTF annotation specified! Needed for STAR and clustering QC."
}

// Check AwsBatch settings

if(workflow.profile.contains("awsbatch")){
    // Check if parameters are specified
    if(!params.awsqueue || !params.awsregion){
        exit 1, ""
    }
    // Check outdir path to be S3 bucket if running on AWS
    if(!params.outdir.startsWith('s3:')){
        exit 1, "Outdir not on a S3"
    }
    // Prevent trace files to be stored on S3 since S3 does not support rolling files
    if(params.tracedir.startsWith('s3:')){
        exit 1. "Specify a local tracedir or run without trace! S3 does not support trace files"
    }
}

/*
Import files 
uses a samplesheet if provided. Otherwise the input is given by the path in the config file
 */
if (params.samplesheet) {
    ch_read_files_fastqc = Channel
        .fromPath(params.samplesheet)

    // check samplesheet for correctness

    SAMPLESHEET_CHECK ( ch_read_files_fastqc )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it, seq_center) }
        .set { reads }

    } else {
    ch_read_files_fastqc = Channel
        .fromPath(params.input)
        .ifEmpty { exit 1, "Cannot find any reads matching: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\n" }
}


workflow {
    
    /*
    /* STEP 1 - FastQC
    */
    
    if(!params.skip_initial_fastqc){

        fastqc_results = FASTQC(ch_read_files_fastqc)
    }
    else {
        fastqc_results = Channel.empty()
    }

    /*
    STEP 2 - Make index files
    */
    // STAR index
    if(params.aligner == 'star' && !params.star_index && params.fasta && !params.skip_alignment){

        star_index = make_STAR_index(fasta_star_index, gtf_make_STAR_index.collect())
    }
    // Bowtie index
    if(params.aligner == 'bowtie1' && !params.bowtie_index && params.fasta && !params.skip_alignment){

        bowtie_index = make_bowtie_index(fasta_bowtie_index, gtf_make_Bowtie_index.collect())
    }

    /*
    STEP 3 - Mapping
    */

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