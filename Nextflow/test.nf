nextflow.enable.dsl = 2

include { SAMPLESHEET_CHECK } from './workflows/check_samplesheet.nf'

workflow INPUT_CHECK {
    take:
    samplesheet  // file: /path/to/samplesheet.csv
    seq_center   // string: sequencing center for read group
    with_control // boolean: samplesheet contains controls


    main:
    SAMPLESHEET_CHECK ( samplesheet, seq_center )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it, seq_center) }
        .set { reads }

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

ch_input = Channel.
    fromPath('./data/samplesheet.csv')

workflow {
    INPUT_CHECK (
        ch_input,
        'file_out',
        false
    )
}
