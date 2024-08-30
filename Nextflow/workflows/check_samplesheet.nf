/*
taken from https://github.com/nf-core/atacseq/blob/2.1.2/modules/local/samplesheet_check.nf

*/

process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_single'


    input:
    path samplesheet

    output:
    path '*.csv'       , emit: csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/atacseq/bin/
    def args = task.ext.args ?: ''
    """
    check_samplesheet.py \\
        $samplesheet \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}