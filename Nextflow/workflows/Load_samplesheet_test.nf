nextflow_workflow {
    /*
    name "Test Workflow SAMPLESHEET_TO_CHANNEL"
    script "subworkflows/local/samplesheet_to_channel/main.nf"
    workflow "SAMPLESHEET_TO_CHANNEL"
    */
    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                skip_tools = 'baserecalibrator'

            }
            workflow {
                """
                // define inputs of the workflow here. Example:
                input[0] = Channel.of([['patient':'test', 'sample':'test',
                                        'status':0, 'lane':'test_L1'],
                                        file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
                                        file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)])
                """
            }
        }

        then {
            assert workflow.success
            assert snapshot(workflow.out).match()
        }

    }

}