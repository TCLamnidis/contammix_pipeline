/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowContammix.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
// ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK      } from '../subworkflows/local/input_check'
include { RUN_CONTAMMIX    } from '../modules/local/run_contammix'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CAT_CAT                     } from '../modules/nf-core/modules/cat/cat/main'
include { SAMTOOLS_FASTQ              } from '../modules/nf-core/modules/samtools/fastq/main'
include { BWA_INDEX                   } from '../modules/nf-core/modules/bwa/index/main'
include { BWA_ALN                     } from '../modules/nf-core/modules/bwa/aln/main'
include { BWA_SAMSE                   } from '../modules/nf-core/modules/bwa/samse/main'
include { IVAR_CONSENSUS              } from '../modules/nf-core/modules/ivar/consensus/main'
include { MAFFT                       } from '../modules/nf-core/modules/mafft/main'
// include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow CONTAMMIX {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    SAMTOOLS_FASTQ (
        INPUT_CHECK.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FASTQ.out.versions)

    IVAR_CONSENSUS (
        INPUT_CHECK.out.bam, false
    )
    ch_versions = ch_versions.mix(IVAR_CONSENSUS.out.versions)

    CAT_CAT (
        IVAR_CONSENSUS.out.fasta
            .combine(Channel.fromList(["${baseDir}/assets/311mt_genomes_MH.fas"]))
            .map{
                    meta, bam, fasta ->

                    [meta, [bam,fasta]]
                }
                // .dump(tag:"cat_inputs")
    )
    ch_versions = ch_versions.mix(CAT_CAT.out.versions)

    MAFFT (
        CAT_CAT.out.file_out
        .dump(tag:"mafft_input")
    )
    ch_versions = ch_versions.mix(MAFFT.out.versions)

    BWA_INDEX (
        IVAR_CONSENSUS.out.fasta
    )
    ch_versions = ch_versions.mix(BWA_INDEX.out.versions)

    //Meta requires 'single_end' attribute for bwa aln
    ch_fastq_updated_meta_with_index = SAMTOOLS_FASTQ.out.fastq
            .join(
                BWA_INDEX.out.index
            )
            // editing the meta map directly is a bad idea. clone it then edit the clone. See: https://midnighter.github.io/nextflow-gotchas/gotchas/shallow-copy/
            .map{
                meta, fastq, index ->
                clone = meta.clone()
                clone['single_end'] = true

                [clone, fastq, index]
            }

    // For BWA_ALN, sinmply split channel in two.
    ch_input_bwa_aln = ch_fastq_updated_meta_with_index
            .multiMap{
                it ->
                fastq: [it[0], it [1]]
                index: it[2]
            }

    BWA_ALN(
        ch_input_bwa_aln.fastq,
        ch_input_bwa_aln.index
    )
    ch_versions = ch_versions.mix(BWA_ALN.out.versions)

    //For BWA_SAMSE add BWA_ALN output then split in two
    ch_input_bwa_samse = ch_fastq_updated_meta_with_index
            .join(BWA_ALN.out.sai)
            .dump(tag:"bwa_samse_input_before multiMap")
            .multiMap{
                meta, fastq, index, sai ->
                fastq: [meta, fastq, sai]
                index: [index]
            }

    //Meta requires 'single_end' attribute for bwa samse
    BWA_SAMSE (
        ch_input_bwa_samse.fastq,
        ch_input_bwa_samse.index
    )
    ch_versions = ch_versions.mix(BWA_SAMSE.out.versions)

    ch_contammix_input = BWA_SAMSE.out.bam
        .join(
            MAFFT.out.fas
                .map{
                meta, fastq ->
                def clone = meta.clone()
                clone['single_end'] = true

                [clone, fastq]
            }
        )
        .dump(tag:"contammix_input")

    RUN_CONTAMMIX (
        ch_contammix_input
    )
    ch_versions = ch_versions.mix(RUN_CONTAMMIX.out.versions)

    // CUSTOM_DUMPSOFTWAREVERSIONS (
    //     ch_versions.unique().collectFile(name: 'collated_versions.yml')
    // )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
