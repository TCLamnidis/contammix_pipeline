//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:'\t' )
        .map { create_bam_channel(it) }
        .dump(tag: "samplesheetcheck_bams")
        .set { reads }

    emit:
    reads                                     // channel: [[ val(meta), reads] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, bam ]
def create_bam_channel(LinkedHashMap row) {
    def meta = [:]
    // TODO create spanning main metadata
    meta.id                 = [ row.sample_id, row.library_id ].join("_").trim()

    meta.sample_id          = row.sample_id

    def array = []
    if (!file(row.bam).exists()) {
        exit 1, "[nf-core/eager] error: Please check input samplesheet. BAM file does not exist!\nFile: ${row.bam}"
    } else {
        array = [ meta, file(row.bam) ]
    }
    return array
}
