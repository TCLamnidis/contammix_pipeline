process RUN_CONTAMMIX {
    tag "$meta.id"
    label 'process_medium'

    // TODO Change this to local copy of container? (once that is ready)
    container "${params.contammix_container_path}"

    input:
    tuple val(meta), path(bam), path(alignment)

    output:
    tuple val(meta), path("*.pdf"), emit: pdf
    tuple val(meta), path("*.txt"), emit: txt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    contammix \\
        --nrThreads ${task.cpus} \\
        --samFn ${bam} \\
        --malnFn ${alignment} \\
        --figure ${prefix}.pdf \\
        ${args} \\
        > ${prefix}.txt

    ## TODO change path to NEWS file if/when a full public container is made.
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        contammix: \$(head -n1 /contamMix/NEWS | cut -f1 -d " " | sed 's/^v//' | sed 's/ (.*)//')
    END_VERSIONS
    """
}
