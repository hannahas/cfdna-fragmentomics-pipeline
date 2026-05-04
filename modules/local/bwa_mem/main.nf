/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BWA_MEM
    Align trimmed reads to reference genome using BWA MEM.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process BWA_MEM {
    tag "$meta.id"
    label 'process_medium'

    container 'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:219b6c272b25e7e642ae3ff0bf0c5c81a5135ab4-0'

    input:
    tuple val(meta), path(reads_1), path(reads_2)
    path bwa_index

    output:
    tuple val(meta), path("${meta.id}.bam"), emit: bam

    script:
    """
    bwa mem \\
        -t ${task.cpus} \\
        -R "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA" \\
        ${bwa_index}/${params.genome}.fa \\
        ${reads_1} \\
        ${reads_2} \\
        | samtools view -bS - > ${meta.id}.bam
    """
}