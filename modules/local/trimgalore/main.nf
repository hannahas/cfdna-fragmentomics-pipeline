/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TRIMGALORE
    Trim adapters and low quality bases from raw reads.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process TRIMGALORE {
    tag "$meta.id"
    label 'process_single'

    container 'quay.io/biocontainers/trim-galore:0.6.7--hdfd78af_0'

    input:
    tuple val(meta), path(reads_1), path(reads_2)

    output:
    tuple val(meta), path("*_val_1.fq.gz"), path("*_val_2.fq.gz"), emit: reads
    tuple val(meta), path("*trimming_report.txt"),                  emit: log

    script:
    """
    trim_galore --paired --gzip ${reads_1} ${reads_2}
    """
}