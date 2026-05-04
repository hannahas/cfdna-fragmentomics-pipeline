/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MULTIQC
    Aggregate QC reports from FastQC and Trim Galore into a single report.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process MULTIQC {
    label 'process_single'

    container 'quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0'

    input:
    path fastqc_zips
    path trimgalore_logs

    output:
    path "multiqc_report.html",     emit: report
    path "multiqc_report_data/",    emit: data

    script:
    """
    multiqc . --filename multiqc_report.html
    """
}