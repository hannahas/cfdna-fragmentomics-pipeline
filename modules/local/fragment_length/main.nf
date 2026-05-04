/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FRAGMENT_LENGTH
    Extract fragment length distributions from cfDNA BAM files.
    Filters by mapping quality and fragment length window.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process FRAGMENT_LENGTH {
    tag "$meta.id"
    label 'process_single'

    publishDir "${params.outdir}/fragment_lengths", mode: 'copy'

    container 'quay.io/biocontainers/samtools:1.17--h00cdaf9_0'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.fragment_lengths.txt"), emit: lengths
    tuple val(meta), path("${meta.id}.fragment_lengths.tsv"), emit: tsv

    script:
    """
    # Extract fragment lengths from paired-end reads
    # Filter by mapping quality and fragment length window
    samtools view -f 0x2 -F 0x904 -q ${params.min_mapq} ${bam} | \
        awk '{
            len = \$9
            if (len < 0) len = -len
            if (len >= ${params.min_frag} && len <= ${params.max_frag})
                print len
        }' > ${meta.id}.fragment_lengths.txt

    # Create TSV with sample metadata for ML
    echo -e "sample\tcondition\tfragment_length" > ${meta.id}.fragment_lengths.tsv
    while read len; do
        echo -e "${meta.id}\t${meta.condition}\t\$len"
    done < ${meta.id}.fragment_lengths.txt >> ${meta.id}.fragment_lengths.tsv
    """
}