/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CFDNA FRAGMENTOMICS WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC          } from '../modules/local/fastqc/main'
include { TRIMGALORE      } from '../modules/local/trimgalore/main'
include { BWA_MEM         } from '../modules/local/bwa_mem/main'
include { SAMTOOLS_SORT   } from '../modules/local/samtools_sort/main'
include { FRAGMENT_LENGTH } from '../modules/local/fragment_length/main'
include { MULTIQC         } from '../modules/local/multiqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CFDNA_WORKFLOW {

    // Read samplesheet
    ch_samplesheet = channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample, condition: row.condition]
            [ meta, file(row.fastq_1), file(row.fastq_2) ]
        }

    // Create channel for BWA index
    ch_bwa_index = channel.fromPath(params.bwa_index, checkIfExists: true).first()

    // QC on raw reads
    FASTQC(ch_samplesheet)

    // Trim adapters
    TRIMGALORE(ch_samplesheet)

    // Align to reference genome
    BWA_MEM(TRIMGALORE.out.reads, ch_bwa_index)

    // Sort and index BAM
    SAMTOOLS_SORT(BWA_MEM.out.bam)

    // Extract fragment lengths
    FRAGMENT_LENGTH(SAMTOOLS_SORT.out.bam)

    // Aggregate QC
    MULTIQC(
        FASTQC.out.zip.map { meta, zip -> zip }.collect(),
        TRIMGALORE.out.log.map { meta, log -> log }.collect()
    )
}