// Process: SeqKit stats on trimmed/processed long reads
process SEQKIT_TRIMMED_LONG {
    container 'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0'
    publishDir "${params.output_dir}/long_read/01_quality_control/trimmed/seqkit", mode: 'copy'

    input:
    tuple val(sample_id), path(input_fastq)

    output:
    path "${sample_id}_seqkit_stats.txt", emit: stats

    script:
    """
    seqkit stats -a ${input_fastq} > ${sample_id}_seqkit_stats.txt
    """
}
