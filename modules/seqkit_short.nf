// Process: Seqkit stats
process SEQKIT_STATS_SHORT {
    container 'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0'
    publishDir "${params.output_dir}/short_read/01_quality_control/raw/seqkit", mode: 'copy'

    input:
    tuple val(sample_id), path(input_fastq)

    output:
    path "${sample_id}_seqkit_stats.txt", emit: stats

    script:
    """
    seqkit stats -a ${input_fastq} > ${sample_id}_seqkit_stats.txt
    """
    
    // stats: simple statistics of FASTA/Q files
    
    // -a: all statistics outputted
    
    // piped into a result file with the "sample_id" as the prefix to the output file name _seqkit_stat.txt
}

