// Process: STAR Alignment
process STAR_ALIGN_SHORT {
    container 'https://depot.galaxyproject.org/singularity/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0'
    publishDir "${params.output_dir}/short_read/03_alignment/bam_files", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path index
    path gtf

    output:
    tuple val(sample_id), path("${sample_id}Aligned.sortedByCoord.out.bam"), emit: bam
    tuple val(sample_id), path("${sample_id}Aligned.sortedByCoord.out.bam.bai"), emit: bam_index
    path "${sample_id}Log.final.out", emit: log_final
    path "${sample_id}Log.out", emit: log_out
    path "${sample_id}Log.progress.out", emit: log_progress
    path "${sample_id}SJ.out.tab", emit: splice_junctions

    script:
    """
    STAR \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${index} \\
        --readFilesIn ${reads[0]} ${reads[1]} \\
        --readFilesCommand zcat \\
        --outFileNamePrefix ${sample_id} \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMunmapped Within \\
        --outSAMattributes Standard \\
        --sjdbGTFfile ${gtf} \\
        --quantMode GeneCounts

    # Index the BAM file
    samtools index ${sample_id}Aligned.sortedByCoord.out.bam
    """
}
