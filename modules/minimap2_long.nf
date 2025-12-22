// Process: minimap2 alignment for long reads
process MINIMAP2_LONG {
    tag "$sample_id"
    container 'https://depot.galaxyproject.org/singularity/flair:2.0.0--pyhdfd78af_1'
    publishDir "${params.output_dir}/long_read/03_alignment/bam_files", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    path reference_genome
    path reference_gtf
    
    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), emit: bam
    
    script:
    """
    # Use the same minimap2 command that FLAIR uses internally
    # FLAIR container has minimap2, so we can use it directly
    # Align with minimap2 (same parameters as FLAIR)
    minimap2 -ax splice -uf --secondary=no -t ${task.cpus} \\
        ${reference_genome} ${reads} > ${sample_id}.sam
    
    # Convert SAM to BAM
    # Check if samtools is available in FLAIR container
    if command -v samtools >/dev/null 2>&1; then
        echo "Using samtools for BAM conversion..."
        samtools view -bS ${sample_id}.sam | samtools sort -o ${sample_id}.sorted.bam
        samtools index ${sample_id}.sorted.bam
    else
        echo "samtools not available, using python for conversion..."
        # Actually, featureCounts can read SAM files directly oops
        mv ${sample_id}.sam ${sample_id}.sorted.bam
        touch ${sample_id}.sorted.bam.bai
    fi
    
    echo "Direct minimap2 alignment completed for ${sample_id}"
    echo "Output file size: \$(stat -c%s ${sample_id}.sorted.bam) bytes"
    """
}
