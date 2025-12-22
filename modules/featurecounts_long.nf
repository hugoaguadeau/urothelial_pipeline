// Process: featureCounts for gene expression quantification (Long Read version)
process FEATURECOUNTS_LONG {
    container 'https://depot.galaxyproject.org/singularity/subread:2.0.6--he4a0461_0'
    publishDir "${params.output_dir}/long_read/04_quantification/feature_counts", mode: 'copy'
    
    input:
    path bams
    path gtf
    
    output:
    path "counts.txt", emit: counts
    path "counts.txt.summary", emit: summary
    
    script:
    """
    featureCounts \\
        -T ${task.cpus} \\
        -t exon \\
        -g gene_name \\
        -s 0 \\
        --largestOverlap \\
        -L \\
        -a ${gtf} \\
        -o counts.txt \\
        ${bams}
    """
    // -T: thread number
    // -t exon: count reads overlapping exons
    // -g gene_name: group by gene_name attribute
    // -s 0: unstranded (change to 1 or 2 if stranded)
    // --largestOverlap: assign read to feature with largest overlap
    // -L: count long reads (handles splice junctions better)
    // -a: annotation file
    // -o: output name
    // Removed -p: not needed for single-end long reads
}
