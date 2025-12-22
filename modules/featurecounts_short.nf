// Process: FeatureCounts
process FEATURECOUNTS_SHORT {
    container 'https://depot.galaxyproject.org/singularity/subread:2.0.6--he4a0461_0'
    publishDir "${params.output_dir}/short_read/04_quantification/feature_counts", mode: 'copy'

    input:
    path bam_files
    path gtf

    output:
    path "counts.txt", emit: counts
    path "counts.txt.summary", emit: summary

    script:
    """
    featureCounts \\
        -T ${task.cpus} \\
        -p \\
        -t exon \\
        -g gene_name \\
        -a ${gtf} \\
        -o counts.txt \\
        ${bam_files}
    """
    // -T: thread number
    // -p: paired-end reads
    // -t exon: count reads overlapping exons
    // -g gene_name: group by gene_name attribute (gives gene symbols like TP53)
    // -a: annotation file
    // -o: output name
}
