// Process: StringTie Re-quantification with merged reference
process STRINGTIE_REQUANT_SHORT {
    container 'https://depot.galaxyproject.org/singularity/stringtie:2.2.1--hecb563c_2'
    publishDir "${params.output_dir}/short_read/04_quantification/stringtie/requantified", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bam_index)
    path merged_gtf

    output:
    tuple val(sample_id), path("${sample_id}_requant.gtf"), emit: gtf
    path "${sample_id}_requant_gene_abundance.tab", emit: gene_abundance
    path "${sample_id}_requant_transcript_abundance.tab", emit: transcript_abundance

    script:
    """
    stringtie \\
        ${bam} \\
        -G ${merged_gtf} \\
        -o ${sample_id}_requant.gtf \\
        -p ${task.cpus} \\
        -A ${sample_id}_requant_gene_abundance.tab \\
        -C ${sample_id}_requant_transcript_abundance.tab \\
        -e \\
        -B \\
        -v
    """

    // This is the same as the initial StringTie run, but now uses
    // the merged GTF as reference instead of the original annotation.
}
