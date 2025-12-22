// Process: StringTie for transcript assembly and quantification
process STRINGTIE_SHORT {
    container 'https://depot.galaxyproject.org/singularity/stringtie:2.2.1--hecb563c_2'
    publishDir "${params.output_dir}/short_read/04_quantification/stringtie/individual_assemblies", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bam_index)
    path reference_gtf

    output:
    tuple val(sample_id), path("${sample_id}.gtf"), emit: gtf
    path "${sample_id}_gene_abundance.tab", emit: gene_abundance
    path "${sample_id}_transcript_abundance.tab", emit: transcript_abundance

    script:
    """
    stringtie \\
        ${bam} \\
        -G ${reference_gtf} \\
        -o ${sample_id}.gtf \\
        -p ${task.cpus} \\
        -A ${sample_id}_gene_abundance.tab \\
        -C ${sample_id}_transcript_abundance.tab \\
        -e \\
        -B \\
        -v
    """

    // -G: reference annotation GTF file
    // -o: output GTF file with assembled transcripts
    // -p: number of threads
    // -A: gene abundance output (gene-level expression)
    // -C: transcript/isoform abundance output  
    // -e: only estimate abundance of given reference transcripts
    //     (limits to known isoforms, more conservative)
    // -B: enable output of ballgown table files (for downstream analysis)
    // -v: verbose output for logging
}
