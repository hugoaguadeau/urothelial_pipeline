// Process: StringTie for transcript quantification (Long Read version)
process STRINGTIE_LONG {
    tag "$sample_id"
    container 'https://depot.galaxyproject.org/singularity/stringtie:2.2.1--hecb563c_2'
    publishDir "${params.output_dir}/long_read/04_quantification/stringtie/individual", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam)
    path gtf
    
    output:
    tuple val(sample_id), path("${sample_id}_transcripts.gtf"), emit: gtf
    tuple val(sample_id), path("${sample_id}_abundances.txt"), emit: abundances
    path "${sample_id}_transcripts.gtf", emit: gtf_for_merge
    
    script:
    """
    stringtie \\
        -L \\
        -p ${task.cpus} \\
        -G ${gtf} \\
        -A ${sample_id}_abundances.txt \\
        -o ${sample_id}_transcripts.gtf \\
        ${bam}
    """
}
