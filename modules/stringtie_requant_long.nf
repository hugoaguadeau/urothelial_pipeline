// Process: Re-quantify with merged GTF and generate Ballgown tables
process STRINGTIE_REQUANT_LONG {
    tag "$sample_id"
    container 'https://depot.galaxyproject.org/singularity/stringtie:2.2.1--hecb563c_2'
    publishDir "${params.output_dir}/long_read/04_quantification/stringtie/ballgown/${sample_id}", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam)
    path merged_gtf
    
    output:
    tuple val(sample_id), path("${sample_id}.gtf"), emit: gtf
    tuple val(sample_id), path("*.ctab"), emit: ballgown_tables
    path "${sample_id}.gtf", emit: gtf_for_prepde
    
    script:
    """
    stringtie \\
        -e \\
        -B \\
        -L \\
        -p ${task.cpus} \\
        -G ${merged_gtf} \\
        -o ${sample_id}.gtf \\
        ${bam}
    """
    // -e: only estimate abundance of given reference transcripts (from merged GTF)
    // -B: enable output of Ballgown table files (*.ctab) - REQUIRED for prepDE.py
    // -L: long reads mode
    // -G: reference annotation (merged GTF)
}
