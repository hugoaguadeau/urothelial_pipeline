// Process: StringTie merge - combine all sample GTFs
process STRINGTIE_MERGE_LONG {
    container 'https://depot.galaxyproject.org/singularity/stringtie:2.2.1--hecb563c_2'
    publishDir "${params.output_dir}/long_read/04_quantification/stringtie", mode: 'copy'
    
    input:
    path gtf_files
    path reference_gtf
    
    output:
    path "merged_transcripts.gtf", emit: merged_gtf
    
    script:
    """
    # Create list of GTF files
    ls *.gtf > gtf_list.txt
    
    # Merge all GTFs
    stringtie --merge \\
        -p ${task.cpus} \\
        -G ${reference_gtf} \\
        -o merged_transcripts.gtf \\
        gtf_list.txt
    """
    // --merge: merge mode
    // -G: reference annotation
    // -o: output merged GTF
}
