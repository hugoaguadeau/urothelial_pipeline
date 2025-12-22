// Process: Generate count matrices from StringTie Ballgown output using prepDE.py
process STRINGTIE_PREPDE_LONG {
    container 'https://depot.galaxyproject.org/singularity/python:3.9'
    publishDir "${params.output_dir}/long_read/04_quantification/stringtie", mode: 'copy'
    
    input:
    path gtf_files
    
    output:
    path "transcript_count_matrix.csv", emit: transcript_counts
    path "gene_count_matrix.csv", emit: gene_counts
    path "sample_list.txt", emit: sample_list
    
    script:
    """
    # Create sample list file (sample_id <TAB> path_to_gtf)
    for gtf in *.gtf; do
        sample=\$(basename \$gtf .gtf)
        echo "\${sample}\t\${gtf}" >> sample_list.txt
    done
    
    # Copy prepDE.py3 from bin directory
    cp ${projectDir}/bin/prepDE.py3 .
    
    # Generate count matrices
    python3 prepDE.py3 \\
        -i sample_list.txt \\
        -g gene_count_matrix.csv \\
        -t transcript_count_matrix.csv \\
        -l 150
    
    # Note: -l 150 is average read length estimate
    # For long reads, this is less critical but required by prepDE
    """
}
