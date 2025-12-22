process PFAM_LONG {
    container 'https://depot.galaxyproject.org/singularity/hmmer:3.3.2--h87f3376_2'
    // CHANGE: Output path set to long_read
    publishDir "${params.output_dir}/long_read/06_biological_analysis/functional_annotation", mode: 'copy'

    input:
    path fasta_seqs
    path pfam_db

    output:
    path "pfam_results.txt", emit: results

    script:
    """
    hmmpress ${pfam_db}

    hmmscan \\
        --domtblout pfam_results.txt \\
        --noali \\
        -E 1e-3 \\
        --domE 1e-3 \\
        --cpu ${task.cpus} \\
        ${pfam_db} \\
        ${fasta_seqs} \\
        > pfam_scan.out

    if [ ! -f pfam_results.txt ]; then
        echo "Error: Pfam results file not created" >&2
        exit 1
    fi
    """
}
