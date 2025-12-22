// Process: HMMER/Pfam (Protein Domain Analysis)
process PFAM_SHORT {
    container 'https://depot.galaxyproject.org/singularity/hmmer:3.3.2--h87f3376_2'
    publishDir "${params.output_dir}/short_read/06_biological_analysis/functional_annotation", mode: 'copy'

    input:
    path fasta_seqs      // isoform_switch_AA.fasta
    path pfam_db         // Pfam-A.hmm

    output:
    path "pfam_results.txt", emit: results

    script:
    """
    # Prepare the Pfam database (required by hmmscan)
    hmmpress ${pfam_db}

    # Run hmmscan as recommended in manual
    # Use --domtblout for domain table output (parseable format)
    hmmscan \\
        --domtblout pfam_results.txt \\
        --noali \\
        -E 1e-3 \\
        --domE 1e-3 \\
        --cpu ${task.cpus} \\
        ${pfam_db} \\
        ${fasta_seqs} \\
        > pfam_scan.out

    # Verify output was created
    if [ ! -f pfam_results.txt ]; then
        echo "Error: Pfam results file not created" >&2
        exit 1
    fi
    """
}
