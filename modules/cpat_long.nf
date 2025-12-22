process CPAT_LONG {
    container 'https://depot.galaxyproject.org/singularity/cpat:3.0.4--py37h77a2a36_0'
    // CHANGE: Path points to long_read
    publishDir "${params.output_dir}/long_read/06_biological_analysis/functional_annotation", mode: 'copy'

    input:
    path fasta_seqs       // isoform_switch_nt.fasta
    path hexamer_tab      // Human_Hexamer.tsv
    path logit_model      // Human_logitModel.RData

    output:
    path "cpat_results.txt", emit: results

    script:
    """
    # Run CPAT
    cpat.py \\
        -g ${fasta_seqs} \\
        -x ${hexamer_tab} \\
        -d ${logit_model} \\
        --top-orf=5 \\
        -o cpat_output

    # Rename to standardised output
    if [ -f cpat_output.ORF_prob.best.tsv ]; then
        mv cpat_output.ORF_prob.best.tsv cpat_results.txt
    else
        echo "Error: CPAT output file not found" >&2
        exit 1
    fi
    """
}
