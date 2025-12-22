// Process: CPAT (Coding Potential Assessment Tool)
process CPAT_SHORT {
    container 'https://depot.galaxyproject.org/singularity/cpat:3.0.4--py37h77a2a36_0'
    publishDir "${params.output_dir}/short_read/06_biological_analysis/functional_annotation", mode: 'copy'

    input:
    path fasta_seqs       // isoform_switch_nt.fasta
    path hexamer_tab      // Human_Hexamer.tsv
    path logit_model      // Human_logitModel.RData

    output:
    path "cpat_results.txt", emit: results

    script:
    """
    # Run CPAT as recommended in manual
    cpat.py \\
        -g ${fasta_seqs} \\
        -x ${hexamer_tab} \\
        -d ${logit_model} \\
        --top-orf=5 \\
        -o cpat_output

    # The manual expects the .ORF_prob.best.tsv file
    # Rename to standardised output
    if [ -f cpat_output.ORF_prob.best.tsv ]; then
        mv cpat_output.ORF_prob.best.tsv cpat_results.txt
    else
        echo "Error: CPAT output file not found" >&2
        exit 1
    fi
    """
}
