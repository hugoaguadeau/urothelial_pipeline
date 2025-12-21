// Process: Bladder Cancer Consensus Classifier (Long Read version)
process CONSENSUS_CLASSIFIER_LONG {
    container 'https://depot.galaxyproject.org/singularity/r-base:4.3.1'
    publishDir "${params.output_dir}/long_read/07_subtyping", mode: 'copy'

    input:
    path tpm_counts

    output:
    path "*_ConsensusClassifier.tsv", emit: subtypes
    path "classifier_output.txt", emit: log

    script:
    """
    #!/usr/bin/env Rscript

    # Load library from mounted path
    library(consensusMIBC, lib.loc="/r_libs")

    cat("Starting Consensus MIBC Classification\\n")
    cat("Input file: ${tpm_counts}\\n\\n")

    # Read TPM counts (CSV format with rownames)
    TPMdf <- read.csv("${tpm_counts}", row.names = 1, check.names = FALSE)
    
    cat("Samples to classify:", ncol(TPMdf), "\\n")
    cat("Genes in input:", nrow(TPMdf), "\\n\\n")

    # Log2 transform
    TPMdflog2 <- log2(TPMdf + 1)

    # Run classifier
    cat("Running getConsensusClass()...\\n")
    newCon <- getConsensusClass(TPMdflog2, gene_id = "hgnc_symbol")

    # Save results
    output_name <- gsub(".csv", "_ConsensusClassifier.tsv", "${tpm_counts}")
    write.table(
        data.frame("ID" = rownames(newCon), newCon),
        output_name,
        sep = "\\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE
    )

    cat("\\nResults saved to:", output_name, "\\n")
    cat("Classification complete\\n")

    # Create log
    sink("classifier_output.txt")
    cat("CONSENSUS MIBC CLASSIFIER LOG\\n")
    cat("=============================\\n\\n")
    cat("Input:", "${tpm_counts}", "\\n")
    cat("Samples classified:", nrow(newCon), "\\n")
    cat("Output:", output_name, "\\n\\n")
    cat("Subtype Summary:\\n")
    print(table(newCon\$consensusClass))
    cat("\\nMean separation level:", mean(newCon\$separationLevel, na.rm = TRUE), "\\n")
    sink()
    """
}
