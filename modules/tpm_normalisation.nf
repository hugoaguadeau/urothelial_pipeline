// Process: TPM Normalisation
// Convert raw counts to TPM (Transcripts Per Million)
// This normalises for both library size and gene length
process TPM_NORMALISATION {
    container 'https://depot.galaxyproject.org/singularity/r-base:4.3.1'
    publishDir "${params.output_dir}/${seq_type}/04_quantification/tpm", mode: 'copy'

    input:
    path counts_file
    val seq_type  // "short_read" or "long_read"

    output:
    path "tpm_normalised_counts.csv", emit: tpm_counts
    path "tpm_normalisation_summary.txt", emit: summary

    script:
    """
    #!/usr/bin/env Rscript

    # Read featureCounts output
    cat("Reading count data...\\n")
    count_data <- read.table("${counts_file}", header = TRUE, row.names = 1,
                             skip = 1, check.names = FALSE)

    # Extract gene lengths (column 5)
    gene_lengths <- count_data[, "Length"]

    # Extract count columns (columns 6 onwards)
    count_matrix <- count_data[, 6:ncol(count_data)]

    # Clean sample names based on sequencing type
    if ("${seq_type}" == "short_read") {
        # Remove STAR alignment suffix for short reads
        colnames(count_matrix) <- gsub("Aligned.sortedByCoord.out.bam\$", "", 
                                       colnames(count_matrix))
    } else if ("${seq_type}" == "long_read") {
        # Remove minimap2 BAM suffix for long reads
        colnames(count_matrix) <- gsub(".sorted.bam\$", "", 
                                       colnames(count_matrix))
    }

    cat("Sample names after cleaning:\\n")
    cat(paste(colnames(count_matrix), collapse = ", "), "\\n\\n")

    # ========================================================================
    # TPM NORMALISATION
    # ========================================================================
    # TPM = Transcripts Per Million
    # Accounts for both sequencing depth AND gene length
    
    # Step 1: Calculate RPK (Reads Per Kilobase)
    # Divide each count by gene length in kb
    rpk_matrix <- count_matrix / (gene_lengths / 1000)
    
    # Step 2: Calculate scaling factors (sum of RPK per sample)
    scaling_factors <- colSums(rpk_matrix, na.rm = TRUE)
    
    # Step 3: Calculate TPM
    # Divide each RPK by the scaling factor and multiply by 1 million
    tpm_matrix <- t(t(rpk_matrix) / scaling_factors) * 1000000
    
    # Handle any NaN or Inf values (shouldn't occur, but safety check)
    tpm_matrix[is.na(tpm_matrix)] <- 0
    tpm_matrix[is.infinite(tpm_matrix)] <- 0

    # ========================================================================
    # SAVE TPM MATRIX
    # ========================================================================
    write.csv(tpm_matrix, "tpm_normalised_counts.csv", row.names = TRUE)

    # ========================================================================
    # CREATE SUMMARY REPORT
    # ========================================================================
    sink("tpm_normalisation_summary.txt")
    
    cat("================================================================================\\n")
    cat("                    TPM NORMALISATION SUMMARY\\n")
    cat("                         (", toupper("${seq_type}"), ")\\n", sep="")
    cat("================================================================================\\n\\n")
    cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n\\n")
    
    cat("INPUT DATA:\\n")
    cat("--------------------------------------------------------------------------------\\n")
    cat(sprintf("Total genes:              %s\\n", 
               format(nrow(count_matrix), big.mark = ",")))
    cat(sprintf("Total samples:            %d\\n", ncol(count_matrix)))
    cat(sprintf("Sample names:             %s\\n", 
               paste(colnames(count_matrix), collapse = ", ")))
    
    cat("\\nRAW COUNT STATISTICS:\\n")
    cat("--------------------------------------------------------------------------------\\n")
    total_reads_per_sample <- colSums(count_matrix, na.rm = TRUE)
    cat(sprintf("Mean reads per sample:    %s\\n", 
               format(round(mean(total_reads_per_sample)), big.mark = ",")))
    cat(sprintf("Min reads per sample:     %s\\n", 
               format(round(min(total_reads_per_sample)), big.mark = ",")))
    cat(sprintf("Max reads per sample:     %s\\n", 
               format(round(max(total_reads_per_sample)), big.mark = ",")))
    
    cat("\\nTPM NORMALISATION:\\n")
    cat("--------------------------------------------------------------------------------\\n")
    cat("Method: TPM (Transcripts Per Million)\\n")
    cat("  1. Divide counts by gene length (kb) → RPK\\n")
    cat("  2. Divide RPK by sum of all RPK → normalise by library size\\n")
    cat("  3. Multiply by 1,000,000 → TPM\\n")
    
    cat("\\nTPM STATISTICS:\\n")
    cat("--------------------------------------------------------------------------------\\n")
    cat(sprintf("Genes with TPM > 0:       %s (%.1f%%)\\n",
               format(sum(rowSums(tpm_matrix) > 0), big.mark = ","),
               (sum(rowSums(tpm_matrix) > 0) / nrow(tpm_matrix)) * 100))
    cat(sprintf("Genes with TPM > 1:       %s (%.1f%%)\\n",
               format(sum(rowSums(tpm_matrix) > 1), big.mark = ","),
               (sum(rowSums(tpm_matrix) > 1) / nrow(tpm_matrix)) * 100))
    cat(sprintf("Genes with TPM > 10:      %s (%.1f%%)\\n",
               format(sum(rowSums(tpm_matrix) > 10), big.mark = ","),
               (sum(rowSums(tpm_matrix) > 10) / nrow(tpm_matrix)) * 100))
    
    cat("\\nPER-SAMPLE TPM SUMMARY:\\n")
    cat("--------------------------------------------------------------------------------\\n")
    for (i in 1:ncol(tpm_matrix)) {
        sample_name <- colnames(tpm_matrix)[i]
        sample_tpm <- tpm_matrix[, i]
        expressed_genes <- sum(sample_tpm > 0)
        
        cat(sprintf("%-20s expressed genes: %s (%.1f%%)\\n",
                   sample_name,
                   format(expressed_genes, big.mark = ","),
                   (expressed_genes / nrow(tpm_matrix)) * 100))
    }
    
    cat("\\n================================================================================\\n")
    cat("TPM normalisation completed successfully\\n")
    cat("Output: tpm_normalised_counts.csv\\n")
    cat("================================================================================\\n")
    
    sink()
    
    cat("TPM normalisation complete\\n")
    """
}
