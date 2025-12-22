// Process: DESeq2 Summary - Extract key statistics (Long Reads)
process DESEQ2_SUMMARY_LONG {
    container 'https://depot.galaxyproject.org/singularity/r-base:4.3.1'
    publishDir "${params.output_dir}/long_read/05_differential_expression/summaries", mode: 'copy'

    input:
    path deseq2_results

    output:
    path "deseq2_summary.csv", emit: summary
    path "deseq2_summary.txt", emit: summary_txt

    script:
    """
    #!/usr/bin/env Rscript

    # Read DESeq2 results
    # The first column contains gene names but gets treated as rownames
    results <- read.csv("${deseq2_results}", stringsAsFactors = FALSE, row.names = 1)

    # Remove any genes with NA padj values for accurate counting
    results_clean <- results[!is.na(results\$padj), ]

    # Calculate summary statistics
    total_genes <- nrow(results_clean)
    
    # Significant genes (padj < 0.05)
    sig_genes <- results_clean[results_clean\$padj < 0.05, ]
    num_sig <- nrow(sig_genes)
    percent_sig <- (num_sig / total_genes) * 100
    
    # Upregulated (positive log2FoldChange, padj < 0.05)
    upregulated <- sig_genes[sig_genes\$log2FoldChange > 0, ]
    num_up <- nrow(upregulated)
    percent_up <- (num_up / total_genes) * 100
    
    # Downregulated (negative log2FoldChange, padj < 0.05)
    downregulated <- sig_genes[sig_genes\$log2FoldChange < 0, ]
    num_down <- nrow(downregulated)
    percent_down <- (num_down / total_genes) * 100
    
    # Quality metrics
    mean_base_mean <- mean(results_clean\$baseMean, na.rm = TRUE)
    median_base_mean <- median(results_clean\$baseMean, na.rm = TRUE)
    
    # Log2FoldChange statistics (for significant genes only)
    if (num_sig > 0) {
        max_lfc_up <- max(sig_genes\$log2FoldChange[sig_genes\$log2FoldChange > 0], na.rm = TRUE)
        max_lfc_down <- min(sig_genes\$log2FoldChange[sig_genes\$log2FoldChange < 0], na.rm = TRUE)
        mean_lfc_up <- mean(upregulated\$log2FoldChange, na.rm = TRUE)
        mean_lfc_down <- mean(downregulated\$log2FoldChange, na.rm = TRUE)
    } else {
        max_lfc_up <- NA
        max_lfc_down <- NA
        mean_lfc_up <- NA
        mean_lfc_down <- NA
    }
    
    # Genes with NA values (removed from analysis)
    genes_with_na <- sum(is.na(read.csv("${deseq2_results}")\$padj))
    
    # Create summary data frame
    summary_df <- data.frame(
        metric = c(
            "total_genes_tested",
            "genes_with_na_removed",
            "significant_genes_padj_0.05",
            "percent_significant",
            "upregulated_genes",
            "percent_upregulated",
            "downregulated_genes",
            "percent_downregulated",
            "mean_base_mean",
            "median_base_mean",
            "max_log2fc_upregulated",
            "max_log2fc_downregulated",
            "mean_log2fc_upregulated",
            "mean_log2fc_downregulated"
        ),
        value = c(
            total_genes,
            genes_with_na,
            num_sig,
            percent_sig,
            num_up,
            percent_up,
            num_down,
            percent_down,
            mean_base_mean,
            median_base_mean,
            max_lfc_up,
            max_lfc_down,
            mean_lfc_up,
            mean_lfc_down
        ),
        stringsAsFactors = FALSE
    )
    
    # Write CSV
    write.csv(summary_df, "deseq2_summary.csv", row.names = FALSE)
    
    # ========================================================================
    # Create detailed text report
    # ========================================================================
    
    sink("deseq2_summary.txt")
    cat("================================================================================\\n")
    cat("                    DESEQ2 DIFFERENTIAL EXPRESSION SUMMARY\\n")
    cat("                              (LONG READS)\\n")
    cat("================================================================================\\n\\n")
    cat("Analysis: Normal Bladder vs Normal Ureter\\n")
    cat("Significance threshold: adjusted p-value < 0.05\\n")
    cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n\\n")
    
    cat("================================================================================\\n")
    cat("OVERALL STATISTICS\\n")
    cat("================================================================================\\n")
    cat(sprintf("Total genes tested:              %s\\n", 
               format(total_genes, big.mark = ",")))
    cat(sprintf("Genes removed (NA values):       %s\\n", 
               format(genes_with_na, big.mark = ",")))
    cat(sprintf("\\nSignificant genes (padj < 0.05): %s (%.2f%%)\\n", 
               format(num_sig, big.mark = ","), percent_sig))
    cat(sprintf("  Upregulated:                   %s (%.2f%%)\\n", 
               format(num_up, big.mark = ","), percent_up))
    cat(sprintf("  Downregulated:                 %s (%.2f%%)\\n", 
               format(num_down, big.mark = ","), percent_down))
    
    cat("\\n================================================================================\\n")
    cat("EXPRESSION LEVEL METRICS\\n")
    cat("================================================================================\\n")
    cat(sprintf("Mean base mean expression:       %.2f\\n", mean_base_mean))
    cat(sprintf("Median base mean expression:     %.2f\\n", median_base_mean))
    
    cat("\\n================================================================================\\n")
    cat("FOLD CHANGE METRICS (SIGNIFICANT GENES ONLY)\\n")
    cat("================================================================================\\n")
    
    if (num_sig > 0) {
        cat("Upregulated genes:\\n")
        cat(sprintf("  Maximum log2 fold change:      %.2f\\n", max_lfc_up))
        cat(sprintf("  Mean log2 fold change:         %.2f\\n", mean_lfc_up))
        
        cat("\\nDownregulated genes:\\n")
        cat(sprintf("  Maximum log2 fold change:      %.2f\\n", max_lfc_down))
        cat(sprintf("  Mean log2 fold change:         %.2f\\n", mean_lfc_down))
    } else {
        cat("No significant genes found.\\n")
    }
    
    cat("\\n================================================================================\\n")
    cat("TOP 10 UPREGULATED GENES (by adjusted p-value)\\n")
    cat("================================================================================\\n")
    
    if (num_up > 0) {
        top_up <- upregulated[order(upregulated\$padj), ]
        top_up <- head(top_up, 10)
        
        cat(sprintf("%-20s %12s %12s %12s\\n", 
                   "Gene", "log2FC", "padj", "baseMean"))
        cat(strrep("-", 60), "\\n")
        
        for (i in 1:nrow(top_up)) {
            cat(sprintf("%-20s %12.2f %12.2e %12.2f\\n",
                       rownames(top_up)[i],
                       top_up\$log2FoldChange[i],
                       top_up\$padj[i],
                       top_up\$baseMean[i]))
        }
    } else {
        cat("No upregulated genes found.\\n")
    }
    
    cat("\\n================================================================================\\n")
    cat("TOP 10 DOWNREGULATED GENES (by adjusted p-value)\\n")
    cat("================================================================================\\n")
    
    if (num_down > 0) {
        top_down <- downregulated[order(downregulated\$padj), ]
        top_down <- head(top_down, 10)
        
        cat(sprintf("%-20s %12s %12s %12s\\n", 
                   "Gene", "log2FC", "padj", "baseMean"))
        cat(strrep("-", 60), "\\n")
        
        for (i in 1:nrow(top_down)) {
            cat(sprintf("%-20s %12.2f %12.2e %12.2f\\n",
                       rownames(top_down)[i],
                       top_down\$log2FoldChange[i],
                       top_down\$padj[i],
                       top_down\$baseMean[i]))
        }
    } else {
        cat("No downregulated genes found.\\n")
    }
    
    cat("================================================================================\\n")
    
    sink()
    
    cat("\\nDESeq2 summary generated successfully\\n")
    cat("Total genes tested:", total_genes, "\\n")
    cat("Significant genes:", num_sig, sprintf("(%.2f%%)\\n", percent_sig))
    """
}
