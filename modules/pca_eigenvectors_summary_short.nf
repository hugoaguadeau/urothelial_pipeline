// Process: PCA Eigenvectors Summary (Short Reads)
// Creates a human-readable text summary of the PCA eigenvector analysis results
process PCA_EIGENVECTORS_SUMMARY_SHORT {
    container 'https://depot.galaxyproject.org/singularity/r-base:4.3.1'
    publishDir "${params.output_dir}/short_read/05_differential_expression/summaries", mode: 'copy'

    input:
    path pca_loadings
    path pca_variance
    path top_genes

    output:
    path "pca_eigenvectors_summary.txt", emit: summary_txt

    script:
    """
    #!/usr/bin/env Rscript

    # ========================================================================
    # PCA EIGENVECTORS SUMMARY - SHORT READS
    # ========================================================================

    # Read the input files
    pca_loadings <- read.csv("${pca_loadings}", row.names = 1)
    pca_variance <- read.csv("${pca_variance}", row.names = 1)
    top_genes <- read.csv("${top_genes}", stringsAsFactors = FALSE)

    # Open output file
    sink("pca_eigenvectors_summary.txt")

    # ========================================================================
    # HEADER
    # ========================================================================
    cat("================================================================================\\n")
    cat("           PCA EIGENVECTORS ANALYSIS SUMMARY (SHORT READS)\\n")
    cat("================================================================================\\n\\n")
    cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n\\n")

    # ========================================================================
    # VARIANCE EXPLAINED
    # ========================================================================
    cat("================================================================================\\n")
    cat("VARIANCE EXPLAINED BY PRINCIPAL COMPONENTS\\n")
    cat("================================================================================\\n\\n")

    # Get number of PCs
    n_pcs <- ncol(pca_loadings)
    cat("Total principal components calculated:", n_pcs, "\\n\\n")

    # Show variance for first 10 PCs (or fewer if less available)
    n_show <- min(10, n_pcs)
    
    cat("Variance explained by first", n_show, "components:\\n")
    cat("--------------------------------------------------------------------------------\\n")
    cat(sprintf("%-6s %15s %20s %20s\\n", "PC", "Std Dev", "Prop. Variance", "Cumulative Prop."))
    cat("--------------------------------------------------------------------------------\\n")

    for (i in 1:n_show) {
        std_dev <- pca_variance[1, i]
        prop_var <- pca_variance[2, i] * 100
        cum_prop <- pca_variance[3, i] * 100
        
        cat(sprintf("%-6s %15.4f %19.2f%% %19.2f%%\\n", 
                   paste0("PC", i), std_dev, prop_var, cum_prop))
    }

    cat("\\n")

    # ========================================================================
    # TOP CONTRIBUTING GENES - PC1
    # ========================================================================
    cat("================================================================================\\n")
    cat("TOP 20 GENES CONTRIBUTING TO PC1\\n")
    cat("================================================================================\\n\\n")

    pc1_genes <- top_genes[top_genes\$PC == "PC1", ]
    
    cat("PC1 explains", sprintf("%.2f%%", pca_variance[2, 1] * 100), "of total variance\\n\\n")
    
    cat(sprintf("%-6s %-20s %15s %15s\\n", "Rank", "Gene", "Loading", "Abs(Loading)"))
    cat("--------------------------------------------------------------------------------\\n")
    
    for (i in 1:nrow(pc1_genes)) {
        cat(sprintf("%-6d %-20s %15.4f %15.4f\\n",
                   pc1_genes\$Rank[i],
                   pc1_genes\$Gene[i],
                   pc1_genes\$Loading[i],
                   pc1_genes\$Abs_Loading[i]))
    }

    cat("\\n")
    cat("Note: Loading values indicate the strength and direction of each gene's\\n")
    cat("contribution to PC1. Positive and negative values indicate opposite\\n")
    cat("directions of contribution.\\n\\n")

    # ========================================================================
    # TOP CONTRIBUTING GENES - PC2
    # ========================================================================
    cat("================================================================================\\n")
    cat("TOP 20 GENES CONTRIBUTING TO PC2\\n")
    cat("================================================================================\\n\\n")

    pc2_genes <- top_genes[top_genes\$PC == "PC2", ]
    
    cat("PC2 explains", sprintf("%.2f%%", pca_variance[2, 2] * 100), "of total variance\\n\\n")
    
    cat(sprintf("%-6s %-20s %15s %15s\\n", "Rank", "Gene", "Loading", "Abs(Loading)"))
    cat("--------------------------------------------------------------------------------\\n")
    
    for (i in 1:nrow(pc2_genes)) {
        cat(sprintf("%-6d %-20s %15.4f %15.4f\\n",
                   pc2_genes\$Rank[i],
                   pc2_genes\$Gene[i],
                   pc2_genes\$Loading[i],
                   pc2_genes\$Abs_Loading[i]))
    }

    cat("\\n")
    cat("Note: Loading values indicate the strength and direction of each gene's\\n")
    cat("contribution to PC2. Positive and negative values indicate opposite\\n")
    cat("directions of contribution.\\n\\n")

    # ========================================================================
    # GENES APPEARING IN BOTH PC1 AND PC2 TOP 20
    # ========================================================================
    cat("================================================================================\\n")
    cat("GENES APPEARING IN TOP 20 FOR BOTH PC1 AND PC2\\n")
    cat("================================================================================\\n\\n")

    common_genes <- intersect(pc1_genes\$Gene, pc2_genes\$Gene)
    
    if (length(common_genes) > 0) {
        cat("Number of genes appearing in both lists:", length(common_genes), "\\n\\n")
        
        cat(sprintf("%-20s %15s %15s\\n", "Gene", "PC1 Loading", "PC2 Loading"))
        cat("--------------------------------------------------------------------------------\\n")
        
        for (gene in common_genes) {
            pc1_loading <- pc1_genes\$Loading[pc1_genes\$Gene == gene]
            pc2_loading <- pc2_genes\$Loading[pc2_genes\$Gene == gene]
            
            cat(sprintf("%-20s %15.4f %15.4f\\n", gene, pc1_loading, pc2_loading))
        }
    } else {
        cat("No genes appear in the top 20 for both PC1 and PC2.\\n")
    }

    cat("\\n")

    # ========================================================================
    # OVERALL STATISTICS
    # ========================================================================
    cat("================================================================================\\n")
    cat("OVERALL STATISTICS\\n")
    cat("================================================================================\\n\\n")

    cat("Total genes analysed:", nrow(pca_loadings), "\\n")
    cat("Total samples analysed:", ncol(pca_variance), "principal components calculated\\n")
    cat("\\n")

    # Distribution of loadings for PC1
    pc1_all_loadings <- pca_loadings[, 1]
    cat("PC1 Loading Distribution:\\n")
    cat("  Minimum:  ", sprintf("%.4f", min(pc1_all_loadings)), "\\n")
    cat("  Q1:       ", sprintf("%.4f", quantile(pc1_all_loadings, 0.25)), "\\n")
    cat("  Median:   ", sprintf("%.4f", median(pc1_all_loadings)), "\\n")
    cat("  Q3:       ", sprintf("%.4f", quantile(pc1_all_loadings, 0.75)), "\\n")
    cat("  Maximum:  ", sprintf("%.4f", max(pc1_all_loadings)), "\\n")
    cat("\\n")

    # Distribution of loadings for PC2
    pc2_all_loadings <- pca_loadings[, 2]
    cat("PC2 Loading Distribution:\\n")
    cat("  Minimum:  ", sprintf("%.4f", min(pc2_all_loadings)), "\\n")
    cat("  Q1:       ", sprintf("%.4f", quantile(pc2_all_loadings, 0.25)), "\\n")
    cat("  Median:   ", sprintf("%.4f", median(pc2_all_loadings)), "\\n")
    cat("  Q3:       ", sprintf("%.4f", quantile(pc2_all_loadings, 0.75)), "\\n")
    cat("  Maximum:  ", sprintf("%.4f", max(pc2_all_loadings)), "\\n")
    cat("\\n")

    # ========================================================================
    # FILES GENERATED
    # ========================================================================
    cat("================================================================================\\n")
    cat("GENERATED FILES\\n")
    cat("================================================================================\\n\\n")

    cat("This summary is based on the following files:\\n")
    cat("  - pca_loadings.csv (all gene loadings for all PCs)\\n")
    cat("  - pca_explained_variance.csv (variance explained by each PC)\\n")
    cat("  - top_contributing_genes.csv (top 20 genes for PC1 and PC2)\\n")
    cat("\\n")

    cat("Additional visualisations generated:\\n")
    cat("  - pca_basic.pdf/png (standard PCA plot)\\n")
    cat("  - pca_labelled.pdf/png (PCA plot with sample labels)\\n")
    cat("  - top_contributing_genes.pdf/png (bar chart of gene contributions)\\n")
    cat("  - biplot_samples_genes.pdf/png (samples and genes together)\\n")
    cat("\\n")

    # ========================================================================
    # FOOTER
    # ========================================================================
    cat("================================================================================\\n")
    cat("END OF PCA EIGENVECTORS SUMMARY\\n")
    cat("================================================================================\\n")

    sink()

    cat("\\nPCA eigenvectors summary generated successfully\\n")
    """
}
