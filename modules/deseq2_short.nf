// Process: DESeq2 differential expression analysis (Short Reads)
process DESEQ2_SHORT {
    container 'https://depot.galaxyproject.org/singularity/bioconductor-deseq2:1.42.0--r43hf17093f_2'
    publishDir "${params.output_dir}/short_read/05_differential_expression/deseq2_results", mode: 'copy'

    input:
    path counts_file
    path sample_info

    output:
    path "differential_expression_results.csv", emit: results
    path "normalised_counts.csv", emit: norm_counts
    path "pca_plot.pdf", emit: pca_plot
    path "ma_plot.pdf", emit: ma_plot
    path "ranked_genes.rnk", emit: ranked_genes

    script:
    """
    #!/usr/bin/env Rscript

    library(DESeq2)
    library(ggplot2)

    # Read count data
    count_data <- read.table("${counts_file}", header = TRUE, row.names = 1,
                            skip = 1, check.names = FALSE)

    # Remove annotation columns (first 5 columns: Chr, Start, End, Strand, Length)
    count_data <- count_data[, -c(1:5)]

    # Clean sample names - remove STAR alignment suffix
    colnames(count_data) <- gsub("Aligned.sortedByCoord.out.bam\$", "",
                                colnames(count_data))

    # Read sample information
    all_samples <- read.csv("${sample_info}", stringsAsFactors = FALSE)

    # Filter to only short-read samples
    all_samples <- all_samples[all_samples\$seq_type == "short", ]

    # Create unique row names using sample_id (now unique after filtering)
    rownames(all_samples) <- all_samples\$sample_id

    # Ensure condition column exists and convert aim1_group to condition
    # For Aim 1: bladder_normal vs ureter_normal becomes bladder vs ureter
    all_samples\$condition <- gsub("_normal", "", all_samples\$aim1_group)
    all_samples\$condition <- factor(all_samples\$condition)

    # Match samples between count data and metadata
    common_samples <- intersect(colnames(count_data), rownames(all_samples))

    # Validation check
    if(length(common_samples) < 2) {
        stop("ERROR: Fewer than 2 samples matched between count data and sample info")
    }

    cat("\\nSample matching summary:\\n")
    cat("Samples in count data:", length(colnames(count_data)), "\\n")
    cat("Samples in metadata:", length(rownames(all_samples)), "\\n")
    cat("Matched samples:", length(common_samples), "\\n")
    cat("Matched sample IDs:", paste(common_samples, collapse=", "), "\\n\\n")

    # Subset to matched samples only
    count_data <- count_data[, common_samples, drop = FALSE]
    sample_info_filtered <- all_samples[common_samples, , drop = FALSE]

    # Verify we have multiple conditions
    if(length(unique(sample_info_filtered\$condition)) < 2) {
        stop("ERROR: Need at least 2 different conditions for differential expression")
    }

    cat("Conditions found:", paste(unique(sample_info_filtered\$condition), collapse=", "), "\\n")
    cat("Samples per condition:\\n")
    print(table(sample_info_filtered\$condition))
    cat("\\n")

    # Create DESeq2 dataset
    dds <- DESeqDataSetFromMatrix(
        countData = count_data,
        colData = sample_info_filtered,
        design = ~ condition
    )

    # Set Reference Level
    # Force "ureter" to be the reference (denominator).
    # This means positive LogFC = upregulated in bladder
    dds\$condition <- relevel(dds\$condition, ref = "ureter")

    # Filter low-count genes (keep genes with at least 10 reads total)
    dds <- dds[rowSums(counts(dds)) >= 10, ]

    cat("Genes retained after filtering:", nrow(dds), "\\n\\n")

    # Run differential expression analysis
    dds <- DESeq(dds)
    res <- results(dds)

    # Sort by adjusted p-value and save
    res_ordered <- res[order(res\$padj), ]
    write.csv(as.data.frame(res_ordered),
              "differential_expression_results.csv")

    # Save normalised counts
    write.csv(counts(dds, normalized = TRUE),
              "normalised_counts.csv")

    # Create ranked gene list for GSEA
    res_df <- as.data.frame(res)
    res_df\$gene <- rownames(res_df)

    # Remove genes with NA values
    res_df <- res_df[!is.na(res_df\$log2FoldChange) & !is.na(res_df\$padj), ]

    # Calculate ranking metric: sign(log2FC) * -log10(padj)
    res_df\$rank_metric <- sign(res_df\$log2FoldChange) * -log10(res_df\$padj)
    res_df <- res_df[order(res_df\$rank_metric, decreasing = TRUE), ]

    write.table(res_df[, c("gene", "rank_metric")],
                file = "ranked_genes.rnk",
                quote = FALSE,
                sep = "\\t",
                row.names = FALSE,
                col.names = FALSE)

    # PCA plot
    # Use variance stabilising transformation for visualisation
    vsd <- vst(dds, blind = FALSE)
    pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
    percent_var <- round(100 * attr(pca_data, "percentVar"))

    pdf("pca_plot.pdf", width = 8, height = 6)
    p <- ggplot(pca_data, aes(x = PC1, y = PC2,
                              color = condition,
                              label = name)) +
        geom_point(size = 4) +
        geom_text(vjust = -1, size = 3) +
        xlab(paste0("PC1: ", percent_var[1], "% variance")) +
        ylab(paste0("PC2: ", percent_var[2], "% variance")) +
        theme_bw() +
        ggtitle("Principal Component Analysis",
                subtitle = "Normal Bladder vs Normal Ureter") +
        theme(plot.title = element_text(hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5))
    print(p)
    dev.off()

    # MA plot
    pdf("ma_plot.pdf", width = 8, height = 6)
    plotMA(res, main = "MA Plot: Bladder vs Ureter",
           ylim = c(-5, 5),
           alpha = 0.05)
    dev.off()

    # Summary statistics
    cat("\\n=== Differential Expression Summary ===\\n")
    cat("Total genes tested:", nrow(res), "\\n")
    cat("Significant genes (padj < 0.05):",
        sum(res\$padj < 0.05, na.rm = TRUE), "\\n")
    cat("Significant upregulated (padj < 0.05, log2FC > 0):",
        sum(res\$padj < 0.05 & res\$log2FoldChange > 0, na.rm = TRUE), "\\n")
    cat("Significant downregulated (padj < 0.05, log2FC < 0):",
        sum(res\$padj < 0.05 & res\$log2FoldChange < 0, na.rm = TRUE), "\\n")
    """
}
