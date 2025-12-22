// Process: DESeq2 differential expression analysis (Long Reads)
// Works for both Aim 1 (normal tissues) and Aim 2 (cancer tissues)
process DESEQ2_LONG {
    container 'https://depot.galaxyproject.org/singularity/bioconductor-deseq2:1.42.0--r43hf17093f_2'
    publishDir "${params.output_dir}/long_read/05_differential_expression/deseq2_results", mode: 'copy'

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

    # Remove annotation columns
    count_data <- count_data[, -c(1:5)]

    # Clean sample names - remove minimap2 BAM suffix
    colnames(count_data) <- gsub(".sorted.bam\$", "", colnames(count_data))

    cat("\\nSamples in count data:", paste(colnames(count_data), collapse=", "), "\\n")

    # Read sample information
    all_samples <- read.csv("${sample_info}", stringsAsFactors = FALSE)

    # Filter to only long-read samples
    all_samples <- all_samples[all_samples\$seq_type == "long", ]

    # Look at which samples are actually in the count data
    # and match them to their aim_group columns

    # Get samples that are in the count data
    samples_in_counts <- colnames(count_data)

    # Filter metadata to only samples present in count data
    all_samples <- all_samples[all_samples\$sample_id %in% samples_in_counts, ]

    if (nrow(all_samples) == 0) {
        stop("ERROR: No samples from metadata matched the count data")
    }

    cat("\\nFiltered to", nrow(all_samples), "samples present in count data\\n")

    # Now check which aim these samples belong to
    aim1_count <- sum(!is.na(all_samples\$aim1_group) & all_samples\$aim1_group != "NA")
    aim2_count <- sum(!is.na(all_samples\$aim2_group) & all_samples\$aim2_group != "NA")

    cat("\\nAuto-detecting analysis aim...\\n")
    cat("Samples with aim1_group:", aim1_count, "\\n")
    cat("Samples with aim2_group:", aim2_count, "\\n")

    # Determine which aim to use and create condition variable
    if (aim1_count > 0 && aim2_count == 0) {
        # Running Aim 1 (Normal comparison)
        cat("Detected: AIM 1 (Normal tissue comparison)\\n")
        all_samples <- all_samples[!is.na(all_samples\$aim1_group) & all_samples\$aim1_group != "NA", ]
        all_samples\$condition <- gsub("_normal", "", all_samples\$aim1_group)
        analysis_title <- "Normal Bladder vs Normal Ureter"
        comparison_type <- "normal"
    } else if (aim2_count > 0 && aim1_count == 0) {
        # Running Aim 2 (Cancer comparison)
        cat("Detected: AIM 2 (Cancer tissue comparison)\\n")
        all_samples <- all_samples[!is.na(all_samples\$aim2_group) & all_samples\$aim2_group != "NA", ]
        all_samples\$condition <- gsub("_cancer", "", all_samples\$aim2_group)
        analysis_title <- "Bladder Cancer vs Ureter Cancer"
        comparison_type <- "cancer"
    } else if (aim1_count > 0 && aim2_count > 0) {
        # Both aims present - shouldn't happen but handle it
        stop("ERROR: Both aim1_group and aim2_group have values in count data. This suggests the count file contains mixed samples from different aims. Please check your alignment/quantification step.")
    } else {
        # Neither aim has values
        stop("ERROR: No valid aim group found for samples in count data. Check sample sheet.")
    }

    # Create unique row names
    rownames(all_samples) <- all_samples\$sample_id
    all_samples\$condition <- factor(all_samples\$condition)
    cat("\\n")

    # Match samples
    common_samples <- intersect(colnames(count_data), rownames(all_samples))

    # Validation
    if(length(common_samples) < 2) {
        stop("ERROR: Fewer than 2 samples matched between count data and filtered metadata")
    }

    cat("\\nSample matching summary:\\n")
    cat("Samples in count data:", length(colnames(count_data)), "\\n")
    cat("Samples in filtered metadata:", length(rownames(all_samples)), "\\n")
    cat("Matched samples:", length(common_samples), "\\n")
    cat("Matched sample IDs:", paste(common_samples, collapse=", "), "\\n\\n")

    # Subset data
    count_data <- count_data[, common_samples, drop = FALSE]
    sample_info_filtered <- all_samples[common_samples, , drop = FALSE]

    if(length(unique(sample_info_filtered\$condition)) < 2) {
        stop("ERROR: Need at least 2 different conditions")
    }

    cat("Conditions:", paste(unique(sample_info_filtered\$condition), collapse=", "), "\\n")
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
    # This means positive LogFC = upregulated in bladder.
    dds\$condition <- relevel(dds\$condition, ref = "ureter")

    dds <- dds[rowSums(counts(dds)) >= 10, ]
    cat("Genes after filtering:", nrow(dds), "\\n\\n")

    # Run analysis
    dds <- DESeq(dds)
    res <- results(dds)

    # Save results
    write.csv(as.data.frame(res[order(res\$padj), ]),
              "differential_expression_results.csv")
    write.csv(counts(dds, normalized = TRUE),
              "normalised_counts.csv")

    # Ranked gene list
    res_df <- as.data.frame(res)
    res_df\$gene <- rownames(res_df)
    res_df <- res_df[!is.na(res_df\$log2FoldChange) & !is.na(res_df\$padj), ]
    res_df\$rank_metric <- sign(res_df\$log2FoldChange) * -log10(res_df\$padj)
    res_df <- res_df[order(res_df\$rank_metric, decreasing = TRUE), ]
    write.table(res_df[, c("gene", "rank_metric")],
                file = "ranked_genes.rnk",
                quote = FALSE, sep = "\\t",
                row.names = FALSE, col.names = FALSE)

    # Plots
    vsd <- vst(dds, blind = FALSE)
    pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
    percent_var <- round(100 * attr(pca_data, "percentVar"))

    pdf("pca_plot.pdf", width = 8, height = 6)
    p <- ggplot(pca_data, aes(x = PC1, y = PC2,
                              color = condition, label = name)) +
        geom_point(size = 4) +
        geom_text(vjust = -1, size = 3) +
        xlab(paste0("PC1: ", percent_var[1], "% variance")) +
        ylab(paste0("PC2: ", percent_var[2], "% variance")) +
        theme_bw() +
        ggtitle("Principal Component Analysis",
                subtitle = analysis_title) +
        theme(plot.title = element_text(hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5))
    print(p)
    dev.off()

    pdf("ma_plot.pdf", width = 8, height = 6)
    plotMA(res, main = paste0("MA Plot: ", analysis_title),
           ylim = c(-5, 5), alpha = 0.05)
    dev.off()

    # Summary
    cat("\\n=== Differential Expression Summary ===\\n")
    cat("Analysis type:", comparison_type, "\\n")
    cat("Comparison:", analysis_title, "\\n")
    cat("Total genes:", nrow(res), "\\n")
    cat("Significant (padj < 0.05):",
        sum(res\$padj < 0.05, na.rm = TRUE), "\\n")
    cat("Upregulated:",
        sum(res\$padj < 0.05 & res\$log2FoldChange > 0, na.rm = TRUE), "\\n")
    cat("Downregulated:",
        sum(res\$padj < 0.05 & res\$log2FoldChange < 0, na.rm = TRUE), "\\n")
    """
}
