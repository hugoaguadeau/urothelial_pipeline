// Process: GSEA - Gene Set Enrichment Analysis (Short Reads)
// Performs GSEA using Hallmark and Reactome gene sets
process GSEA_SHORT {
    container 'https://depot.galaxyproject.org/singularity/bioconductor-fgsea:1.28.0--r43hf17093f_0'
    publishDir "${params.output_dir}/short_read/05_differential_expression/gsea", mode: 'copy'

    input:
    path ranked_genes    // .rnk file from DESeq2
    path hallmark_gmt    // Hallmark gene sets
    path reactome_gmt    // Reactome pathways

    output:
    path "gsea_hallmark_results.csv", emit: hallmark_results
    path "gsea_reactome_results.csv", emit: reactome_results
    path "gsea_hallmark_top20_plot.pdf", emit: hallmark_plot
    path "gsea_reactome_top20_plot.pdf", emit: reactome_plot
    path "gsea_combined_summary.txt", emit: summary

    script:
    """
    #!/usr/bin/env Rscript

    # Load required libraries
    suppressPackageStartupMessages({
        library(fgsea)
        library(data.table)
        library(ggplot2)
    })

    # Read ranked gene list
    ranked_list <- read.table("${ranked_genes}",
                             header = FALSE,
                             col.names = c("gene", "score"),
                             stringsAsFactors = FALSE)

    # Create named vector (required format for fgsea)
    gene_ranks <- setNames(ranked_list\$score, ranked_list\$gene)

    # Remove any duplicates (keep the one with highest absolute score)
    gene_ranks <- gene_ranks[!duplicated(names(gene_ranks))]

    # Sort by score (most positive to most negative)
    gene_ranks <- sort(gene_ranks, decreasing = TRUE)

    cat("Total genes in ranked list:", length(gene_ranks), "\\n")
    cat("Score range:", range(gene_ranks), "\\n\\n")

    # Function to read GMT files
    read_gmt <- function(gmt_file) {
        pathways <- gmtPathways(gmt_file)
        return(pathways)
    }

    # Function to run GSEA
    run_gsea <- function(gene_ranks, pathways, pathway_name) {

        cat("Running GSEA for", pathway_name, "...\\n")
        cat("Number of gene sets:", length(pathways), "\\n")

        # Run fgsea
        set.seed(42)
        gsea_results <- fgsea(pathways = pathways,
                             stats = gene_ranks,
                             minSize = 15,
                             maxSize = 500)

        # Convert to data frame and sort by adjusted p-value
        gsea_df <- as.data.frame(gsea_results)
        gsea_df <- gsea_df[order(gsea_df\$padj), ]

        # FIX: Convert 'leadingEdge' list column to string
        if ("leadingEdge" %in% colnames(gsea_df)) {
            gsea_df\$leadingEdge <- sapply(gsea_df\$leadingEdge, paste, collapse = ";")
        }

        # Add useful columns
        gsea_df\$direction <- ifelse(gsea_df\$NES > 0, "Upregulated", "Downregulated")
        gsea_df\$significant <- gsea_df\$padj < 0.05

        # Log basic stats to stdout
        cat("Significant pathways (padj < 0.05):", sum(gsea_df\$significant, na.rm = TRUE), "\\n")

        return(gsea_df)
    }

    # Function to create enrichment plot
    create_enrichment_plot <- function(gsea_df, pathway_name, output_file) {
        # Select top 20 pathways (10 up, 10 down if available)
        sig_results <- gsea_df[gsea_df\$significant, ]

        if (nrow(sig_results) == 0) {
            cat("No significant pathways found for", pathway_name, "\\n")
            pdf(output_file, width = 10, height = 6)
            plot.new()
            text(0.5, 0.5, paste("No significant pathways found\\n(padj < 0.05)"), cex = 1.5)
            dev.off()
            return()
        }

        up_pathways <- sig_results[sig_results\$NES > 0, ]
        down_pathways <- sig_results[sig_results\$NES < 0, ]

        top_up <- head(up_pathways, 10)
        top_down <- head(down_pathways, 10)

        top_pathways <- rbind(top_up, top_down)
        top_pathways <- top_pathways[order(top_pathways\$NES), ]

        # Shorten pathway names
        top_pathways\$pathway_short <- sapply(top_pathways\$pathway, function(x) {
            if (nchar(x) > 60) { paste0(substr(x, 1, 57), "...") } else { x }
        })

        p <- ggplot(top_pathways, aes(x = NES, y = reorder(pathway_short, NES), fill = direction)) +
            geom_bar(stat = "identity") +
            scale_fill_manual(values = c("Upregulated" = "#E74C3C",
                                        "Downregulated" = "#3498DB")) +
            theme_minimal() +
            theme(
                axis.text.y = element_text(size = 9),
                axis.title = element_text(size = 12, face = "bold"),
                legend.position = "bottom"
            ) +
            labs(
                title = paste("Top GSEA Results:", pathway_name),
                x = "Normalised Enrichment Score (NES)",
                y = "Pathway"
            ) +
            geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50")

        ggsave(output_file, plot = p, width = 10, height = 6)
    }

    # Run GSEA Analysis
    hallmark_pathways <- read_gmt("${hallmark_gmt}")
    hallmark_results <- run_gsea(gene_ranks, hallmark_pathways, "Hallmark")
    write.csv(hallmark_results, "gsea_hallmark_results.csv", row.names = FALSE)
    create_enrichment_plot(hallmark_results, "Hallmark Gene Sets", "gsea_hallmark_top20_plot.pdf")

    reactome_pathways <- read_gmt("${reactome_gmt}")
    reactome_results <- run_gsea(gene_ranks, reactome_pathways, "Reactome")
    write.csv(reactome_results, "gsea_reactome_results.csv", row.names = FALSE)
    create_enrichment_plot(reactome_results, "Reactome Pathways", "gsea_reactome_top20_plot.pdf")

    # --- GENERATE SUMMARY REPORT ---
    sink("gsea_combined_summary.txt")
    cat("=================================================================\\n")
    cat("     GENE SET ENRICHMENT ANALYSIS (GSEA) SUMMARY\\n")
    cat("=================================================================\\n\\n")

    cat("Analysis Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n")
    cat("Input genes:", length(gene_ranks), "\\n")
    cat("Thresholds: FDR (padj) < 0.05 and Nominal (pval) < 0.05\\n\\n")

    # --- HALLMARK SUMMARY ---
    cat("-----------------------------------------------------------------\\n")
    cat("HALLMARK GENE SETS\\n")
    cat("-----------------------------------------------------------------\\n")
    cat("Total gene sets tested:", nrow(hallmark_results), "\\n\\n")

    # 1. FDR Significant
    cat("[ FDR SIGNIFICANT (padj < 0.05) ]\\n")
    cat("Total:", sum(hallmark_results\$padj < 0.05, na.rm = TRUE), "\\n")
    cat("  - Enriched in condition 1 (NES>0):", sum(hallmark_results\$NES > 0 & hallmark_results\$padj < 0.05, na.rm = TRUE), "\\n")
    cat("  - Enriched in condition 2 (NES<0):", sum(hallmark_results\$NES < 0 & hallmark_results\$padj < 0.05, na.rm = TRUE), "\\n\\n")

    # 2. Nominal P-value Significant
    cat("[ NOMINALLY SIGNIFICANT (pval < 0.05) ]\\n")
    cat("Total:", sum(hallmark_results\$pval < 0.05, na.rm = TRUE), "\\n")
    cat("  - Enriched in condition 1 (NES>0):", sum(hallmark_results\$NES > 0 & hallmark_results\$pval < 0.05, na.rm = TRUE), "\\n")
    cat("  - Enriched in condition 2 (NES<0):", sum(hallmark_results\$NES < 0 & hallmark_results\$pval < 0.05, na.rm = TRUE), "\\n\\n")

    # Top 5 List
    sig_hallmark <- hallmark_results[hallmark_results\$padj < 0.05, ]
    if (nrow(sig_hallmark) > 0) {
        cat("Top 5 most significantly enriched Hallmark pathways (sorted by padj):\\n")
        top5 <- head(sig_hallmark[order(sig_hallmark\$padj), ], 5)
        for (i in 1:nrow(top5)) {
            cat(sprintf("  %d. %s\\n", i, top5\$pathway[i]))
            cat(sprintf("     NES = %.2f, pval = %.2e, padj = %.2e\\n",
                       top5\$NES[i], top5\$pval[i], top5\$padj[i]))
        }
    } else {
        cat("No pathways passed FDR < 0.05. Checking raw p-values...\\n")
        sig_hallmark_raw <- hallmark_results[hallmark_results\$pval < 0.05, ]
        if(nrow(sig_hallmark_raw) > 0) {
             cat("Top 5 (by raw p-value):\\n")
             top5 <- head(sig_hallmark_raw[order(sig_hallmark_raw\$pval), ], 5)
             for (i in 1:nrow(top5)) {
                cat(sprintf("  %d. %s\\n", i, top5\$pathway[i]))
                cat(sprintf("     NES = %.2f, pval = %.2e, padj = %.2e\\n",
                           top5\$NES[i], top5\$pval[i], top5\$padj[i]))
            }
        }
    }

    # --- REACTOME SUMMARY ---
    cat("\\n-----------------------------------------------------------------\\n")
    cat("REACTOME PATHWAYS\\n")
    cat("-----------------------------------------------------------------\\n")
    cat("Total pathways tested:", nrow(reactome_results), "\\n\\n")

    # 1. FDR Significant
    cat("[ FDR SIGNIFICANT (padj < 0.05) ]\\n")
    cat("Total:", sum(reactome_results\$padj < 0.05, na.rm = TRUE), "\\n")
    cat("  - Enriched in condition 1 (NES>0):", sum(reactome_results\$NES > 0 & reactome_results\$padj < 0.05, na.rm = TRUE), "\\n")
    cat("  - Enriched in condition 2 (NES<0):", sum(reactome_results\$NES < 0 & reactome_results\$padj < 0.05, na.rm = TRUE), "\\n\\n")

    # 2. Nominal P-value Significant
    cat("[ NOMINALLY SIGNIFICANT (pval < 0.05) ]\\n")
    cat("Total:", sum(reactome_results\$pval < 0.05, na.rm = TRUE), "\\n")
    cat("  - Enriched in condition 1 (NES>0):", sum(reactome_results\$NES > 0 & reactome_results\$pval < 0.05, na.rm = TRUE), "\\n")
    cat("  - Enriched in condition 2 (NES<0):", sum(reactome_results\$NES < 0 & reactome_results\$pval < 0.05, na.rm = TRUE), "\\n\\n")

    # Top 5 List
    sig_reactome <- reactome_results[reactome_results\$padj < 0.05, ]
    if (nrow(sig_reactome) > 0) {
        cat("Top 5 most significantly enriched Reactome pathways (sorted by padj):\\n")
        top5 <- head(sig_reactome[order(sig_reactome\$padj), ], 5)
        for (i in 1:nrow(top5)) {
            cat(sprintf("  %d. %s\\n", i, top5\$pathway[i]))
            cat(sprintf("     NES = %.2f, pval = %.2e, padj = %.2e\\n",
                       top5\$NES[i], top5\$pval[i], top5\$padj[i]))
        }
    } else {
        cat("No pathways passed FDR < 0.05. Checking raw p-values...\\n")
        sig_reactome_raw <- reactome_results[reactome_results\$pval < 0.05, ]
        if(nrow(sig_reactome_raw) > 0) {
             cat("Top 5 (by raw p-value):\\n")
             top5 <- head(sig_reactome_raw[order(sig_reactome_raw\$pval), ], 5)
             for (i in 1:nrow(top5)) {
                cat(sprintf("  %d. %s\\n", i, top5\$pathway[i]))
                cat(sprintf("     NES = %.2f, pval = %.2e, padj = %.2e\\n",
                           top5\$NES[i], top5\$pval[i], top5\$padj[i]))
            }
        }
    }

    cat("\\n=================================================================\\n")
    cat("END OF GSEA SUMMARY\\n")
    cat("=================================================================\\n")
    sink()

    cat("\\nGSEA analysis complete!\\n")
    """
}
