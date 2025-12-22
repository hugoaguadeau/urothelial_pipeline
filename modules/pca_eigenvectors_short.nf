// Process: PCA Eigenvectors Analysis (Short Reads)
process PCA_EIGENVECTORS_SHORT {
    container 'https://depot.galaxyproject.org/singularity/mulled-v2-9e1a016dbabab74e6ff1620b7d53af6985ff2247:76760ec587e610aa88f11f1f45ff16abd8789105-0'
    publishDir "${params.output_dir}/short_read/05_differential_expression/pca_analysis", mode: 'copy'

    input:
    path normalised_counts    // From DESeq2
    path sample_info          // samples.csv

    output:
    path "pca_loadings.csv", emit: pca_loadings
    path "pca_explained_variance.csv", emit: pca_var
    path "top_contributing_genes.csv", emit: top_genes
    path "*.pdf", emit: plots_pdf
    path "*.png", emit: plots_png

    script:
    """
    #!/usr/bin/env Rscript

    # ========================================================================
    # PCA EIGENVECTORS ANALYSIS - SHORT READS
    # ========================================================================
   
    # Load required packages
    library(ggplot2)
    library(ggrepel)

    # ========================================================================
    # SECTION 1: READ AND PREPARE DATA
    # ========================================================================

    cat("Reading normalised counts from DESeq2...\\n")
    norm_counts <- read.csv("${normalised_counts}", row.names = 1, check.names = FALSE)
    cat("Normalised counts dimensions:", nrow(norm_counts), "genes x",
        ncol(norm_counts), "samples\\n")

    # Read sample information
    cat("\\nReading sample information...\\n")
    sample_info_full <- read.csv("${sample_info}", stringsAsFactors = FALSE)

    # Normalise column names to lowercase to handle RIN/rin/Rin
    colnames(sample_info_full) <- tolower(colnames(sample_info_full))

    # Filter to short read samples only (avoid duplicate sample IDs)
    sample_info <- sample_info_full[sample_info_full\$seq_type == "short", ]
    cat("Filtered to", nrow(sample_info), "short read samples\\n")

    # Create a tissue mapping for each sample
    tissue_map <- setNames(sample_info\$tissue, sample_info\$sample_id)

    # Check if RIN column exists
    has_rin <- "rin" %in% colnames(sample_info)
    if(has_rin) {
        cat("RIN column found. Will generate RIN-based PCA plots.\\n")
        rin_map <- setNames(sample_info\$rin, sample_info\$sample_id)
    } else {
        cat("No 'rin' column found in samples.csv. Skipping RIN plots.\\n")
    }

    # ========================================================================
    # SECTION 2: PREPARE DATA FOR PCA
    # ========================================================================

    cat("\\nPreparing data for PCA...\\n")
    norm_counts_t <- t(norm_counts)

    # Remove genes with zero or near-zero variance
    gene_vars <- apply(norm_counts_t, 2, var)
    genes_to_keep <- gene_vars > 1e-8
    norm_counts_t <- norm_counts_t[, genes_to_keep]

    # Check we have enough data
    if (nrow(norm_counts_t) < 2) stop("ERROR: Need at least 2 samples for PCA")
    if (ncol(norm_counts_t) < 2) stop("ERROR: Need at least 2 genes with variance for PCA")

    # ========================================================================
    # SECTION 3: PERFORM PCA
    # ========================================================================

    cat("\\nPerforming PCA...\\n")
    pca_result <- prcomp(norm_counts_t, scale. = TRUE)
    pca_loadings <- pca_result\$rotation
    pca_summary <- summary(pca_result)
    pca_variance <- pca_summary\$importance
    pca_scores <- pca_result\$x

    pc1_var <- round(pca_variance[2, 1] * 100, 2)
    pc2_var <- round(pca_variance[2, 2] * 100, 2)

    # ========================================================================
    # SECTION 4: SAVE RAW RESULTS
    # ========================================================================

    write.csv(pca_loadings, "pca_loadings.csv")
    write.csv(pca_variance, "pca_explained_variance.csv")

    # ========================================================================
    # SECTION 5: IDENTIFY TOP CONTRIBUTING GENES
    # ========================================================================

    # Get top 20 genes for PC1
    pc1_loadings_abs <- abs(pca_loadings[, 1])
    pc1_top_indices <- order(pc1_loadings_abs, decreasing = TRUE)[1:20]
    pc1_top_genes <- rownames(pca_loadings)[pc1_top_indices]

    pc1_top_df <- data.frame(
        Gene = pc1_top_genes,
        Loading = pca_loadings[pc1_top_genes, 1],
        Abs_Loading = abs(pca_loadings[pc1_top_genes, 1]),
        PC = "PC1",
        Rank = 1:20,
        stringsAsFactors = FALSE
    )

    # Get top 20 genes for PC2
    pc2_loadings_abs <- abs(pca_loadings[, 2])
    pc2_top_indices <- order(pc2_loadings_abs, decreasing = TRUE)[1:20]
    pc2_top_genes <- rownames(pca_loadings)[pc2_top_indices]

    pc2_top_df <- data.frame(
        Gene = pc2_top_genes,
        Loading = pca_loadings[pc2_top_genes, 2],
        Abs_Loading = abs(pca_loadings[pc2_top_genes, 2]),
        PC = "PC2",
        Rank = 1:20,
        stringsAsFactors = FALSE
    )

    top_genes <- rbind(pc1_top_df, pc2_top_df)
    write.csv(top_genes, "top_contributing_genes.csv", row.names = FALSE)

    # ========================================================================
    # SECTION 6: PREPARE DATA FOR PLOTTING
    # ========================================================================

    pca_data <- as.data.frame(pca_scores)
    pca_data\$Sample <- rownames(pca_data)
    pca_data\$Tissue <- tissue_map[rownames(pca_data)]

    if(has_rin) {
        pca_data\$RIN <- rin_map[rownames(pca_data)]
    }

    # ========================================================================
    # SECTION 7: CREATE PLOTS 
    # ========================================================================

    cat("\\nCreating visualisations...\\n")

    # Define a 'Large' Publication Theme
    # Sizes are significantly increased to be legible when resized in papers
    pub_theme <- theme_bw() +
        theme(
            plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 24),
            axis.title = element_text(size = 24, face = "bold"),
            axis.text = element_text(size = 20, color = "black"),
            legend.title = element_text(size = 22, face = "bold"),
            legend.text = element_text(size = 20),
            legend.position = "top",
            strip.text = element_text(size = 22, face = "bold")
        )

    # Increased dimensions to allow space for large text
    save_both_formats <- function(plot_obj, filename_base, width = 12, height = 10) {
        pdf(paste0(filename_base, ".pdf"), width = width, height = height)
        print(plot_obj)
        dev.off()
        png(paste0(filename_base, ".png"), width = width * 100, height = height * 100, res = 100)
        print(plot_obj)
        dev.off()
    }

    # Plot 1: Basic PCA
    p_basic <- ggplot(pca_data, aes(x = PC1, y = PC2, colour = Tissue)) +
        geom_point(size = 8) +
        scale_colour_manual(values = c("bladder" = "#E41A1C", "ureter" = "#377EB8")) +
        labs(
            title = "Sample clustering by principal component analysis",
            x = paste0("Principal component 1 (", pc1_var, "%)"),
            y = paste0("Principal component 2 (", pc2_var, "%)"),
            colour = "Tissue"
        ) +
        pub_theme

    save_both_formats(p_basic, "pca_basic")

    # Plot 2: PCA with sample labels
    p_labelled <- ggplot(pca_data, aes(x = PC1, y = PC2, colour = Tissue)) +
        geom_point(size = 8) +
        geom_text_repel(
            aes(label = Sample),
            size = 8,
            box.padding = 1.0,
            point.padding = 0.5,
            force = 2,
            max.overlaps = 20,
            show.legend = FALSE
        ) +
        scale_colour_manual(values = c("bladder" = "#E41A1C", "ureter" = "#377EB8")) +
        labs(
            title = "Sample clustering by principal component analysis",
            x = paste0("Principal component 1 (", pc1_var, "%)"),
            y = paste0("Principal component 2 (", pc2_var, "%)"),
            colour = "Tissue"
        ) +
        pub_theme

    save_both_formats(p_labelled, "pca_labelled")

    # Plot 2b: RIN PCA (If data exists)
    if(has_rin) {
        p_rin <- ggplot(pca_data, aes(x = PC1, y = PC2, colour = RIN)) +
            geom_point(size = 9) +
            geom_text_repel(
                aes(label=Sample),
                color="black",
                size=8,
                box.padding = 1.0,
                point.padding = 0.5
            ) +
            scale_colour_gradient(low = "red", high = "blue", name = "RIN Score") +
            labs(
                title = "Principal component analysis coloured by RNA integrity",
                x = paste0("Principal component 1 (", pc1_var, "%)"),
                y = paste0("Principal component 2 (", pc2_var, "%)")
            ) +
            pub_theme

        save_both_formats(p_rin, "pca_by_rin_quality")
    }

    # Plot 3: Top contributing genes
    p_genes <- ggplot(top_genes, aes(x = reorder(Gene, Abs_Loading), y = Loading, fill = PC)) +
        geom_bar(stat = "identity") +
        facet_wrap(~ PC, scales = "free") +
        coord_flip() +
        scale_fill_manual(values = c("PC1" = "#66C2A5", "PC2" = "#FC8D62")) +
        labs(
            title = "Principal component gene loadings",
            subtitle = "Top 20 genes contributing to PC1 and PC2",
            x = "Gene",
            y = "Loading value"
        ) +
        pub_theme +
        theme(legend.position = "none")

    save_both_formats(p_genes, "top_contributing_genes", width = 14, height = 12)

    # Plot 4: Biplot
    gene_importance <- sqrt(pca_loadings[, 1]^2 + pca_loadings[, 2]^2)
    top_50_genes <- names(sort(gene_importance, decreasing = TRUE)[1:50])
    gene_coords <- as.data.frame(pca_loadings[top_50_genes, 1:2])
    gene_coords\$Gene <- rownames(gene_coords)

    scaling_factor <- 0.8 * max(abs(pca_data[, c("PC1", "PC2")])) /
                      max(abs(gene_coords[, c("PC1", "PC2")]))
    gene_coords\$PC1 <- gene_coords\$PC1 * scaling_factor
    gene_coords\$PC2 <- gene_coords\$PC2 * scaling_factor

    p_biplot <- ggplot() +
        geom_point(data = pca_data, aes(x = PC1, y = PC2, colour = Tissue), size = 8) +
        geom_text_repel(
            data = pca_data,
            aes(x = PC1, y = PC2, label = Sample, colour = Tissue),
            size = 8,
            box.padding = 1.0,
            max.overlaps = 20
        ) +
        geom_segment(
            data = gene_coords,
            aes(x = 0, y = 0, xend = PC1, yend = PC2),
            arrow = arrow(length = unit(0.3, "cm")),
            colour = "grey50", alpha = 0.6,
            size = 1
        ) +
        geom_text_repel(
            data = gene_coords,
            aes(x = PC1, y = PC2, label = Gene),
            size = 6,
            colour = "darkblue",
            max.overlaps = 30,
            box.padding = 0.5
        ) +
        scale_colour_manual(values = c("bladder" = "#E41A1C", "ureter" = "#377EB8")) +
        labs(
            title = "Principal component analysis biplot",
            subtitle = "Samples and top 50 genes with highest loadings",
            x = paste0("Principal component 1 (", pc1_var, "%)"),
            y = paste0("Principal component 2 (", pc2_var, "%)"),
            colour = "Tissue"
        ) +
        pub_theme

    save_both_formats(p_biplot, "biplot_samples_genes", width = 14, height = 12)

    cat("\\nPCA ANALYSIS COMPLETE\\n")
    """
}

