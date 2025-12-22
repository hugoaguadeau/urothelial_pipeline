// Process: StringTie Summary - Analyse transcript quantification results
process STRINGTIE_SUMMARY_SHORT {
    container 'https://depot.galaxyproject.org/singularity/r-base:4.3.1'
    publishDir "${params.output_dir}/short_read/04_quantification/summaries", mode: 'copy'

    input:
    path gene_abundance_files
    path transcript_abundance_files
    path gene_count_matrix
    path transcript_count_matrix

    output:
    path "stringtie_summary.csv", emit: summary
    path "stringtie_summary.txt", emit: summary_txt
    path "stringtie_detailed_stats.csv", emit: detailed_stats

    script:
    """
    #!/usr/bin/env Rscript

    # ============================================================================
    # StringTie Quantification Summary
    # ============================================================================

    # Read gene count matrix
    gene_counts <- read.csv("${gene_count_matrix}", row.names = 1, check.names = FALSE)
    
    # Read transcript count matrix
    transcript_counts <- read.csv("${transcript_count_matrix}", row.names = 1, check.names = FALSE)

    # Get sample names
    sample_ids <- colnames(gene_counts)

    # Initialise results data frame
    results <- data.frame(
        sample_id = character(),
        total_genes_detected = numeric(),
        genes_low_expression = numeric(),
        genes_medium_expression = numeric(),
        genes_high_expression = numeric(),
        total_transcripts_detected = numeric(),
        mean_transcripts_per_gene = numeric(),
        median_gene_count = numeric(),
        median_transcript_count = numeric(),
        percent_genes_expressed = numeric(),
        stringsAsFactors = FALSE
    )

    # Process each sample
    for (sample_id in sample_ids) {
        # Gene-level statistics
        gene_sample <- gene_counts[[sample_id]]
        genes_detected <- sum(gene_sample > 0)
        
        # Expression level categories (arbitrary thresholds for educational purposes)
        # Low: 1-10 counts, Medium: 11-100 counts, High: >100 counts
        genes_low <- sum(gene_sample > 0 & gene_sample <= 10)
        genes_medium <- sum(gene_sample > 10 & gene_sample <= 100)
        genes_high <- sum(gene_sample > 100)
        
        # Transcript-level statistics
        transcript_sample <- transcript_counts[[sample_id]]
        transcripts_detected <- sum(transcript_sample > 0)
        
        # Calculate transcripts per gene (approximation)
        # This is a rough estimate as we don't have direct gene-transcript mapping here
        mean_transcripts_per_gene <- transcripts_detected / max(genes_detected, 1)
        
        # Median counts
        median_gene <- median(gene_sample[gene_sample > 0])
        if (is.na(median_gene)) median_gene <- 0
        
        median_transcript <- median(transcript_sample[transcript_sample > 0])
        if (is.na(median_transcript)) median_transcript <- 0
        
        # Percentage of genes expressed (out of total genes in count matrix)
        percent_expressed <- (genes_detected / nrow(gene_counts)) * 100
        
        # Add to results
        row <- data.frame(
            sample_id = sample_id,
            total_genes_detected = genes_detected,
            genes_low_expression = genes_low,
            genes_medium_expression = genes_medium,
            genes_high_expression = genes_high,
            total_transcripts_detected = transcripts_detected,
            mean_transcripts_per_gene = mean_transcripts_per_gene,
            median_gene_count = median_gene,
            median_transcript_count = median_transcript,
            percent_genes_expressed = percent_expressed,
            stringsAsFactors = FALSE
        )
        
        results <- rbind(results, row)
    }

    # Sort results by sample ID
    results <- results[order(results\$sample_id), ]
    rownames(results) <- NULL

    # Write CSV output
    write.csv(results, "stringtie_summary.csv", row.names = FALSE)

    # Create detailed statistics for comparison
    detailed_stats <- data.frame(
        metric = c(
            "Mean genes detected",
            "Range of genes detected",
            "Mean transcripts detected",
            "Range of transcripts detected",
            "Mean transcripts per gene",
            "Mean percentage genes expressed"
        ),
        value = c(
            sprintf("%.0f", mean(results\$total_genes_detected)),
            sprintf("%.0f - %.0f", min(results\$total_genes_detected), max(results\$total_genes_detected)),
            sprintf("%.0f", mean(results\$total_transcripts_detected)),
            sprintf("%.0f - %.0f", min(results\$total_transcripts_detected), max(results\$total_transcripts_detected)),
            sprintf("%.2f", mean(results\$mean_transcripts_per_gene)),
            sprintf("%.2f%%", mean(results\$percent_genes_expressed))
        ),
        stringsAsFactors = FALSE
    )
    
    write.csv(detailed_stats, "stringtie_detailed_stats.csv", row.names = FALSE)

    # ============================================================================
    # Generate human-readable text summary
    # ============================================================================
    
    sink("stringtie_summary.txt")
    
    cat("================================================================================\\n")
    cat("                   STRINGTIE QUANTIFICATION SUMMARY\\n")
    cat("================================================================================\\n\\n")
    cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n")
    cat("================================================================================\\n\\n")

    # Per-sample summaries
    for (i in 1:nrow(results)) {
        cat("--------------------------------------------------------------------------------\\n")
        cat(sprintf("SAMPLE: %s\\n", results\$sample_id[i]))
        cat("--------------------------------------------------------------------------------\\n\\n")
        
        cat("GENE-LEVEL QUANTIFICATION:\\n")
        cat(sprintf("  Total genes detected:           %s\\n", 
                   format(results\$total_genes_detected[i], big.mark = ",")))
        cat(sprintf("  Percentage of genes expressed:  %.2f%%\\n", 
                   results\$percent_genes_expressed[i]))
        cat(sprintf("  Median gene count:              %s\\n\\n",
                   format(round(results\$median_gene_count[i]), big.mark = ",")))
        
        cat("  Expression level distribution:\\n")
        cat(sprintf("    Low expression (1-10 counts):    %s (%.1f%%)\\n",
                   format(results\$genes_low_expression[i], big.mark = ","),
                   (results\$genes_low_expression[i]/results\$total_genes_detected[i])*100))
        cat(sprintf("    Medium expression (11-100):      %s (%.1f%%)\\n",
                   format(results\$genes_medium_expression[i], big.mark = ","),
                   (results\$genes_medium_expression[i]/results\$total_genes_detected[i])*100))
        cat(sprintf("    High expression (>100 counts):   %s (%.1f%%)\\n\\n",
                   format(results\$genes_high_expression[i], big.mark = ","),
                   (results\$genes_high_expression[i]/results\$total_genes_detected[i])*100))
        
        cat("TRANSCRIPT-LEVEL QUANTIFICATION:\\n")
        cat(sprintf("  Total transcripts detected:     %s\\n",
                   format(results\$total_transcripts_detected[i], big.mark = ",")))
        cat(sprintf("  Mean transcripts per gene:      %.2f\\n",
                   results\$mean_transcripts_per_gene[i]))
        cat(sprintf("  Median transcript count:        %s\\n\\n",
                   format(round(results\$median_transcript_count[i]), big.mark = ",")))
    }

    # Overall summary
    cat("================================================================================\\n")
    cat("OVERALL COHORT SUMMARY\\n")
    cat("================================================================================\\n\\n")
    
    cat(sprintf("Total samples analysed:           %d\\n\\n", nrow(results)))
    
    cat("GENE DETECTION:\\n")
    cat(sprintf("  Mean genes detected:            %s\\n",
               format(round(mean(results\$total_genes_detected)), big.mark = ",")))
    cat(sprintf("  Standard deviation:             %s\\n",
               format(round(sd(results\$total_genes_detected)), big.mark = ",")))
    cat(sprintf("  Range:                          %s - %s\\n",
               format(min(results\$total_genes_detected), big.mark = ","),
               format(max(results\$total_genes_detected), big.mark = ",")))
    cat(sprintf("  Mean percentage expressed:      %.2f%%\\n",
               mean(results\$percent_genes_expressed)))
    cat(sprintf("  Median percentage expressed:    %.2f%%\\n\\n",
               median(results\$percent_genes_expressed)))
    
    cat("TRANSCRIPT DETECTION:\\n")
    cat(sprintf("  Mean transcripts detected:      %s\\n",
               format(round(mean(results\$total_transcripts_detected)), big.mark = ",")))
    cat(sprintf("  Standard deviation:             %s\\n",
               format(round(sd(results\$total_transcripts_detected)), big.mark = ",")))
    cat(sprintf("  Range:                          %s - %s\\n",
               format(min(results\$total_transcripts_detected), big.mark = ","),
               format(max(results\$total_transcripts_detected), big.mark = ",")))
    cat(sprintf("  Mean transcripts per gene:      %.2f\\n",
               mean(results\$mean_transcripts_per_gene)))
    cat(sprintf("  Median transcripts per gene:    %.2f\\n\\n",
               median(results\$mean_transcripts_per_gene)))
    
    cat("EXPRESSION LEVEL PATTERNS:\\n")
    cat(sprintf("  Mean low expression genes:      %s (%.1f%%)\\n",
               format(round(mean(results\$genes_low_expression)), big.mark = ","),
               mean(results\$genes_low_expression/results\$total_genes_detected)*100))
    cat(sprintf("  Mean medium expression genes:   %s (%.1f%%)\\n",
               format(round(mean(results\$genes_medium_expression)), big.mark = ","),
               mean(results\$genes_medium_expression/results\$total_genes_detected)*100))
    cat(sprintf("  Mean high expression genes:     %s (%.1f%%)\\n",
               format(round(mean(results\$genes_high_expression)), big.mark = ","),
               mean(results\$genes_high_expression/results\$total_genes_detected)*100))
    
    cat("================================================================================\\n")
    
    sink()
    """
}

