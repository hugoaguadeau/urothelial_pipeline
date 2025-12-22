// Process: SeqKit Statistical Analysis - Comprehensive stats on QC data
process SEQKIT_STATS_ANALYSIS_SHORT {
    container 'https://depot.galaxyproject.org/singularity/r-base:4.3.1'
    publishDir "${params.output_dir}/short_read/01_quality_control/raw/summaries", mode: 'copy'

    input:
    path seqkit_summary  // The summary CSV we created earlier
    path sample_info     // samples.csv to get tissue groups

    output:
    path "sample_statistics.csv", emit: sample_stats
    path "sample_statistics.txt", emit: sample_stats_txt
    path "read_pair_comparison.csv", emit: pair_comparison
    path "read_pair_comparison.txt", emit: pair_comparison_txt
    path "group_statistics.csv", emit: group_stats
    path "group_statistics.txt", emit: group_stats_txt

    script:
    """
    #!/usr/bin/env Rscript

    # Read the seqkit summary
    seqkit <- read.csv("${seqkit_summary}", stringsAsFactors = FALSE)
    
    # Read sample information
    samples <- read.csv("${sample_info}", stringsAsFactors = FALSE)
    
    # Filter to short read samples only
    samples <- samples[samples\$seq_type == "short", ]
    
    # ========================================================================
    # OUTPUT 1: PER-SAMPLE DESCRIPTIVE STATISTICS
    # ========================================================================
    
    # Get unique sample IDs
    sample_ids <- unique(seqkit\$sample_id)
    
    # Create data frame for sample statistics
    sample_stats <- data.frame(
        sample_id = character(),
        tissue = character(),
        total_read_pairs = numeric(),
        total_bases = numeric(),
        mean_read_length = numeric(),
        mean_Q20_percent = numeric(),
        mean_Q30_percent = numeric(),
        mean_avg_qual = numeric(),
        mean_gc_percent = numeric(),
        Q20_range = numeric(),
        Q30_range = numeric(),
        avgqual_range = numeric(),
        gc_range = numeric(),
        Q20_sd = numeric(),
        Q30_sd = numeric(),
        avgqual_sd = numeric(),
        gc_sd = numeric(),
        stringsAsFactors = FALSE
    )
    
    # Calculate statistics for each sample
    for (sample_id in sample_ids) {
        
        # Get R1 and R2 data for this sample
        sample_data <- seqkit[seqkit\$sample_id == sample_id, ]
        
        # Extract base sample ID without tissue suffix for matching
        # Y2391-U -> Y2391, Y2631-B -> Y2631
        base_sample_id <- sub("-[BU]\$", "", sample_id)
        
        # Get tissue type from sample info using base ID
        # Use the "tissue" column directly (bladder/ureter)
        tissue <- samples\$tissue[samples\$sample_id == base_sample_id]
        if (length(tissue) == 0) tissue <- NA
        else tissue <- tissue[1]
        
        # Convert character numbers with commas to numeric
        num_seqs <- as.numeric(gsub(",", "", sample_data\$num_seqs))
        sum_len <- as.numeric(gsub(",", "", sample_data\$sum_len))
        avg_len <- as.numeric(sample_data\$avg_len)
        Q20 <- as.numeric(sample_data\$Q20_percent)
        Q30 <- as.numeric(sample_data\$Q30_percent)
        avgqual <- as.numeric(sample_data\$avg_qual)
        gc <- as.numeric(sample_data\$gc_percent)
        
        # Calculate statistics
        row <- data.frame(
            sample_id = sample_id,
            tissue = tissue,
            total_read_pairs = num_seqs[1],  # R1 and R2 have same count
            total_bases = sum(sum_len),
            mean_read_length = mean(avg_len),
            mean_Q20_percent = mean(Q20),
            mean_Q30_percent = mean(Q30),
            mean_avg_qual = mean(avgqual),
            mean_gc_percent = mean(gc),
            Q20_range = max(Q20) - min(Q20),
            Q30_range = max(Q30) - min(Q30),
            avgqual_range = max(avgqual) - min(avgqual),
            gc_range = max(gc) - min(gc),
            Q20_sd = sd(Q20),
            Q30_sd = sd(Q30),
            avgqual_sd = sd(avgqual),
            gc_sd = sd(gc),
            stringsAsFactors = FALSE
        )
        
        sample_stats <- rbind(sample_stats, row)
    }
    
    # Sort by sample_id
    sample_stats <- sample_stats[order(sample_stats\$sample_id), ]
    rownames(sample_stats) <- NULL
    
    # Write CSV output
    write.csv(sample_stats, "sample_statistics.csv", row.names = FALSE)
    
    # Write human-readable text report
    sink("sample_statistics.txt")
    cat("================================================================================\\n")
    cat("                    PER-SAMPLE DESCRIPTIVE STATISTICS\\n")
    cat("================================================================================\\n\\n")
    cat("Analysis of quality metrics averaged across R1 and R2 reads\\n")
    cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n\\n")
    
    for (i in 1:nrow(sample_stats)) {
        cat("--------------------------------------------------------------------------------\\n")
        cat(sprintf("SAMPLE: %s (%s)\\n", sample_stats\$sample_id[i], 
                   toupper(sample_stats\$tissue[i])))
        cat("--------------------------------------------------------------------------------\\n")
        cat(sprintf("Total read pairs:     %s\\n", 
                   format(sample_stats\$total_read_pairs[i], big.mark = ",")))
        cat(sprintf("Total bases:          %s\\n", 
                   format(sample_stats\$total_bases[i], big.mark = ",")))
        cat(sprintf("Mean read length:     %.0f bp\\n", 
                   sample_stats\$mean_read_length[i]))
        cat("\\nQuality Metrics (mean ± SD):\\n")
        cat(sprintf("  Q20%%:               %.2f%% ± %.2f%%\\n", 
                   sample_stats\$mean_Q20_percent[i], sample_stats\$Q20_sd[i]))
        cat(sprintf("  Q30%%:               %.2f%% ± %.2f%%\\n", 
                   sample_stats\$mean_Q30_percent[i], sample_stats\$Q30_sd[i]))
        cat(sprintf("  Average quality:    %.2f ± %.2f\\n", 
                   sample_stats\$mean_avg_qual[i], sample_stats\$avgqual_sd[i]))
        cat(sprintf("  GC content:         %.2f%% ± %.2f%%\\n", 
                   sample_stats\$mean_gc_percent[i], sample_stats\$gc_sd[i]))
        cat("\\nR1-R2 Variability (range):\\n")
        cat(sprintf("  Q20%% range:         %.2f%%\\n", sample_stats\$Q20_range[i]))
        cat(sprintf("  Q30%% range:         %.2f%%\\n", sample_stats\$Q30_range[i]))
        cat(sprintf("  Quality range:      %.2f\\n", sample_stats\$avgqual_range[i]))
        cat(sprintf("  GC range:           %.2f%%\\n", sample_stats\$gc_range[i]))
        cat("\\n")
    }
    
    cat("================================================================================\\n")
    cat("Summary across all samples:\\n")
    cat(sprintf("  Total samples:        %d\\n", nrow(sample_stats)))
    cat(sprintf("  Bladder samples:      %d\\n", 
               sum(sample_stats\$tissue == "bladder", na.rm = TRUE)))
    cat(sprintf("  Ureter samples:       %d\\n", 
               sum(sample_stats\$tissue == "ureter", na.rm = TRUE)))
    cat(sprintf("  Mean Q30%%:            %.2f%% (range: %.2f%% - %.2f%%)\\n",
               mean(sample_stats\$mean_Q30_percent),
               min(sample_stats\$mean_Q30_percent),
               max(sample_stats\$mean_Q30_percent)))
    cat(sprintf("  Mean GC content:      %.2f%% (range: %.2f%% - %.2f%%)\\n",
               mean(sample_stats\$mean_gc_percent),
               min(sample_stats\$mean_gc_percent),
               max(sample_stats\$mean_gc_percent)))
    cat("================================================================================\\n")
    sink()
    
    # ========================================================================
    # OUTPUT 2: READ PAIR COMPARISON (R1 vs R2)
    # ========================================================================
    
    pair_comparison <- data.frame(
        sample_id = character(),
        tissue = character(),
        Q20_R1 = numeric(),
        Q20_R2 = numeric(),
        Q20_difference = numeric(),
        Q20_percent_drop = numeric(),
        Q30_R1 = numeric(),
        Q30_R2 = numeric(),
        Q30_difference = numeric(),
        Q30_percent_drop = numeric(),
        avgqual_R1 = numeric(),
        avgqual_R2 = numeric(),
        avgqual_difference = numeric(),
        avgqual_percent_drop = numeric(),
        gc_R1 = numeric(),
        gc_R2 = numeric(),
        gc_difference = numeric(),
        stringsAsFactors = FALSE
    )
    
    # Calculate pair comparisons for each sample
    for (sample_id in sample_ids) {
        
        sample_data <- seqkit[seqkit\$sample_id == sample_id, ]
        
        # Get R1 and R2 rows
        r1 <- sample_data[sample_data\$read == "R1", ]
        r2 <- sample_data[sample_data\$read == "R2", ]
        
        # Extract base sample ID for matching
        base_sample_id <- sub("-[BU]\$", "", sample_id)
        
        # Get tissue type from sample info - use "tissue" column
        tissue <- samples\$tissue[samples\$sample_id == base_sample_id]
        if (length(tissue) == 0) tissue <- NA
        else tissue <- tissue[1]
        
        # Extract values
        Q20_R1 <- as.numeric(r1\$Q20_percent)
        Q20_R2 <- as.numeric(r2\$Q20_percent)
        Q30_R1 <- as.numeric(r1\$Q30_percent)
        Q30_R2 <- as.numeric(r2\$Q30_percent)
        avgqual_R1 <- as.numeric(r1\$avg_qual)
        avgqual_R2 <- as.numeric(r2\$avg_qual)
        gc_R1 <- as.numeric(r1\$gc_percent)
        gc_R2 <- as.numeric(r2\$gc_percent)
        
        # Calculate differences and percent drops
        row <- data.frame(
            sample_id = sample_id,
            tissue = tissue,
            Q20_R1 = Q20_R1,
            Q20_R2 = Q20_R2,
            Q20_difference = Q20_R2 - Q20_R1,
            Q20_percent_drop = ((Q20_R2 - Q20_R1) / Q20_R1) * 100,
            Q30_R1 = Q30_R1,
            Q30_R2 = Q30_R2,
            Q30_difference = Q30_R2 - Q30_R1,
            Q30_percent_drop = ((Q30_R2 - Q30_R1) / Q30_R1) * 100,
            avgqual_R1 = avgqual_R1,
            avgqual_R2 = avgqual_R2,
            avgqual_difference = avgqual_R2 - avgqual_R1,
            avgqual_percent_drop = ((avgqual_R2 - avgqual_R1) / avgqual_R1) * 100,
            gc_R1 = gc_R1,
            gc_R2 = gc_R2,
            gc_difference = gc_R2 - gc_R1,
            stringsAsFactors = FALSE
        )
        
        pair_comparison <- rbind(pair_comparison, row)
    }
    
    # Sort by sample_id
    pair_comparison <- pair_comparison[order(pair_comparison\$sample_id), ]
    rownames(pair_comparison) <- NULL
    
    # Write CSV output
    write.csv(pair_comparison, "read_pair_comparison.csv", row.names = FALSE)
    
    # Write human-readable text report
    sink("read_pair_comparison.txt")
    cat("================================================================================\\n")
    cat("                      READ PAIR COMPARISON (R1 vs R2)\\n")
    cat("================================================================================\\n\\n")
    cat("Analysis of quality degradation from R1 to R2\\n")
    cat("Note: Some quality drop in R2 is normal for Illumina sequencing\\n")
    cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n\\n")
    
    for (i in 1:nrow(pair_comparison)) {
        cat("--------------------------------------------------------------------------------\\n")
        cat(sprintf("SAMPLE: %s (%s)\\n", pair_comparison\$sample_id[i], 
                   toupper(pair_comparison\$tissue[i])))
        cat("--------------------------------------------------------------------------------\\n")
        cat("Q20%% (bases with quality ≥20):\\n")
        cat(sprintf("  R1: %.2f%%  |  R2: %.2f%%  |  Drop: %.2f%% (%.1f%% relative)\\n",
                   pair_comparison\$Q20_R1[i], pair_comparison\$Q20_R2[i],
                   pair_comparison\$Q20_difference[i], pair_comparison\$Q20_percent_drop[i]))
        
        cat("\\nQ30%% (bases with quality ≥30):\\n")
        cat(sprintf("  R1: %.2f%%  |  R2: %.2f%%  |  Drop: %.2f%% (%.1f%% relative)\\n",
                   pair_comparison\$Q30_R1[i], pair_comparison\$Q30_R2[i],
                   pair_comparison\$Q30_difference[i], pair_comparison\$Q30_percent_drop[i]))
        
        cat("\\nAverage Quality Score:\\n")
        cat(sprintf("  R1: %.2f  |  R2: %.2f  |  Drop: %.2f (%.1f%% relative)\\n",
                   pair_comparison\$avgqual_R1[i], pair_comparison\$avgqual_R2[i],
                   pair_comparison\$avgqual_difference[i], pair_comparison\$avgqual_percent_drop[i]))
        
        cat("\\nGC Content:\\n")
        cat(sprintf("  R1: %.2f%%  |  R2: %.2f%%  |  Difference: %.2f%%\\n",
                   pair_comparison\$gc_R1[i], pair_comparison\$gc_R2[i],
                   pair_comparison\$gc_difference[i]))
        cat("\\n")
    }
    
    cat("================================================================================\\n")
    cat("Summary of R2 degradation:\\n")
    cat(sprintf("  Mean Q30%% drop:       %.2f%%\\n", 
               mean(pair_comparison\$Q30_difference)))
    cat(sprintf("  Mean relative Q30%% drop: %.1f%%\\n", 
               mean(pair_comparison\$Q30_percent_drop)))
    cat(sprintf("  Largest Q30%% drop:    %.2f%% (%s)\\n",
               min(pair_comparison\$Q30_difference),
               pair_comparison\$sample_id[which.min(pair_comparison\$Q30_difference)]))
    cat(sprintf("  Smallest Q30%% drop:   %.2f%% (%s)\\n",
               max(pair_comparison\$Q30_difference),
               pair_comparison\$sample_id[which.max(pair_comparison\$Q30_difference)]))
    cat("================================================================================\\n")
    sink()
    
    # ========================================================================
    # OUTPUT 3: GROUP STATISTICS (Bladder vs Ureter)
    # ========================================================================
    
    # Split samples by tissue type using the "tissue" column directly
    bladder_samples <- sample_stats[sample_stats\$tissue == "bladder", ]
    ureter_samples <- sample_stats[sample_stats\$tissue == "ureter", ]
    
    # Check if we have enough samples for statistical tests
    if (nrow(bladder_samples) < 2 || nrow(ureter_samples) < 2) {
        # Not enough samples for t-test, create empty output with message
        group_stats <- data.frame(
            note = "Not enough samples in one or both groups for statistical testing",
            bladder_n = nrow(bladder_samples),
            ureter_n = nrow(ureter_samples),
            stringsAsFactors = FALSE
        )
        
        write.csv(group_stats, "group_statistics.csv", row.names = FALSE)
        
        sink("group_statistics.txt")
        cat("================================================================================\\n")
        cat("                   GROUP COMPARISON: BLADDER vs URETER\\n")
        cat("================================================================================\\n\\n")
        cat("ERROR: Not enough samples for statistical testing\\n")
        cat(sprintf("  Bladder samples: %d\\n", nrow(bladder_samples)))
        cat(sprintf("  Ureter samples:  %d\\n", nrow(ureter_samples)))
        cat("\\nAt least 2 samples per group are required for t-test and Wilcoxon test.\\n")
        cat("================================================================================\\n")
        sink()
        
    } else {
        # Metrics to compare
        metrics <- c("total_read_pairs", "mean_Q20_percent", "mean_Q30_percent", 
                     "mean_avg_qual", "mean_gc_percent")
        
        group_stats <- data.frame(
            metric = character(),
            bladder_n = integer(),
            bladder_mean = numeric(),
            bladder_sd = numeric(),
            bladder_min = numeric(),
            bladder_max = numeric(),
            bladder_range = numeric(),
            ureter_n = integer(),
            ureter_mean = numeric(),
            ureter_sd = numeric(),
            ureter_min = numeric(),
            ureter_max = numeric(),
            ureter_range = numeric(),
            difference = numeric(),
            ttest_pvalue = numeric(),
            ttest_significant = character(),
            wilcoxon_pvalue = numeric(),
            wilcoxon_significant = character(),
            stringsAsFactors = FALSE
        )
        
        # Calculate statistics for each metric
        for (metric in metrics) {
            
            bladder_vals <- bladder_samples[[metric]]
            ureter_vals <- ureter_samples[[metric]]
            
            # T-test
            ttest_result <- t.test(bladder_vals, ureter_vals)
            ttest_pval <- ttest_result\$p.value
            ttest_sig <- ifelse(ttest_pval < 0.05, "YES", "NO")
            
            # Wilcoxon test
            wilcox_result <- wilcox.test(bladder_vals, ureter_vals, exact = FALSE)
            wilcox_pval <- wilcox_result\$p.value
            wilcox_sig <- ifelse(wilcox_pval < 0.05, "YES", "NO")
            
            # Create row
            row <- data.frame(
                metric = metric,
                bladder_n = length(bladder_vals),
                bladder_mean = mean(bladder_vals),
                bladder_sd = sd(bladder_vals),
                bladder_min = min(bladder_vals),
                bladder_max = max(bladder_vals),
                bladder_range = max(bladder_vals) - min(bladder_vals),
                ureter_n = length(ureter_vals),
                ureter_mean = mean(ureter_vals),
                ureter_sd = sd(ureter_vals),
                ureter_min = min(ureter_vals),
                ureter_max = max(ureter_vals),
                ureter_range = max(ureter_vals) - min(ureter_vals),
                difference = mean(bladder_vals) - mean(ureter_vals),
                ttest_pvalue = ttest_pval,
                ttest_significant = ttest_sig,
                wilcoxon_pvalue = wilcox_pval,
                wilcoxon_significant = wilcox_sig,
                stringsAsFactors = FALSE
            )
            
            group_stats <- rbind(group_stats, row)
        }
        
        # Write CSV output
        write.csv(group_stats, "group_statistics.csv", row.names = FALSE)
        
        # Write human-readable text report
        sink("group_statistics.txt")
        cat("================================================================================\\n")
        cat("                   GROUP COMPARISON: BLADDER vs URETER\\n")
        cat("================================================================================\\n\\n")
        cat("Statistical comparison of quality metrics between tissue types\\n")
        cat("Tests performed: Student's t-test and Wilcoxon rank-sum test\\n")
        cat("Significance threshold: p < 0.05\\n")
        cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n\\n")
        
        cat(sprintf("Sample sizes: Bladder n=%d, Ureter n=%d\\n\\n", 
                   nrow(bladder_samples), nrow(ureter_samples)))
        
        for (i in 1:nrow(group_stats)) {
            cat("--------------------------------------------------------------------------------\\n")
            
            # Format metric name nicely
            metric_name <- group_stats\$metric[i]
            metric_display <- switch(metric_name,
                "total_read_pairs" = "Total Read Pairs",
                "mean_Q20_percent" = "Mean Q20%",
                "mean_Q30_percent" = "Mean Q30%",
                "mean_avg_qual" = "Mean Average Quality",
                "mean_gc_percent" = "Mean GC Content",
                metric_name)
            
            cat(sprintf("%s:\\n", metric_display))
            cat("--------------------------------------------------------------------------------\\n")
            
            # Bladder stats
            cat("BLADDER:\\n")
            if (metric_name == "total_read_pairs") {
                cat(sprintf("  Mean: %s\\n", format(round(group_stats\$bladder_mean[i]), big.mark = ",")))
                cat(sprintf("  SD:   %s\\n", format(round(group_stats\$bladder_sd[i]), big.mark = ",")))
                cat(sprintf("  Range: %s - %s\\n", 
                           format(round(group_stats\$bladder_min[i]), big.mark = ","),
                           format(round(group_stats\$bladder_max[i]), big.mark = ",")))
            } else {
                cat(sprintf("  Mean: %.2f\\n", group_stats\$bladder_mean[i]))
                cat(sprintf("  SD:   %.2f\\n", group_stats\$bladder_sd[i]))
                cat(sprintf("  Range: %.2f - %.2f\\n", 
                           group_stats\$bladder_min[i], group_stats\$bladder_max[i]))
            }
            
            cat("\\nURETER:\\n")
            if (metric_name == "total_read_pairs") {
                cat(sprintf("  Mean: %s\\n", format(round(group_stats\$ureter_mean[i]), big.mark = ",")))
                cat(sprintf("  SD:   %s\\n", format(round(group_stats\$ureter_sd[i]), big.mark = ",")))
                cat(sprintf("  Range: %s - %s\\n", 
                           format(round(group_stats\$ureter_min[i]), big.mark = ","),
                           format(round(group_stats\$ureter_max[i]), big.mark = ",")))
            } else {
                cat(sprintf("  Mean: %.2f\\n", group_stats\$ureter_mean[i]))
                cat(sprintf("  SD:   %.2f\\n", group_stats\$ureter_sd[i]))
                cat(sprintf("  Range: %.2f - %.2f\\n", 
                           group_stats\$ureter_min[i], group_stats\$ureter_max[i]))
            }
            
            # Difference
            cat("\\nDIFFERENCE (Bladder - Ureter):\\n")
            if (metric_name == "total_read_pairs") {
                cat(sprintf("  %s\\n", format(round(group_stats\$difference[i]), big.mark = ",")))
            } else {
                cat(sprintf("  %.2f\\n", group_stats\$difference[i]))
            }
            
            # Statistical tests
            cat("\\nSTATISTICAL TESTS:\\n")
            cat(sprintf("  T-test p-value:        %.4f  %s\\n", 
                       group_stats\$ttest_pvalue[i],
                       ifelse(group_stats\$ttest_significant[i] == "YES", 
                             "[SIGNIFICANT]", "[not significant]")))
            cat(sprintf("  Wilcoxon p-value:      %.4f  %s\\n", 
                       group_stats\$wilcoxon_pvalue[i],
                       ifelse(group_stats\$wilcoxon_significant[i] == "YES", 
                             "[SIGNIFICANT]", "[not significant]")))
            cat("\\n")
        }
        
        cat("================================================================================\\n")
        cat("SUMMARY:\\n")
        cat("--------------------------------------------------------------------------------\\n")
        sig_count <- sum(group_stats\$ttest_significant == "YES" | 
                        group_stats\$wilcoxon_significant == "YES")
        cat(sprintf("Metrics tested: %d\\n", nrow(group_stats)))
        cat(sprintf("Significant differences found: %d\\n", sig_count))
        
        if (sig_count > 0) {
            cat("\\nSignificant metrics:\\n")
            for (i in 1:nrow(group_stats)) {
                if (group_stats\$ttest_significant[i] == "YES" || 
                    group_stats\$wilcoxon_significant[i] == "YES") {
                    metric_name <- group_stats\$metric[i]
                    metric_display <- switch(metric_name,
                        "total_read_pairs" = "Total Read Pairs",
                        "mean_Q20_percent" = "Mean Q20%",
                        "mean_Q30_percent" = "Mean Q30%",
                        "mean_avg_qual" = "Mean Average Quality",
                        "mean_gc_percent" = "Mean GC Content",
                        metric_name)
                    cat(sprintf("  - %s\\n", metric_display))
                }
            }
        }
        cat("================================================================================\\n")
        sink()
    }
    """
}
