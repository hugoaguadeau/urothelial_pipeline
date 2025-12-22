// Process: SeqKit Statistical Analysis - Comprehensive stats on QC data (Long Reads)
process SEQKIT_STATS_ANALYSIS_LONG {
    container 'https://depot.galaxyproject.org/singularity/r-base:4.3.1'
    publishDir "${params.output_dir}/long_read/01_quality_control/raw/summaries", mode: 'copy'

    input:
    path seqkit_summary  // The summary CSV we created earlier
    path sample_info     // samples.csv to get tissue groups

    output:
    path "sample_statistics.csv", emit: sample_stats
    path "sample_statistics.txt", emit: sample_stats_txt
    path "group_statistics.csv", emit: group_stats
    path "group_statistics.txt", emit: group_stats_txt

    script:
    """
    #!/usr/bin/env Rscript

    # Read the seqkit summary
    seqkit <- read.csv("${seqkit_summary}", stringsAsFactors = FALSE)

    # Read sample information
    samples <- read.csv("${sample_info}", stringsAsFactors = FALSE)

    # Filter to long read samples only
    samples <- samples[samples\$seq_type == "long", ]

    # ========================================================================
    # OUTPUT 1: PER-SAMPLE DESCRIPTIVE STATISTICS
    # ========================================================================

    # Get unique sample IDs
    sample_ids <- unique(seqkit\$sample_id)

    # Create data frame for sample statistics
    sample_stats <- data.frame(
        sample_id = character(),
        tissue = character(),
        condition = character(),
        total_reads = numeric(),
        total_bases = numeric(),
        mean_read_length = numeric(),
        min_read_length = numeric(),
        max_read_length = numeric(),
        read_length_range = numeric(),
        Q20_percent = numeric(),
        Q30_percent = numeric(),
        mean_avg_qual = numeric(),
        gc_percent = numeric(),
        stringsAsFactors = FALSE
    )

    # Calculate statistics for each sample
    for (sample_id in sample_ids) {

        # Get data for this sample
        sample_data <- seqkit[seqkit\$sample_id == sample_id, ]

        # Get tissue type and condition from sample info
        tissue <- samples\$tissue[samples\$sample_id == sample_id]
        if (length(tissue) == 0) tissue <- NA
        else tissue <- tissue[1]

        condition <- samples\$condition[samples\$sample_id == sample_id]
        if (length(condition) == 0) condition <- NA
        else condition <- condition[1]

        # Convert character numbers with commas to numeric
        num_seqs <- as.numeric(gsub(",", "", sample_data\$num_seqs))
        sum_len <- as.numeric(gsub(",", "", sample_data\$sum_len))
        min_len <- as.numeric(gsub(",", "", sample_data\$min_len))
        avg_len <- as.numeric(sample_data\$avg_len)
        max_len <- as.numeric(gsub(",", "", sample_data\$max_len))
        Q20 <- as.numeric(sample_data\$Q20_percent)
        Q30 <- as.numeric(sample_data\$Q30_percent)
        avgqual <- as.numeric(sample_data\$avg_qual)
        gc <- as.numeric(sample_data\$gc_percent)

        # Calculate statistics
        row <- data.frame(
            sample_id = sample_id,
            tissue = tissue,
            condition = condition,
            total_reads = num_seqs,
            total_bases = sum_len,
            mean_read_length = avg_len,
            min_read_length = min_len,
            max_read_length = max_len,
            read_length_range = max_len - min_len,
            Q20_percent = Q20,
            Q30_percent = Q30,
            mean_avg_qual = avgqual,
            gc_percent = gc,
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
    cat("                              LONG READS (RAW)\\n")
    cat("================================================================================\\n\\n")
    cat("Analysis of quality metrics for Oxford Nanopore long-read sequencing data\\n")
    cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n\\n")

    for (i in 1:nrow(sample_stats)) {
        cat("--------------------------------------------------------------------------------\\n")
        cat(sprintf("SAMPLE: %s\\n", sample_stats\$sample_id[i]))
        cat("--------------------------------------------------------------------------------\\n")
        cat(sprintf("Tissue:               %s\\n", toupper(sample_stats\$tissue[i])))
        cat(sprintf("Condition:            %s\\n", sample_stats\$condition[i]))
        cat(sprintf("Total reads:          %s\\n",
                   format(sample_stats\$total_reads[i], big.mark = ",")))
        cat(sprintf("Total bases:          %s Gb\\n",
                   format(sample_stats\$total_bases[i] / 1e9, digits = 3)))
        cat("\\nRead Length Statistics:\\n")
        cat(sprintf("  Mean length:        %s bp\\n",
                   format(round(sample_stats\$mean_read_length[i]), big.mark = ",")))
        cat(sprintf("  Minimum length:     %s bp\\n",
                   format(sample_stats\$min_read_length[i], big.mark = ",")))
        cat(sprintf("  Maximum length:     %s bp\\n",
                   format(sample_stats\$max_read_length[i], big.mark = ",")))
        cat(sprintf("  Range:              %s bp\\n",
                   format(sample_stats\$read_length_range[i], big.mark = ",")))
        cat("\\nQuality Metrics:\\n")
        cat(sprintf("  Q20%%:               %.2f%%\\n", sample_stats\$Q20_percent[i]))
        cat(sprintf("  Q30%%:               %.2f%%\\n", sample_stats\$Q30_percent[i]))
        cat(sprintf("  Average quality:    %.2f\\n", sample_stats\$mean_avg_qual[i]))
        cat(sprintf("  GC content:         %.2f%%\\n", sample_stats\$gc_percent[i]))
        cat("\\n")
    }

    cat("================================================================================\\n")
    cat("SUMMARY ACROSS ALL SAMPLES:\\n")
    cat("--------------------------------------------------------------------------------\\n")
    cat(sprintf("Total samples:          %d\\n", nrow(sample_stats)))
    cat(sprintf("Bladder samples:        %d\\n",
               sum(sample_stats\$tissue == "bladder", na.rm = TRUE)))
    cat(sprintf("Ureter samples:         %d\\n",
               sum(sample_stats\$tissue == "ureter", na.rm = TRUE)))
    cat("\\n")
    cat(sprintf("Mean total reads:       %s (range: %s - %s)\\n",
               format(round(mean(sample_stats\$total_reads)), big.mark = ","),
               format(round(min(sample_stats\$total_reads)), big.mark = ","),
               format(round(max(sample_stats\$total_reads)), big.mark = ",")))
    cat(sprintf("Mean read length:       %s bp (range: %s - %s bp)\\n",
               format(round(mean(sample_stats\$mean_read_length)), big.mark = ","),
               format(round(min(sample_stats\$mean_read_length)), big.mark = ","),
               format(round(max(sample_stats\$mean_read_length)), big.mark = ",")))
    cat(sprintf("Mean Q30%%:              %.2f%% (range: %.2f%% - %.2f%%)\\n",
               mean(sample_stats\$Q30_percent),
               min(sample_stats\$Q30_percent),
               max(sample_stats\$Q30_percent)))
    cat(sprintf("Mean GC content:        %.2f%% (range: %.2f%% - %.2f%%)\\n",
               mean(sample_stats\$gc_percent),
               min(sample_stats\$gc_percent),
               max(sample_stats\$gc_percent)))
    cat("================================================================================\\n")
    sink()

    # ========================================================================
    # OUTPUT 2: GROUP STATISTICS (Bladder vs Ureter)
    # ========================================================================

    # Split samples by tissue type
    bladder_samples <- sample_stats[sample_stats\$tissue == "bladder", ]
    ureter_samples <- sample_stats[sample_stats\$tissue == "ureter", ]

    # Check if we have enough samples for statistical tests
    if (nrow(bladder_samples) < 2 || nrow(ureter_samples) < 2) {
        # Not enough samples for t-test, create output with message
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
        cat("                              LONG READS (RAW)\\n")
        cat("================================================================================\\n\\n")
        cat("ERROR: Not enough samples for statistical testing\\n")
        cat(sprintf("  Bladder samples: %d\\n", nrow(bladder_samples)))
        cat(sprintf("  Ureter samples:  %d\\n", nrow(ureter_samples)))
        cat("\\nAt least 2 samples per group are required for t-test and Wilcoxon test.\\n")
        cat("================================================================================\\n")
        sink()

    } else {
        # Metrics to compare
        metrics <- c("total_reads", "mean_read_length", "Q20_percent", 
                     "Q30_percent", "mean_avg_qual", "gc_percent")

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
        cat("                              LONG READS (RAW)\\n")
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
                "total_reads" = "Total Reads",
                "mean_read_length" = "Mean Read Length",
                "Q20_percent" = "Q20%",
                "Q30_percent" = "Q30%",
                "mean_avg_qual" = "Mean Average Quality",
                "gc_percent" = "GC Content",
                metric_name)

            cat(sprintf("%s:\\n", metric_display))
            cat("--------------------------------------------------------------------------------\\n")

            # Bladder stats
            cat("BLADDER:\\n")
            if (metric_name %in% c("total_reads", "mean_read_length")) {
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
            if (metric_name %in% c("total_reads", "mean_read_length")) {
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
            if (metric_name %in% c("total_reads", "mean_read_length")) {
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
                        "total_reads" = "Total Reads",
                        "mean_read_length" = "Mean Read Length",
                        "Q20_percent" = "Q20%",
                        "Q30_percent" = "Q30%",
                        "mean_avg_qual" = "Mean Average Quality",
                        "gc_percent" = "GC Content",
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
