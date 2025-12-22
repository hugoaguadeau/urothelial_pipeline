// Process: QC Comparison - Compare raw vs processed long read quality
process QC_COMPARISON_LONG {
    container 'https://depot.galaxyproject.org/singularity/r-base:4.3.1'
    publishDir "${params.output_dir}/long_read/01_quality_control/comparisons", mode: 'copy'

    input:
    path raw_fastqc_summary
    path trimmed_fastqc_summary
    path raw_seqkit_summary
    path trimmed_seqkit_summary
    path fastplong_summary

    output:
    path "qc_comparison.csv", emit: comparison
    path "qc_comparison_detailed.txt", emit: comparison_txt
    path "qc_comparison_summary.txt", emit: summary_txt

    script:
    """
    #!/usr/bin/env Rscript

    # Read all input files
    raw_fastqc <- read.csv("${raw_fastqc_summary}", stringsAsFactors = FALSE)
    trimmed_fastqc <- read.csv("${trimmed_fastqc_summary}", stringsAsFactors = FALSE)
    raw_seqkit <- read.csv("${raw_seqkit_summary}", stringsAsFactors = FALSE)
    trimmed_seqkit <- read.csv("${trimmed_seqkit_summary}", stringsAsFactors = FALSE)
    fastplong <- read.csv("${fastplong_summary}", stringsAsFactors = FALSE)

    # Get unique sample IDs
    sample_ids <- unique(raw_seqkit\$sample_id)

    # Create comprehensive comparison data frame
    comparison <- data.frame(
        sample_id = character(),
        raw_total_reads = numeric(),
        processed_total_reads = numeric(),
        reads_removed = numeric(),
        reads_removed_percent = numeric(),
        reads_retained_percent = numeric(),
        raw_total_bases = numeric(),
        processed_total_bases = numeric(),
        bases_removed = numeric(),
        bases_removed_percent = numeric(),
        raw_mean_length = numeric(),
        processed_mean_length = numeric(),
        length_change = numeric(),
        raw_Q20_percent = numeric(),
        processed_Q20_percent = numeric(),
        Q20_improvement = numeric(),
        raw_Q30_percent = numeric(),
        processed_Q30_percent = numeric(),
        Q30_improvement = numeric(),
        raw_avgqual = numeric(),
        processed_avgqual = numeric(),
        avgqual_improvement = numeric(),
        raw_gc_mean = numeric(),
        processed_gc_mean = numeric(),
        gc_change = numeric(),
        raw_fastqc_passes = numeric(),
        processed_fastqc_passes = numeric(),
        fastqc_improvement = numeric(),
        stringsAsFactors = FALSE
    )

    # Process each sample
    for (sample_id in sample_ids) {

        # Get data for this sample
        raw_sk <- raw_seqkit[raw_seqkit\$sample_id == sample_id, ]
        proc_sk <- trimmed_seqkit[trimmed_seqkit\$sample_id == sample_id, ]

        # Check if we have data
        if (nrow(raw_sk) == 0 || nrow(proc_sk) == 0) {
            next
        }

        # Read counts
        raw_reads <- as.numeric(gsub(",", "", raw_sk\$num_seqs))
        proc_reads <- as.numeric(gsub(",", "", proc_sk\$num_seqs))

        if (length(raw_reads) == 0 || length(proc_reads) == 0) {
            next
        }

        reads_removed <- raw_reads - proc_reads
        reads_removed_pct <- (reads_removed / raw_reads) * 100
        reads_retained_pct <- (proc_reads / raw_reads) * 100

        # Base counts
        raw_bases <- as.numeric(gsub(",", "", raw_sk\$sum_len))
        proc_bases <- as.numeric(gsub(",", "", proc_sk\$sum_len))
        bases_removed <- raw_bases - proc_bases
        bases_removed_pct <- (bases_removed / raw_bases) * 100

        # Read length (remove commas before converting)
        raw_length <- as.numeric(gsub(",", "", raw_sk\$avg_len))
        proc_length <- as.numeric(gsub(",", "", proc_sk\$avg_len))
        length_change <- proc_length - raw_length

        # Quality metrics
        raw_q20 <- as.numeric(raw_sk\$Q20_percent)
        proc_q20 <- as.numeric(proc_sk\$Q20_percent)
        q20_improvement <- proc_q20 - raw_q20

        raw_q30 <- as.numeric(raw_sk\$Q30_percent)
        proc_q30 <- as.numeric(proc_sk\$Q30_percent)
        q30_improvement <- proc_q30 - raw_q30

        raw_qual <- as.numeric(raw_sk\$avg_qual)
        proc_qual <- as.numeric(proc_sk\$avg_qual)
        qual_improvement <- proc_qual - raw_qual

        # GC content
        raw_gc <- as.numeric(raw_sk\$gc_percent)
        proc_gc <- as.numeric(proc_sk\$gc_percent)
        gc_change <- proc_gc - raw_gc

        # FastQC pass rates
        raw_fqc <- raw_fastqc[raw_fastqc\$sample_id == sample_id, ]
        proc_fqc <- trimmed_fastqc[trimmed_fastqc\$sample_id == sample_id, ]

        # Count passes
        qc_columns <- c("per_base_quality", "per_sequence_quality",
                       "per_base_content", "duplication_levels", "adapter_content")

        raw_passes <- 0
        proc_passes <- 0

        for (col in qc_columns) {
            if (col %in% colnames(raw_fqc)) {
                raw_passes <- raw_passes + sum(raw_fqc[[col]] == "pass", na.rm = TRUE)
            }
            if (col %in% colnames(proc_fqc)) {
                proc_passes <- proc_passes + sum(proc_fqc[[col]] == "pass", na.rm = TRUE)
            }
        }

        fastqc_improvement <- proc_passes - raw_passes

        # Create row
        row <- data.frame(
            sample_id = sample_id,
            raw_total_reads = raw_reads,
            processed_total_reads = proc_reads,
            reads_removed = reads_removed,
            reads_removed_percent = reads_removed_pct,
            reads_retained_percent = reads_retained_pct,
            raw_total_bases = raw_bases,
            processed_total_bases = proc_bases,
            bases_removed = bases_removed,
            bases_removed_percent = bases_removed_pct,
            raw_mean_length = raw_length,
            processed_mean_length = proc_length,
            length_change = length_change,
            raw_Q20_percent = raw_q20,
            processed_Q20_percent = proc_q20,
            Q20_improvement = q20_improvement,
            raw_Q30_percent = raw_q30,
            processed_Q30_percent = proc_q30,
            Q30_improvement = q30_improvement,
            raw_avgqual = raw_qual,
            processed_avgqual = proc_qual,
            avgqual_improvement = qual_improvement,
            raw_gc_mean = raw_gc,
            processed_gc_mean = proc_gc,
            gc_change = gc_change,
            raw_fastqc_passes = raw_passes,
            processed_fastqc_passes = proc_passes,
            fastqc_improvement = fastqc_improvement,
            stringsAsFactors = FALSE
        )

        comparison <- rbind(comparison, row)
    }

    # Sort by sample_id
    comparison <- comparison[order(comparison\$sample_id), ]
    rownames(comparison) <- NULL

    # Write CSV
    write.csv(comparison, "qc_comparison.csv", row.names = FALSE)

    # Detailed text report
    sink("qc_comparison_detailed.txt")
    cat("================================================================================\\n")
    cat("         QC COMPARISON: RAW vs PROCESSED LONG READS\\n")
    cat("================================================================================\\n\\n")
    cat("Detailed comparison of read quality before and after Pychopper + Fastplong\\n")
    cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n\\n")

    for (i in 1:nrow(comparison)) {
        cat("================================================================================\\n")
        cat(sprintf("SAMPLE: %s\\n", comparison\$sample_id[i]))
        cat("================================================================================\\n\\n")

        cat("READ COUNTS:\\n")
        cat("--------------------------------------------------------------------------------\\n")
        cat(sprintf("  Raw total reads:          %s\\n",
                   format(comparison\$raw_total_reads[i], big.mark = ",")))
        cat(sprintf("  Processed total reads:    %s\\n",
                   format(comparison\$processed_total_reads[i], big.mark = ",")))
        cat(sprintf("  Reads removed:            %s (%.2f%%)\\n",
                   format(comparison\$reads_removed[i], big.mark = ","),
                   comparison\$reads_removed_percent[i]))
        cat(sprintf("  Reads retained:           %.2f%%\\n",
                   comparison\$reads_retained_percent[i]))

        cat("\\nBASE COUNTS:\\n")
        cat("--------------------------------------------------------------------------------\\n")
        cat(sprintf("  Raw total bases:          %s\\n",
                   format(comparison\$raw_total_bases[i], big.mark = ",")))
        cat(sprintf("  Processed total bases:    %s\\n",
                   format(comparison\$processed_total_bases[i], big.mark = ",")))
        cat(sprintf("  Bases removed:            %s (%.2f%%)\\n",
                   format(comparison\$bases_removed[i], big.mark = ","),
                   comparison\$bases_removed_percent[i]))
        cat(sprintf("  Mean read length change:  %.1f bp\\n",
                   comparison\$length_change[i]))

        cat("\\nQUALITY METRICS (Q30%%):\\n")
        cat("--------------------------------------------------------------------------------\\n")
        cat(sprintf("  Raw Q30%%:                 %.2f%%\\n",
                   comparison\$raw_Q30_percent[i]))
        cat(sprintf("  Processed Q30%%:           %.2f%%\\n",
                   comparison\$processed_Q30_percent[i]))
        cat(sprintf("  Q30 improvement:          %+.2f%%\\n",
                   comparison\$Q30_improvement[i]))

        cat("\\nQUALITY METRICS (Q20%%):\\n")
        cat("--------------------------------------------------------------------------------\\n")
        cat(sprintf("  Raw Q20%%:                 %.2f%%\\n",
                   comparison\$raw_Q20_percent[i]))
        cat(sprintf("  Processed Q20%%:           %.2f%%\\n",
                   comparison\$processed_Q20_percent[i]))
        cat(sprintf("  Q20 improvement:          %+.2f%%\\n",
                   comparison\$Q20_improvement[i]))

        cat("\\nAVERAGE QUALITY SCORE:\\n")
        cat("--------------------------------------------------------------------------------\\n")
        cat(sprintf("  Raw avg qual:             %.2f\\n",
                   comparison\$raw_avgqual[i]))
        cat(sprintf("  Processed avg qual:       %.2f\\n",
                   comparison\$processed_avgqual[i]))
        cat(sprintf("  Quality improvement:      %+.2f\\n",
                   comparison\$avgqual_improvement[i]))

        cat("\\nGC CONTENT:\\n")
        cat("--------------------------------------------------------------------------------\\n")
        cat(sprintf("  Raw GC%%:                  %.2f%%\\n", comparison\$raw_gc_mean[i]))
        cat(sprintf("  Processed GC%%:            %.2f%%\\n", comparison\$processed_gc_mean[i]))
        cat(sprintf("  GC change:                %+.2f%%\\n", comparison\$gc_change[i]))

        cat("\\nFASTQC QUALITY CHECKS:\\n")
        cat("--------------------------------------------------------------------------------\\n")
        cat(sprintf("  Raw passes:               %d/5\\n", comparison\$raw_fastqc_passes[i]))
        cat(sprintf("  Processed passes:         %d/5\\n", comparison\$processed_fastqc_passes[i]))
        cat(sprintf("  Improvement:              %+d checks\\n",
                   comparison\$fastqc_improvement[i]))

        cat("\\n")
    }

    cat("================================================================================\\n")
    sink()

    # Summary report
    sink("qc_comparison_summary.txt")
    cat("================================================================================\\n")
    cat("    QC COMPARISON SUMMARY: PROCESSING EFFECTIVENESS (LONG READS)\\n")
    cat("================================================================================\\n\\n")
    cat("Overall statistics for Pychopper + Fastplong processing\\n")
    cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n\\n")

    cat("OVERALL STATISTICS ACROSS ALL SAMPLES:\\n")
    cat("================================================================================\\n\\n")

    cat("READ RETENTION:\\n")
    cat("--------------------------------------------------------------------------------\\n")
    cat(sprintf("  Mean reads retained:      %.2f%%\\n",
               mean(comparison\$reads_retained_percent)))
    cat(sprintf("  Range:                    %.2f%% - %.2f%%\\n",
               min(comparison\$reads_retained_percent),
               max(comparison\$reads_retained_percent)))
    cat(sprintf("  Mean reads removed:       %.2f%%\\n",
               mean(comparison\$reads_removed_percent)))
    cat(sprintf("  Total reads removed:      %s\\n",
               format(sum(comparison\$reads_removed), big.mark = ",")))

    cat("\\nBASE RETENTION:\\n")
    cat("--------------------------------------------------------------------------------\\n")
    cat(sprintf("  Mean bases removed:       %.2f%%\\n",
               mean(comparison\$bases_removed_percent)))
    cat(sprintf("  Total bases removed:      %s\\n",
               format(sum(comparison\$bases_removed), big.mark = ",")))
    cat(sprintf("  Mean length change:       %.1f bp\\n",
               mean(comparison\$length_change)))

    cat("\\nQUALITY IMPROVEMENTS:\\n")
    cat("--------------------------------------------------------------------------------\\n")
    cat(sprintf("  Mean Q30%% improvement:    %+.2f%%\\n",
               mean(comparison\$Q30_improvement)))
    cat(sprintf("  Mean Q20%% improvement:    %+.2f%%\\n",
               mean(comparison\$Q20_improvement)))
    cat(sprintf("  Mean avgqual improvement: %+.2f\\n",
               mean(comparison\$avgqual_improvement)))

    cat("\\nQUALITY BEFORE AND AFTER:\\n")
    cat("--------------------------------------------------------------------------------\\n")
    cat(sprintf("  Raw mean Q30%%:            %.2f%%\\n", mean(comparison\$raw_Q30_percent)))
    cat(sprintf("  Processed mean Q30%%:      %.2f%%\\n", mean(comparison\$processed_Q30_percent)))
    cat(sprintf("  Raw mean avg quality:     %.2f\\n", mean(comparison\$raw_avgqual)))
    cat(sprintf("  Processed mean avg qual:  %.2f\\n", mean(comparison\$processed_avgqual)))

    cat("\\nFASTQC CHECKS:\\n")
    cat("--------------------------------------------------------------------------------\\n")
    cat(sprintf("  Mean raw passes:          %.1f/5\\n",
               mean(comparison\$raw_fastqc_passes)))
    cat(sprintf("  Mean processed passes:    %.1f/5\\n",
               mean(comparison\$processed_fastqc_passes)))
    cat(sprintf("  Mean improvement:         %+.1f checks\\n",
               mean(comparison\$fastqc_improvement)))

    cat("\\n================================================================================\\n")
    sink()
    """
}
