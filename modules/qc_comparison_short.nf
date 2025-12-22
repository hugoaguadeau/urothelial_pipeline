// Process: QC Comparison - Compare raw vs trimmed short read quality
process QC_COMPARISON_SHORT {
    container 'https://depot.galaxyproject.org/singularity/r-base:4.3.1'
    publishDir "${params.output_dir}/short_read/01_quality_control/comparisons", mode: 'copy'

    input:
    path raw_fastqc_summary
    path trimmed_fastqc_summary
    path raw_seqkit_summary
    path trimmed_seqkit_summary
    path fastp_summary

    output:
    path "qc_comparison.csv", emit: comparison
    path "qc_comparison_detailed.txt", emit: comparison_txt
    path "qc_comparison_summary.txt", emit: summary_txt

    script:
    """
    #!/usr/bin/env Rscript

    # Load all input files
    raw_fastqc <- read.csv("${raw_fastqc_summary}", stringsAsFactors = FALSE)
    trimmed_fastqc <- read.csv("${trimmed_fastqc_summary}", stringsAsFactors = FALSE)
    raw_seqkit <- read.csv("${raw_seqkit_summary}", stringsAsFactors = FALSE)
    trimmed_seqkit <- read.csv("${trimmed_seqkit_summary}", stringsAsFactors = FALSE)
    fastp <- read.csv("${fastp_summary}", stringsAsFactors = FALSE)

    # Extract base sample identifiers
    raw_seqkit\$base_id <- sub("-[BU]\$", "", raw_seqkit\$sample_id)
    trimmed_seqkit\$base_id <- sub("-[BU]\$", "", trimmed_seqkit\$sample_id)
    fastp\$base_id <- sub("-[BU]\$", "", fastp\$sample_id)
    raw_fastqc\$base_id <- sub("-[BU]\$", "", raw_fastqc\$sample_id)
    trimmed_fastqc\$base_id <- sub("-[BU]\$", "", trimmed_fastqc\$sample_id)
    
    # Get the unique base sample identifiers from trimmed data
    sample_ids <- unique(trimmed_seqkit\$base_id)

    # Create comparison data frame
    comparison <- data.frame(
        sample_id = character(),
        raw_total_reads = numeric(),
        trimmed_total_reads = numeric(),
        reads_removed = numeric(),
        reads_removed_percent = numeric(),
        reads_retained_percent = numeric(),
        raw_total_bases = numeric(),
        trimmed_total_bases = numeric(),
        bases_removed = numeric(),
        bases_removed_percent = numeric(),
        raw_mean_length = numeric(),
        trimmed_mean_length = numeric(),
        length_change = numeric(),
        raw_Q20_R1 = numeric(),
        raw_Q20_R2 = numeric(),
        raw_Q20_mean = numeric(),
        trimmed_Q20_R1 = numeric(),
        trimmed_Q20_R2 = numeric(),
        trimmed_Q20_mean = numeric(),
        Q20_improvement = numeric(),
        raw_Q30_R1 = numeric(),
        raw_Q30_R2 = numeric(),
        raw_Q30_mean = numeric(),
        trimmed_Q30_R1 = numeric(),
        trimmed_Q30_R2 = numeric(),
        trimmed_Q30_mean = numeric(),
        Q30_improvement = numeric(),
        raw_avgqual_R1 = numeric(),
        raw_avgqual_R2 = numeric(),
        raw_avgqual_mean = numeric(),
        trimmed_avgqual_R1 = numeric(),
        trimmed_avgqual_R2 = numeric(),
        trimmed_avgqual_mean = numeric(),
        avgqual_improvement = numeric(),
        raw_gc_mean = numeric(),
        trimmed_gc_mean = numeric(),
        gc_change = numeric(),
        raw_fastqc_passes = numeric(),
        trimmed_fastqc_passes = numeric(),
        fastqc_improvement = numeric(),
        stringsAsFactors = FALSE
    )

    # Process each sample
    for (sample_id in sample_ids) {
        
        # Retrieve data for this sample using base identifier
        raw_sk <- raw_seqkit[raw_seqkit\$base_id == sample_id, ]
        trim_sk <- trimmed_seqkit[trimmed_seqkit\$base_id == sample_id, ]
        
        # Verify data availability
        if (nrow(raw_sk) == 0 || nrow(trim_sk) == 0) {
            next
        }
        
        # Calculate read counts using R1 as total read pairs
        raw_reads <- as.numeric(gsub(",", "", raw_sk\$num_seqs[raw_sk\$read == "R1"]))
        trim_reads <- as.numeric(gsub(",", "", trim_sk\$num_seqs[trim_sk\$read == "R1"]))
        
        if (length(raw_reads) == 0 || length(trim_reads) == 0) {
            next
        }
        
        reads_removed <- raw_reads - trim_reads
        reads_removed_pct <- (reads_removed / raw_reads) * 100
        reads_retained_pct <- (trim_reads / raw_reads) * 100
        
        # Calculate base counts
        raw_bases <- sum(as.numeric(gsub(",", "", raw_sk\$sum_len)))
        trim_bases <- sum(as.numeric(gsub(",", "", trim_sk\$sum_len)))
        bases_removed <- raw_bases - trim_bases
        bases_removed_pct <- (bases_removed / raw_bases) * 100
        
        # Calculate read length
        raw_length <- mean(as.numeric(raw_sk\$avg_len))
        trim_length <- mean(as.numeric(trim_sk\$avg_len))
        length_change <- trim_length - raw_length
        
        # Extract quality metrics
        raw_q20_r1 <- as.numeric(raw_sk\$Q20_percent[raw_sk\$read == "R1"])
        raw_q20_r2 <- as.numeric(raw_sk\$Q20_percent[raw_sk\$read == "R2"])
        raw_q20_mean <- mean(c(raw_q20_r1, raw_q20_r2))
        
        trim_q20_r1 <- as.numeric(trim_sk\$Q20_percent[trim_sk\$read == "R1"])
        trim_q20_r2 <- as.numeric(trim_sk\$Q20_percent[trim_sk\$read == "R2"])
        trim_q20_mean <- mean(c(trim_q20_r1, trim_q20_r2))
        
        q20_improvement <- trim_q20_mean - raw_q20_mean
        
        raw_q30_r1 <- as.numeric(raw_sk\$Q30_percent[raw_sk\$read == "R1"])
        raw_q30_r2 <- as.numeric(raw_sk\$Q30_percent[raw_sk\$read == "R2"])
        raw_q30_mean <- mean(c(raw_q30_r1, raw_q30_r2))
        
        trim_q30_r1 <- as.numeric(trim_sk\$Q30_percent[trim_sk\$read == "R1"])
        trim_q30_r2 <- as.numeric(trim_sk\$Q30_percent[trim_sk\$read == "R2"])
        trim_q30_mean <- mean(c(trim_q30_r1, trim_q30_r2))
        
        q30_improvement <- trim_q30_mean - raw_q30_mean
        
        raw_qual_r1 <- as.numeric(raw_sk\$avg_qual[raw_sk\$read == "R1"])
        raw_qual_r2 <- as.numeric(raw_sk\$avg_qual[raw_sk\$read == "R2"])
        raw_qual_mean <- mean(c(raw_qual_r1, raw_qual_r2))
        
        trim_qual_r1 <- as.numeric(trim_sk\$avg_qual[trim_sk\$read == "R1"])
        trim_qual_r2 <- as.numeric(trim_sk\$avg_qual[trim_sk\$read == "R2"])
        trim_qual_mean <- mean(c(trim_qual_r1, trim_qual_r2))
        
        qual_improvement <- trim_qual_mean - raw_qual_mean
        
        # Calculate GC content
        raw_gc <- mean(as.numeric(raw_sk\$gc_percent))
        trim_gc <- mean(as.numeric(trim_sk\$gc_percent))
        gc_change <- trim_gc - raw_gc
        
        # Calculate FastQC pass rates
        raw_fqc <- raw_fastqc[raw_fastqc\$base_id == sample_id, ]
        trim_fqc <- trimmed_fastqc[trimmed_fastqc\$base_id == sample_id, ]
        
        # Count passes across quality checks
        qc_columns <- c("per_base_quality", "per_sequence_quality", 
                       "per_base_content", "duplication_levels", "adapter_content")
        
        raw_passes <- 0
        trim_passes <- 0
        
        for (col in qc_columns) {
            if (col %in% colnames(raw_fqc)) {
                raw_passes <- raw_passes + sum(raw_fqc[[col]] == "pass", na.rm = TRUE)
            }
            if (col %in% colnames(trim_fqc)) {
                trim_passes <- trim_passes + sum(trim_fqc[[col]] == "pass", na.rm = TRUE)
            }
        }
        
        fastqc_improvement <- trim_passes - raw_passes
        
        # Construct row for this sample
        row <- data.frame(
            sample_id = sample_id,
            raw_total_reads = raw_reads,
            trimmed_total_reads = trim_reads,
            reads_removed = reads_removed,
            reads_removed_percent = reads_removed_pct,
            reads_retained_percent = reads_retained_pct,
            raw_total_bases = raw_bases,
            trimmed_total_bases = trim_bases,
            bases_removed = bases_removed,
            bases_removed_percent = bases_removed_pct,
            raw_mean_length = raw_length,
            trimmed_mean_length = trim_length,
            length_change = length_change,
            raw_Q20_R1 = raw_q20_r1,
            raw_Q20_R2 = raw_q20_r2,
            raw_Q20_mean = raw_q20_mean,
            trimmed_Q20_R1 = trim_q20_r1,
            trimmed_Q20_R2 = trim_q20_r2,
            trimmed_Q20_mean = trim_q20_mean,
            Q20_improvement = q20_improvement,
            raw_Q30_R1 = raw_q30_r1,
            raw_Q30_R2 = raw_q30_r2,
            raw_Q30_mean = raw_q30_mean,
            trimmed_Q30_R1 = trim_q30_r1,
            trimmed_Q30_R2 = trim_q30_r2,
            trimmed_Q30_mean = trim_q30_mean,
            Q30_improvement = q30_improvement,
            raw_avgqual_R1 = raw_qual_r1,
            raw_avgqual_R2 = raw_qual_r2,
            raw_avgqual_mean = raw_qual_mean,
            trimmed_avgqual_R1 = trim_qual_r1,
            trimmed_avgqual_R2 = trim_qual_r2,
            trimmed_avgqual_mean = trim_qual_mean,
            avgqual_improvement = qual_improvement,
            raw_gc_mean = raw_gc,
            trimmed_gc_mean = trim_gc,
            gc_change = gc_change,
            raw_fastqc_passes = raw_passes,
            trimmed_fastqc_passes = trim_passes,
            fastqc_improvement = fastqc_improvement,
            stringsAsFactors = FALSE
        )
        
        comparison <- rbind(comparison, row)
    }

    # Sort by sample identifier
    comparison <- comparison[order(comparison\$sample_id), ]
    rownames(comparison) <- NULL

    # Write CSV output
    write.csv(comparison, "qc_comparison.csv", row.names = FALSE)

    # Generate detailed text report
    sink("qc_comparison_detailed.txt")
    cat("================================================================================\\n")
    cat("         QC COMPARISON: RAW vs TRIMMED SHORT READS\\n")
    cat("================================================================================\\n\\n")
    cat("Detailed comparison of read quality before and after Fastp trimming\\n")
    cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n\\n")

    for (i in 1:nrow(comparison)) {
        cat("================================================================================\\n")
        cat(sprintf("SAMPLE: %s\\n", comparison\$sample_id[i]))
        cat("================================================================================\\n\\n")
        
        cat("READ COUNTS:\\n")
        cat("--------------------------------------------------------------------------------\\n")
        cat(sprintf("  Raw total reads:          %s\\n",
                   format(comparison\$raw_total_reads[i], big.mark = ",")))
        cat(sprintf("  Trimmed total reads:      %s\\n",
                   format(comparison\$trimmed_total_reads[i], big.mark = ",")))
        cat(sprintf("  Reads removed:            %s (%.2f%%)\\n",
                   format(comparison\$reads_removed[i], big.mark = ","),
                   comparison\$reads_removed_percent[i]))
        cat(sprintf("  Reads retained:           %.2f%%\\n",
                   comparison\$reads_retained_percent[i]))
        
        cat("\\nBASE COUNTS:\\n")
        cat("--------------------------------------------------------------------------------\\n")
        cat(sprintf("  Raw total bases:          %s\\n",
                   format(comparison\$raw_total_bases[i], big.mark = ",")))
        cat(sprintf("  Trimmed total bases:      %s\\n",
                   format(comparison\$trimmed_total_bases[i], big.mark = ",")))
        cat(sprintf("  Bases removed:            %s (%.2f%%)\\n",
                   format(comparison\$bases_removed[i], big.mark = ","),
                   comparison\$bases_removed_percent[i]))
        cat(sprintf("  Mean read length change:  %.1f bp\\n",
                   comparison\$length_change[i]))
        
        cat("\\nQUALITY METRICS (Q30%%):\\n")
        cat("--------------------------------------------------------------------------------\\n")
        cat(sprintf("  Raw Q30%%:     R1=%.2f%%, R2=%.2f%%, Mean=%.2f%%\\n",
                   comparison\$raw_Q30_R1[i], comparison\$raw_Q30_R2[i], 
                   comparison\$raw_Q30_mean[i]))
        cat(sprintf("  Trimmed Q30%%: R1=%.2f%%, R2=%.2f%%, Mean=%.2f%%\\n",
                   comparison\$trimmed_Q30_R1[i], comparison\$trimmed_Q30_R2[i], 
                   comparison\$trimmed_Q30_mean[i]))
        cat(sprintf("  Q30 improvement:          %+.2f%%\\n", 
                   comparison\$Q30_improvement[i]))
        
        cat("\\nQUALITY METRICS (Q20%%):\\n")
        cat("--------------------------------------------------------------------------------\\n")
        cat(sprintf("  Raw Q20%%:     R1=%.2f%%, R2=%.2f%%, Mean=%.2f%%\\n",
                   comparison\$raw_Q20_R1[i], comparison\$raw_Q20_R2[i], 
                   comparison\$raw_Q20_mean[i]))
        cat(sprintf("  Trimmed Q20%%: R1=%.2f%%, R2=%.2f%%, Mean=%.2f%%\\n",
                   comparison\$trimmed_Q20_R1[i], comparison\$trimmed_Q20_R2[i], 
                   comparison\$trimmed_Q20_mean[i]))
        cat(sprintf("  Q20 improvement:          %+.2f%%\\n", 
                   comparison\$Q20_improvement[i]))
        
        cat("\\nAVERAGE QUALITY SCORE:\\n")
        cat("--------------------------------------------------------------------------------\\n")
        cat(sprintf("  Raw avg qual:   R1=%.2f, R2=%.2f, Mean=%.2f\\n",
                   comparison\$raw_avgqual_R1[i], comparison\$raw_avgqual_R2[i], 
                   comparison\$raw_avgqual_mean[i]))
        cat(sprintf("  Trimmed avg qual: R1=%.2f, R2=%.2f, Mean=%.2f\\n",
                   comparison\$trimmed_avgqual_R1[i], comparison\$trimmed_avgqual_R2[i], 
                   comparison\$trimmed_avgqual_mean[i]))
        cat(sprintf("  Quality improvement:      %+.2f\\n", 
                   comparison\$avgqual_improvement[i]))
        
        cat("\\nGC CONTENT:\\n")
        cat("--------------------------------------------------------------------------------\\n")
        cat(sprintf("  Raw GC%%:                  %.2f%%\\n", comparison\$raw_gc_mean[i]))
        cat(sprintf("  Trimmed GC%%:              %.2f%%\\n", comparison\$trimmed_gc_mean[i]))
        cat(sprintf("  GC change:                %+.2f%%\\n", comparison\$gc_change[i]))
        
        cat("\\nFASTQC QUALITY CHECKS:\\n")
        cat("--------------------------------------------------------------------------------\\n")
        cat(sprintf("  Raw passes:               %d/10\\n", comparison\$raw_fastqc_passes[i]))
        cat(sprintf("  Trimmed passes:           %d/10\\n", comparison\$trimmed_fastqc_passes[i]))
        cat(sprintf("  Improvement:              %+d checks\\n",
                   comparison\$fastqc_improvement[i]))
        
        cat("\\n")
    }
    
    cat("================================================================================\\n")
    sink()

    # Generate summary report
    sink("qc_comparison_summary.txt")
    cat("================================================================================\\n")
    cat("    QC COMPARISON SUMMARY: PROCESSING EFFECTIVENESS (SHORT READS)\\n")
    cat("================================================================================\\n\\n")
    cat("Overall statistics for Fastp processing\\n")
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
    cat(sprintf("  Raw mean Q30%%:            %.2f%%\\n", mean(comparison\$raw_Q30_mean)))
    cat(sprintf("  Trimmed mean Q30%%:        %.2f%%\\n", mean(comparison\$trimmed_Q30_mean)))
    cat(sprintf("  Raw mean avg quality:     %.2f\\n", mean(comparison\$raw_avgqual_mean)))
    cat(sprintf("  Trimmed mean avg quality: %.2f\\n", mean(comparison\$trimmed_avgqual_mean)))
    
    cat("\\nFASTQC CHECKS:\\n")
    cat("--------------------------------------------------------------------------------\\n")
    cat(sprintf("  Mean raw passes:          %.1f/10\\n", 
               mean(comparison\$raw_fastqc_passes)))
    cat(sprintf("  Mean trimmed passes:      %.1f/10\\n", 
               mean(comparison\$trimmed_fastqc_passes)))
    cat(sprintf("  Mean improvement:         %+.1f checks\\n", 
               mean(comparison\$fastqc_improvement)))
    
    cat("\\n================================================================================\\n")
    sink()
    """
}