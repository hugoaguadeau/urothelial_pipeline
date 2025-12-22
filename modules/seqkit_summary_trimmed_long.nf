// Process: SeqKit Summary - Combine trimmed read stats into one CSV (Long Reads)
process SEQKIT_SUMMARY_TRIMMED_LONG {
    container 'https://depot.galaxyproject.org/singularity/r-base:4.3.1'
    publishDir "${params.output_dir}/long_read/01_quality_control/trimmed/summaries", mode: 'copy'

    input:
    path seqkit_files  // All the individual seqkit stats files from trimmed reads

    output:
    path "trimmed_seqkit_summary.csv", emit: summary
    path "trimmed_seqkit_summary.txt", emit: summary_txt

    script:
    """
    #!/usr/bin/env Rscript

    # Get list of all seqkit stats files
    stats_files <- list.files(pattern = "*_seqkit_stats.txt", full.names = TRUE)

    # Create empty data frame to store results
    results <- data.frame(
        sample_id = character(),
        num_seqs = character(),
        sum_len = character(),
        min_len = character(),
        avg_len = character(),
        max_len = character(),
        Q20_percent = character(),
        Q30_percent = character(),
        avg_qual = character(),
        gc_percent = character(),
        stringsAsFactors = FALSE
    )

    # Process each stats file
    for (stats_file in stats_files) {

        # Read the file
        data <- read.table(stats_file, header = TRUE,
                          stringsAsFactors = FALSE,
                          comment.char = "",
                          check.names = FALSE)

        # Should have exactly 1 data row (plus header)
        if (nrow(data) < 1) {
            next
        }

        # Parse filename to get sample_id
        # Format: Y2391_seqkit_stats.txt -> Y2391
        sample_id <- sub("_seqkit_stats.txt\$", "", basename(stats_file))

        # Extract values
        row <- data.frame(
            sample_id = sample_id,
            num_seqs = as.character(data[1, "num_seqs"]),
            sum_len = as.character(data[1, "sum_len"]),
            min_len = as.character(data[1, "min_len"]),
            avg_len = as.character(data[1, "avg_len"]),
            max_len = as.character(data[1, "max_len"]),
            Q20_percent = as.character(data[1, "Q20(%)"]),
            Q30_percent = as.character(data[1, "Q30(%)"]),
            avg_qual = as.character(data[1, "AvgQual"]),
            gc_percent = as.character(data[1, "GC(%)"]),
            stringsAsFactors = FALSE
        )

        # Add this row to results
        results <- rbind(results, row)
    }

    # Sort by sample_id
    results <- results[order(results\$sample_id), ]
    rownames(results) <- NULL

    # Write CSV
    write.csv(results, "trimmed_seqkit_summary.csv", row.names = FALSE)

    # Write text summary
    sink("trimmed_seqkit_summary.txt")
    cat("================================================================================\\n")
    cat("                   SEQKIT SUMMARY - TRIMMED LONG READS\\n")
    cat("================================================================================\\n\\n")
    cat("Sequence statistics after Pychopper and Fastplong processing\\n")
    cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n\\n")

    for (i in 1:nrow(results)) {
        cat("--------------------------------------------------------------------------------\\n")
        cat(sprintf("SAMPLE: %s\\n", results\$sample_id[i]))
        cat("--------------------------------------------------------------------------------\\n")
        cat(sprintf("  Number of sequences:   %s\\n", results\$num_seqs[i]))
        cat(sprintf("  Total bases:           %s\\n", results\$sum_len[i]))
        cat(sprintf("  Min length:            %s bp\\n", results\$min_len[i]))
        cat(sprintf("  Average length:        %s bp\\n", results\$avg_len[i]))
        cat(sprintf("  Max length:            %s bp\\n", results\$max_len[i]))
        cat(sprintf("  Q20%%:                  %s%%\\n", results\$Q20_percent[i]))
        cat(sprintf("  Q30%%:                  %s%%\\n", results\$Q30_percent[i]))
        cat(sprintf("  Average quality:       %s\\n", results\$avg_qual[i]))
        cat(sprintf("  GC content:            %s%%\\n", results\$gc_percent[i]))
        cat("\\n")
    }

    cat("================================================================================\\n")
    cat("SUMMARY:\\n")
    cat("--------------------------------------------------------------------------------\\n")
    cat(sprintf("Total samples: %d\\n", nrow(results)))
    cat("================================================================================\\n")
    sink()
    """
}
