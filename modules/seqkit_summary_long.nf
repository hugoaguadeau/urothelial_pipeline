// Process: SeqKit Summary - Combine all stats into one CSV (Long Reads - Raw)
process SEQKIT_SUMMARY_LONG {
    container 'https://depot.galaxyproject.org/singularity/r-base:4.3.1'
    publishDir "${params.output_dir}/long_read/01_quality_control/raw/summaries", mode: 'copy'

    input:
    path seqkit_files  // All the individual seqkit stats files

    output:
    path "seqkit_summary.csv", emit: summary

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
        filename <- as.character(data[1, "file"])
        
        # Extract sample ID from the filename in the stats
        # The file column contains the original fastq name
        sample_id <- sub("_seqkit_stats.txt\$", "", basename(stats_file))

        # Extract values - use column names from seqkit output
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

    # Reset row names
    rownames(results) <- NULL

    # Write to CSV
    write.csv(results, "seqkit_summary.csv", row.names = FALSE)
    """
}
