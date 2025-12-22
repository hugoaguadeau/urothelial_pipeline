// Process: SeqKit Summary - Combine trimmed read stats into one CSV
process SEQKIT_SUMMARY_TRIMMED_SHORT {
    container 'https://depot.galaxyproject.org/singularity/r-base:4.3.1'
    publishDir "${params.output_dir}/short_read/01_quality_control/trimmed/summaries", mode: 'copy'

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
        read = character(),
        num_seqs = character(),
        sum_len = character(),
        avg_len = character(),
        Q20_percent = character(),
        Q30_percent = character(),
        avg_qual = character(),
        gc_percent = character(),
        stringsAsFactors = FALSE
    )
    
    # Process each stats file
    for (stats_file in stats_files) {
        
        # Read the file - let R auto-detect whitespace
        data <- read.table(stats_file, header = TRUE, 
                          stringsAsFactors = FALSE, 
                          comment.char = "",
                          check.names = FALSE)
        
        # Each file has 2 rows (R1 and R2)
        for (i in 1:nrow(data)) {
            
            # Parse filename to get sample_id and read
            # Format: Y2391_trimmed_R1.fastq.gz
            filename <- as.character(data[i, "file"])
            
            # Check if filename is valid
            if (is.na(filename) || length(filename) == 0 || filename == "") {
                next
            }
            
            # Remove .fastq.gz extension
            filename <- sub("[.]fastq[.]gz\$", "", filename)
            
            # Determine read number
            if (grepl("R1", filename)) {
                read_num <- "R1"
                sample_id <- sub("_trimmed_R1\$", "", filename)
            } else if (grepl("R2", filename)) {
                read_num <- "R2"
                sample_id <- sub("_trimmed_R2\$", "", filename)
            } else {
                read_num <- "unknown"
                sample_id <- filename
            }
            
            # Create a row with selected metrics
            row <- data.frame(
                sample_id = sample_id,
                read = read_num,
                num_seqs = as.character(data[i, "num_seqs"]),
                sum_len = as.character(data[i, "sum_len"]),
                avg_len = as.character(data[i, "avg_len"]),
                Q20_percent = as.character(data[i, "Q20(%)"]),
                Q30_percent = as.character(data[i, "Q30(%)"]),
                avg_qual = as.character(data[i, "AvgQual"]),
                gc_percent = as.character(data[i, "GC(%)"]),
                stringsAsFactors = FALSE
            )
            
            # Add this row to results
            results <- rbind(results, row)
        }
    }
    
    # Sort by sample_id and read
    results <- results[order(results\$sample_id, results\$read), ]
    
    # Reset row names
    rownames(results) <- NULL
    
    # Write to CSV
    write.csv(results, "trimmed_seqkit_summary.csv", row.names = FALSE)
    
    # Write text summary
    sink("trimmed_seqkit_summary.txt")
    cat("================================================================================\\n")
    cat("                   SEQKIT SUMMARY - TRIMMED READS\\n")
    cat("================================================================================\\n\\n")
    cat("Sequence statistics after adapter trimming and filtering\\n\\n")
    
    cat(sprintf("Total samples: %d\\n", length(unique(results\$sample_id))))
    cat(sprintf("Total files analyzed: %d\\n\\n", nrow(results)))
    
    # Calculate summary statistics
    q30_values <- as.numeric(results\$Q30_percent)
    q20_values <- as.numeric(results\$Q20_percent)
    
    cat("Quality Summary:\\n")
    cat(sprintf("  Mean Q30%%: %.2f%%\\n", mean(q30_values, na.rm = TRUE)))
    cat(sprintf("  Mean Q20%%: %.2f%%\\n", mean(q20_values, na.rm = TRUE)))
    cat(sprintf("  Q30%% range: %.2f%% - %.2f%%\\n", 
               min(q30_values, na.rm = TRUE), max(q30_values, na.rm = TRUE)))
    
    cat("\\n================================================================================\\n")
    sink()
    """
}
