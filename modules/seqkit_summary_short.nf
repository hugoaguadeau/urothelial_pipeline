// Process: SeqKit Summary - Combine all stats into one CSV
process SEQKIT_SUMMARY_SHORT {
    container 'https://depot.galaxyproject.org/singularity/r-base:4.3.1'
    publishDir "${params.output_dir}/short_read/01_quality_control/raw/summaries", mode: 'copy'

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
        
        # Read the file - use sep="" to let R figure out whitespace
        # This handles both tabs and multiple spaces
        data <- read.table(stats_file, header = TRUE, 
                          stringsAsFactors = FALSE, 
                          comment.char = "",
                          check.names = FALSE)
        
        # Each file has 2 rows (R1 and R2)
        for (i in 1:nrow(data)) {
            
            # Parse filename to get sample_id and read
            # Format: Y2391-U_read1.fastq.gz
            filename <- as.character(data[i, "file"])
            
            # Check if filename is valid (not NA or empty)
            if (is.na(filename) || length(filename) == 0 || filename == "") {
                next  # Skip this row
            }
            
            # Remove .fastq.gz extension
            filename <- sub("[.]fastq[.]gz\$", "", filename)
            
            # Determine read number
            if (grepl("read1", filename)) {
                read_num <- "R1"
                sample_id <- sub("_read1\$", "", filename)
            } else if (grepl("read2", filename)) {
                read_num <- "R2"
                sample_id <- sub("_read2\$", "", filename)
            } else {
                read_num <- "unknown"
                sample_id <- filename
            }
            
            # Access columns - using column names directly
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
    write.csv(results, "seqkit_summary.csv", row.names = FALSE)
    """
}
