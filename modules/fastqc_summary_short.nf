// Process: FastQC Summary - Extract key metrics into CSV
process FASTQC_SUMMARY_SHORT {
    container 'https://depot.galaxyproject.org/singularity/r-base:4.3.1'
    publishDir "${params.output_dir}/short_read/01_quality_control/raw/summaries", mode: 'copy'

    input:
    path fastqc_zips  // All the .zip files from FastQC

    output:
    path "fastqc_summary.csv", emit: summary

    script:
    """
    #!/usr/bin/env Rscript

    # Get list of all FastQC zip files in current directory
    zip_files <- list.files(pattern = "*_fastqc.zip", full.names = TRUE)
    
    # Create empty data frame to store results
    results <- data.frame(
        sample_id = character(),
        read = character(),
        total_sequences = character(),
        sequence_length = character(),
        percent_gc = character(),
        basic_statistics = character(),
        per_base_quality = character(),
        per_sequence_quality = character(),
        per_base_content = character(),
        duplication_levels = character(),
        adapter_content = character(),
        stringsAsFactors = FALSE
    )
    
    # Process each zip file
    for (zip_file in zip_files) {
        
        # Extract the folder name (remove .zip extension)
        # Using [.] instead of \\. to avoid escape character issues
        folder_name <- sub("[.]zip\$", "", basename(zip_file))
        
        # The data file inside the zip
        data_file <- file.path(folder_name, "fastqc_data.txt")
        
        # Read the fastqc_data.txt from inside the zip file
        # unz() creates a connection to a file inside a zip archive
        con <- unz(zip_file, data_file)
        lines <- readLines(con)
        close(con)
        
        # Parse the filename to get sample ID and read number
        # Format: Y2391-U_read1_fastqc.zip
        filename <- sub("_fastqc[.]zip\$", "", basename(zip_file))
        
        # Determine read number and sample ID
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
        
        # Create a row for this sample
        row <- data.frame(
            sample_id = sample_id,
            read = read_num,
            total_sequences = NA,
            sequence_length = NA,
            percent_gc = NA,
            basic_statistics = NA,
            per_base_quality = NA,
            per_sequence_quality = NA,
            per_base_content = NA,
            duplication_levels = NA,
            adapter_content = NA,
            stringsAsFactors = FALSE
        )
        
        # Extract metrics from the lines of fastqc_data.txt
        for (i in seq_along(lines)) {
            line <- lines[i]
            
            # Split line by tabs
            parts <- strsplit(line, "\\t")[[1]]
            
            # Extract different metrics based on line content
            if (grepl("^Total Sequences", line)) {
                row\$total_sequences <- parts[2]
            }
            else if (grepl("^Sequence length", line)) {
                row\$sequence_length <- parts[2]
            }
            else if (grepl("^%GC", line)) {
                row\$percent_gc <- parts[2]
            }
            # Quality check results - these start with >>
            else if (grepl("^>>Basic Statistics", line)) {
                # Last element is the status (pass/warn/fail)
                row\$basic_statistics <- trimws(parts[length(parts)])
            }
            else if (grepl("^>>Per base sequence quality", line)) {
                row\$per_base_quality <- trimws(parts[length(parts)])
            }
            else if (grepl("^>>Per sequence quality scores", line)) {
                row\$per_sequence_quality <- trimws(parts[length(parts)])
            }
            else if (grepl("^>>Per base sequence content", line)) {
                row\$per_base_content <- trimws(parts[length(parts)])
            }
            else if (grepl("^>>Sequence Duplication Levels", line)) {
                row\$duplication_levels <- trimws(parts[length(parts)])
            }
            else if (grepl("^>>Adapter Content", line)) {
                row\$adapter_content <- trimws(parts[length(parts)])
            }
        }
        
        # Add this row to our results
        results <- rbind(results, row)
    }
    
    # Sort results by sample_id and read for consistent ordering
    results <- results[order(results\$sample_id, results\$read), ]
    
    # Reset row names to 1, 2, 3... instead of keeping old indices
    rownames(results) <- NULL
    
    # Write to CSV - this is the only output that matters
    write.csv(results, "fastqc_summary.csv", row.names = FALSE)
    """
}
