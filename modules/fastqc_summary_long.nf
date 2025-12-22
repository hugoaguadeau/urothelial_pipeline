// Process: FastQC Summary - Extract key metrics into CSV (Long Reads)
process FASTQC_SUMMARY_LONG {
    container 'https://depot.galaxyproject.org/singularity/r-base:4.3.1'
    publishDir "${params.output_dir}/long_read/01_quality_control/raw/summaries", mode: 'copy'

    input:
    path fastqc_zips

    output:
    path "fastqc_summary.csv", emit: summary
    path "fastqc_summary.txt", emit: summary_txt

    script:
    """
    #!/usr/bin/env Rscript

    zip_files <- list.files(pattern = "*_fastqc.zip", full.names = TRUE)

    results <- data.frame(
        sample_id = character(),
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

    for (zip_file in zip_files) {
        folder_name <- sub("[.]zip\$", "", basename(zip_file))
        data_file <- file.path(folder_name, "fastqc_data.txt")

        con <- unz(zip_file, data_file)
        lines <- readLines(con)
        close(con)

        sample_id <- sub("_fastqc[.]zip\$", "", basename(zip_file))
        sample_id <- sub("-[BU]\$", "", sample_id)

        row <- data.frame(
            sample_id = sample_id,
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

        for (i in seq_along(lines)) {
            line <- lines[i]
            parts <- strsplit(line, "\\t")[[1]]

            if (grepl("^Total Sequences", line)) {
                row\$total_sequences <- parts[2]
            }
            else if (grepl("^Sequence length", line)) {
                row\$sequence_length <- parts[2]
            }
            else if (grepl("^%GC", line)) {
                row\$percent_gc <- parts[2]
            }
            else if (grepl("^>>Basic Statistics", line)) {
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

        results <- rbind(results, row)
    }

    results <- results[order(results\$sample_id), ]
    rownames(results) <- NULL

    write.csv(results, "fastqc_summary.csv", row.names = FALSE)

    sink("fastqc_summary.txt")
    cat("================================================================================\\n")
    cat("                   FASTQC SUMMARY - LONG READS (RAW)\\n")
    cat("================================================================================\\n\\n")
    cat("Quality control metrics for raw long read data\\n")
    cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n\\n")

    for (i in 1:nrow(results)) {
        cat("--------------------------------------------------------------------------------\\n")
        cat(sprintf("SAMPLE: %s\\n", results\$sample_id[i]))
        cat("--------------------------------------------------------------------------------\\n")
        cat(sprintf("  Total sequences:       %s\\n", results\$total_sequences[i]))
        cat(sprintf("  Sequence length:       %s bp\\n", results\$sequence_length[i]))
        cat(sprintf("  GC content:            %s%%\\n", results\$percent_gc[i]))
        cat("\\n  Quality Check Results:\\n")
        cat(sprintf("    Basic Statistics:         %s\\n",
                   toupper(results\$basic_statistics[i])))
        cat(sprintf("    Per base quality:         %s\\n",
                   toupper(results\$per_base_quality[i])))
        cat(sprintf("    Per sequence quality:     %s\\n",
                   toupper(results\$per_sequence_quality[i])))
        cat(sprintf("    Per base content:         %s\\n",
                   toupper(results\$per_base_content[i])))
        cat(sprintf("    Duplication levels:       %s\\n",
                   toupper(results\$duplication_levels[i])))
        cat(sprintf("    Adapter content:          %s\\n",
                   toupper(results\$adapter_content[i])))
        cat("\\n")
    }

    cat("================================================================================\\n")
    cat("OVERALL SUMMARY:\\n")
    cat("--------------------------------------------------------------------------------\\n")
    cat(sprintf("Total samples: %d\\n\\n", nrow(results)))

    cat("Quality Check Summary:\\n")
    cat(sprintf("  Basic Statistics PASS:     %d/%d\\n",
               sum(results\$basic_statistics == "pass"), nrow(results)))
    cat(sprintf("  Per base quality PASS:     %d/%d\\n",
               sum(results\$per_base_quality == "pass"), nrow(results)))
    cat(sprintf("  Per sequence quality PASS: %d/%d\\n",
               sum(results\$per_sequence_quality == "pass"), nrow(results)))
    cat(sprintf("  Adapter content PASS:      %d/%d\\n",
               sum(results\$adapter_content == "pass"), nrow(results)))

    total_warnings <- sum(results\$basic_statistics == "warn") +
                     sum(results\$per_base_quality == "warn") +
                     sum(results\$per_sequence_quality == "warn") +
                     sum(results\$per_base_content == "warn") +
                     sum(results\$duplication_levels == "warn") +
                     sum(results\$adapter_content == "warn")

    total_failures <- sum(results\$basic_statistics == "fail") +
                     sum(results\$per_base_quality == "fail") +
                     sum(results\$per_sequence_quality == "fail") +
                     sum(results\$per_base_content == "fail") +
                     sum(results\$duplication_levels == "fail") +
                     sum(results\$adapter_content == "fail")

    cat(sprintf("\\nTotal warnings across all tests:  %d\\n", total_warnings))
    cat(sprintf("Total failures across all tests:  %d\\n", total_failures))
    cat("================================================================================\\n")
    sink()
    """
}
