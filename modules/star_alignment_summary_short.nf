// Process: STAR Alignment Summary - Extract mapping statistics
process STAR_ALIGNMENT_SUMMARY_SHORT {
    container 'https://depot.galaxyproject.org/singularity/r-base:4.3.1'
    publishDir "${params.output_dir}/short_read/03_alignment/summaries", mode: 'copy'

    input:
    path log_files

    output:
    path "star_alignment_summary.csv", emit: summary
    path "star_alignment_summary.txt", emit: summary_txt

    script:
    """
    #!/usr/bin/env Rscript

    log_files <- list.files(pattern = "Log.final.out\$", full.names = TRUE)

    results <- data.frame(
        sample_id = character(),
        input_reads = numeric(),
        uniquely_mapped = numeric(),
        uniquely_mapped_percent = numeric(),
        multi_mapped = numeric(),
        multi_mapped_percent = numeric(),
        multi_mapped_toomany = numeric(),
        multi_mapped_toomany_percent = numeric(),
        unmapped_too_short = numeric(),
        unmapped_too_short_percent = numeric(),
        unmapped_other = numeric(),
        unmapped_other_percent = numeric(),
        total_mapped = numeric(),
        total_mapped_percent = numeric(),
        stringsAsFactors = FALSE
    )

    for (log_file in log_files) {
        sample_id <- sub("Log.final.out\$", "", basename(log_file))

        lines <- readLines(log_file)

        input_reads <- NA
        uniquely_mapped <- NA
        uniquely_mapped_pct <- NA
        multi_mapped <- NA
        multi_mapped_pct <- NA
        multi_toomany <- NA
        multi_toomany_pct <- NA
        unmapped_short <- NA
        unmapped_short_pct <- NA
        unmapped_other <- NA
        unmapped_other_pct <- NA

        for (line in lines) {
            line <- trimws(line)
            parts <- strsplit(line, "|", fixed = TRUE)[[1]]

            if (length(parts) >= 2) {
                key <- trimws(parts[1])
                value <- trimws(parts[2])

                if (grepl("Number of input reads", key)) {
                    input_reads <- as.numeric(value)
                }
                else if (grepl("Uniquely mapped reads number", key)) {
                    uniquely_mapped <- as.numeric(value)
                }
                else if (grepl("Uniquely mapped reads %", key)) {
                    uniquely_mapped_pct <- as.numeric(sub("%", "", value))
                }
                else if (grepl("Number of reads mapped to multiple loci", key)) {
                    multi_mapped <- as.numeric(value)
                }
                else if (grepl("% of reads mapped to multiple loci", key)) {
                    multi_mapped_pct <- as.numeric(sub("%", "", value))
                }
                else if (grepl("Number of reads mapped to too many loci", key)) {
                    multi_toomany <- as.numeric(value)
                }
                else if (grepl("% of reads mapped to too many loci", key)) {
                    multi_toomany_pct <- as.numeric(sub("%", "", value))
                }
                else if (grepl("% of reads unmapped: too short", key)) {
                    unmapped_short_pct <- as.numeric(sub("%", "", value))
                }
                else if (grepl("% of reads unmapped: other", key)) {
                    unmapped_other_pct <- as.numeric(sub("%", "", value))
                }
            }
        }

        if (!is.na(input_reads) && !is.na(unmapped_short_pct)) {
            unmapped_short <- round(input_reads * (unmapped_short_pct / 100))
        }

        if (!is.na(input_reads) && !is.na(unmapped_other_pct)) {
            unmapped_other <- round(input_reads * (unmapped_other_pct / 100))
        }

        if (!is.na(uniquely_mapped) && !is.na(multi_mapped)) {
            total_mapped <- uniquely_mapped + multi_mapped
            if (!is.na(input_reads) && input_reads > 0) {
                total_mapped_percent <- (total_mapped / input_reads) * 100
            } else {
                total_mapped_percent <- NA
            }
        } else {
            total_mapped <- NA
            total_mapped_percent <- NA
        }

        row <- data.frame(
            sample_id = sample_id,
            input_reads = input_reads,
            uniquely_mapped = uniquely_mapped,
            uniquely_mapped_percent = uniquely_mapped_pct,
            multi_mapped = multi_mapped,
            multi_mapped_percent = multi_mapped_pct,
            multi_mapped_toomany = multi_toomany,
            multi_mapped_toomany_percent = multi_toomany_pct,
            unmapped_too_short = unmapped_short,
            unmapped_too_short_percent = unmapped_short_pct,
            unmapped_other = unmapped_other,
            unmapped_other_percent = unmapped_other_pct,
            total_mapped = total_mapped,
            total_mapped_percent = total_mapped_percent,
            stringsAsFactors = FALSE
        )

        results <- rbind(results, row)
    }

    results <- results[order(results\$sample_id), ]
    rownames(results) <- NULL

    write.csv(results, "star_alignment_summary.csv", row.names = FALSE)

    sink("star_alignment_summary.txt")
    cat("================================================================================\\n")
    cat("                        STAR ALIGNMENT SUMMARY\\n")
    cat("================================================================================\\n\\n")
    cat("Mapping statistics from STAR aligner\\n")
    cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n\\n")

    for (i in 1:nrow(results)) {
        cat("--------------------------------------------------------------------------------\\n")
        cat(sprintf("SAMPLE: %s\\n", results\$sample_id[i]))
        cat("--------------------------------------------------------------------------------\\n")

        cat("INPUT:\\n")
        cat(sprintf("  Total reads:              %s\\n",
                   format(results\$input_reads[i], big.mark = ",")))

        cat("\\nMAPPED READS:\\n")
        cat(sprintf("  Uniquely mapped:          %s (%.2f%%)\\n",
                   format(results\$uniquely_mapped[i], big.mark = ","),
                   results\$uniquely_mapped_percent[i]))
        cat(sprintf("  Multi-mapped:             %s (%.2f%%)\\n",
                   format(results\$multi_mapped[i], big.mark = ","),
                   results\$multi_mapped_percent[i]))
        cat(sprintf("  Mapped to too many loci:  %s (%.2f%%)\\n",
                   format(results\$multi_mapped_toomany[i], big.mark = ","),
                   results\$multi_mapped_toomany_percent[i]))
        cat(sprintf("  TOTAL MAPPED:             %s (%.2f%%)\\n",
                   format(results\$total_mapped[i], big.mark = ","),
                   results\$total_mapped_percent[i]))

        cat("\\nUNMAPPED READS:\\n")
        cat(sprintf("  Too short:                %s (%.2f%%)\\n",
                   format(results\$unmapped_too_short[i], big.mark = ","),
                   results\$unmapped_too_short_percent[i]))
        cat(sprintf("  Other reasons:            %s (%.2f%%)\\n",
                   format(results\$unmapped_other[i], big.mark = ","),
                   results\$unmapped_other_percent[i]))

        cat("\\n")
    }

    cat("================================================================================\\n")
    cat("OVERALL SUMMARY:\\n")
    cat("--------------------------------------------------------------------------------\\n")
    cat(sprintf("Total samples:              %d\\n", nrow(results)))
    cat(sprintf("Mean unique mapping rate:   %.2f%%\\n",
               mean(results\$uniquely_mapped_percent, na.rm = TRUE)))
    cat(sprintf("Mean total mapping rate:    %.2f%%\\n",
               mean(results\$total_mapped_percent, na.rm = TRUE)))
    cat(sprintf("Range of unique mapping:    %.2f%% - %.2f%%\\n",
               min(results\$uniquely_mapped_percent, na.rm = TRUE),
               max(results\$uniquely_mapped_percent, na.rm = TRUE)))
    cat("================================================================================\\n")
    sink()
    """
}
