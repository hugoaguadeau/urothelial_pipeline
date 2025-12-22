// Process: Pychopper Summary - Extract classification statistics
process PYCHOPPER_SUMMARY_LONG {
    container 'https://depot.galaxyproject.org/singularity/r-base:4.3.1'
    publishDir "${params.output_dir}/long_read/02_preprocessing/summaries", mode: 'copy'

    input:
    path pychopper_stats

    output:
    path "pychopper_summary.csv", emit: summary
    path "pychopper_summary.txt", emit: summary_txt

    script:
    """
    #!/usr/bin/env Rscript

    stats_files <- list.files(pattern = "*_stats.tsv\$", full.names = TRUE)

    results <- data.frame(
        sample_id = character(),
        input_reads = numeric(),
        pass_reads = numeric(),
        pass_reads_percent = numeric(),
        len_fail = numeric(),
        len_fail_percent = numeric(),
        qc_fail = numeric(),
        qc_fail_percent = numeric(),
        processed_reads = numeric(),
        primers_found = numeric(),
        primers_found_percent = numeric(),
        rescue = numeric(),
        rescue_percent = numeric(),
        unusable = numeric(),
        unusable_percent = numeric(),
        strand_plus = numeric(),
        strand_minus = numeric(),
        rescue_strand_plus = numeric(),
        rescue_strand_minus = numeric(),
        stringsAsFactors = FALSE
    )

    for (stats_file in stats_files) {
        sample_id <- sub("_stats.tsv\$", "", basename(stats_file))

        tryCatch({
            data <- read.table(stats_file, header = TRUE, sep = "\\t",
                              stringsAsFactors = FALSE)

            get_value <- function(category, name) {
                val <- data[data\$Category == category & data\$Name == name, "Value"]
                if (length(val) == 0) return(0)
                return(as.numeric(val[1]))
            }

            pass_reads <- get_value("ReadStats", "PassReads")
            len_fail <- get_value("ReadStats", "LenFail")
            qc_fail <- get_value("ReadStats", "QcFail")
            
            # CORRECTED: Input reads is what went INTO Pychopper (from raw FASTQ)
            # This should equal PassReads + QcFail (reads that were actually processed)
            # LenFail represents reads filtered out before processing
            input_reads <- pass_reads + qc_fail
            processed_reads <- pass_reads + qc_fail

            primers_found <- get_value("Classification", "Primers_found")
            rescue <- get_value("Classification", "Rescue")
            unusable <- get_value("Classification", "Unusable")

            strand_plus <- get_value("Strand", "+")
            strand_minus <- get_value("Strand", "-")
            rescue_strand_plus <- get_value("RescueStrand", "+")
            rescue_strand_minus <- get_value("RescueStrand", "-")

            if (input_reads > 0) {
                # Percentages relative to input reads
                pass_reads_pct <- (pass_reads / input_reads) * 100
                qc_fail_pct <- (qc_fail / input_reads) * 100
                # LenFail percentage relative to total that went through length filter
                len_fail_pct <- (len_fail / (input_reads + len_fail)) * 100
            } else {
                pass_reads_pct <- qc_fail_pct <- len_fail_pct <- 0
            }
            
            if (processed_reads > 0) {
                # Classification percentages relative to processed reads
                primers_found_pct <- (primers_found / processed_reads) * 100
                rescue_pct <- (rescue / processed_reads) * 100
                unusable_pct <- (unusable / processed_reads) * 100
            } else {
                primers_found_pct <- rescue_pct <- unusable_pct <- 0
            }

            row <- data.frame(
                sample_id = sample_id,
                input_reads = input_reads,
                pass_reads = pass_reads,
                pass_reads_percent = pass_reads_pct,
                len_fail = len_fail,
                len_fail_percent = len_fail_pct,
                qc_fail = qc_fail,
                qc_fail_percent = qc_fail_pct,
                processed_reads = processed_reads,
                primers_found = primers_found,
                primers_found_percent = primers_found_pct,
                rescue = rescue,
                rescue_percent = rescue_pct,
                unusable = unusable,
                unusable_percent = unusable_pct,
                strand_plus = strand_plus,
                strand_minus = strand_minus,
                rescue_strand_plus = rescue_strand_plus,
                rescue_strand_minus = rescue_strand_minus,
                stringsAsFactors = FALSE
            )

            results <- rbind(results, row)

        }, error = function(e) {
            cat("Warning: Could not parse", stats_file, ":", conditionMessage(e), "\\n")
        })
    }

    results <- results[order(results\$sample_id), ]
    rownames(results) <- NULL

    write.csv(results, "pychopper_summary.csv", row.names = FALSE)

    sink("pychopper_summary.txt")
    cat("================================================================================\\n")
    cat("                        PYCHOPPER CLASSIFICATION SUMMARY\\n")
    cat("================================================================================\\n\\n")
    cat("Full-length cDNA read identification and classification\\n")
    cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n\\n")

    for (i in 1:nrow(results)) {
        cat("--------------------------------------------------------------------------------\\n")
        cat(sprintf("SAMPLE: %s\\n", results\$sample_id[i]))
        cat("--------------------------------------------------------------------------------\\n")

        cat("  READ STATISTICS:\\n")
        cat(sprintf("    Input reads:          %s\\n",
                   format(results\$input_reads[i], big.mark = ",")))
        cat(sprintf("    Pass reads:           %s (%.2f%%)\\n",
                   format(results\$pass_reads[i], big.mark = ","),
                   results\$pass_reads_percent[i]))
        cat(sprintf("    Length fail:          %s (%.2f%%)\\n",
                   format(results\$len_fail[i], big.mark = ","),
                   results\$len_fail_percent[i]))
        cat(sprintf("    QC fail:              %s (%.2f%%)\\n",
                   format(results\$qc_fail[i], big.mark = ","),
                   results\$qc_fail_percent[i]))
        cat(sprintf("    Processed reads:      %s\\n",
                   format(results\$processed_reads[i], big.mark = ",")))

        cat("\\n  CLASSIFICATION:\\n")
        cat(sprintf("    Primers found:        %s (%.2f%%)\\n",
                   format(results\$primers_found[i], big.mark = ","),
                   results\$primers_found_percent[i]))
        cat(sprintf("    Rescued:              %s (%.2f%%)\\n",
                   format(results\$rescue[i], big.mark = ","),
                   results\$rescue_percent[i]))
        cat(sprintf("    Unusable:             %s (%.2f%%)\\n",
                   format(results\$unusable[i], big.mark = ","),
                   results\$unusable_percent[i]))

        cat("\\n  STRAND ORIENTATION:\\n")
        cat(sprintf("    Plus strand:          %s\\n",
                   format(results\$strand_plus[i], big.mark = ",")))
        cat(sprintf("    Minus strand:         %s\\n",
                   format(results\$strand_minus[i], big.mark = ",")))
        cat(sprintf("    Rescue plus:          %s\\n",
                   format(results\$rescue_strand_plus[i], big.mark = ",")))
        cat(sprintf("    Rescue minus:         %s\\n",
                   format(results\$rescue_strand_minus[i], big.mark = ",")))

        full_length <- results\$primers_found[i] + results\$rescue[i]
        full_length_pct <- (full_length / results\$input_reads[i]) * 100
        cat("\\n  FULL-LENGTH IDENTIFICATION:\\n")
        cat(sprintf("    Total full-length:    %s (%.2f%%)\\n",
                   format(full_length, big.mark = ","),
                   full_length_pct))
        cat("\\n")
    }

    cat("================================================================================\\n")
    cat("OVERALL SUMMARY:\\n")
    cat("--------------------------------------------------------------------------------\\n")
    cat(sprintf("Total samples processed: %d\\n", nrow(results)))
    cat(sprintf("Total input reads:       %s\\n",
               format(sum(results\$input_reads), big.mark = ",")))
    cat(sprintf("Total length fail:       %s\\n",
               format(sum(results\$len_fail), big.mark = ",")))
    cat(sprintf("Mean pass rate:          %.2f%%\\n",
               mean(results\$pass_reads_percent)))
    cat(sprintf("Mean primers found:      %.2f%%\\n",
               mean(results\$primers_found_percent)))
    cat(sprintf("Mean rescue rate:        %.2f%%\\n",
               mean(results\$rescue_percent)))

    total_full_length <- sum(results\$primers_found + results\$rescue)
    total_input_reads <- sum(results\$input_reads)
    overall_fl_pct <- (total_full_length / total_input_reads) * 100

    cat(sprintf("Overall full-length ID:  %.2f%%\\n", overall_fl_pct))
    cat("================================================================================\\n")
    sink()
    """
}
