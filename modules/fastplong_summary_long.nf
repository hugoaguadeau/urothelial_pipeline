// Process: Fastplong Summary - Extract key metrics from JSON reports
process FASTPLONG_SUMMARY_LONG {
    container 'https://depot.galaxyproject.org/singularity/r-base:4.3.1'
    publishDir "${params.output_dir}/long_read/02_preprocessing/summaries", mode: 'copy'

    input:
    path fastplong_jsons

    output:
    path "fastplong_summary.csv", emit: summary
    path "fastplong_summary.txt", emit: summary_txt

    script:
    """
    #!/usr/bin/env Rscript

    local_lib <- "./r_packages"
    dir.create(local_lib, showWarnings = FALSE, recursive = TRUE)
    .libPaths(c(local_lib, .libPaths()))

    if (!requireNamespace("jsonlite", quietly = TRUE)) {
        install.packages("jsonlite", lib = local_lib,
                        repos = "https://cloud.r-project.org/",
                        quiet = TRUE, dependencies = TRUE)
    }

    library(jsonlite)

    json_files <- list.files(pattern = "*_fastplong[.]json\$", full.names = TRUE)

    results <- data.frame(
        sample_id = character(),
        before_total_reads = numeric(),
        before_total_bases = numeric(),
        before_q20_rate = numeric(),
        before_q30_rate = numeric(),
        before_mean_length = numeric(),
        before_gc_content = numeric(),
        after_total_reads = numeric(),
        after_total_bases = numeric(),
        after_q20_rate = numeric(),
        after_q30_rate = numeric(),
        after_mean_length = numeric(),
        after_gc_content = numeric(),
        reads_passed = numeric(),
        reads_low_quality = numeric(),
        reads_too_short = numeric(),
        reads_filtered = numeric(),
        reads_filtered_percent = numeric(),
        adapter_trimmed_reads = numeric(),
        adapter_trimmed_bases = numeric(),
        polyx_trimmed_reads = numeric(),
        polyx_trimmed_bases = numeric(),
        bases_removed = numeric(),
        stringsAsFactors = FALSE
    )

    for (json_file in json_files) {
        sample_id <- sub("_fastplong[.]json\$", "", basename(json_file))

        tryCatch({
            data <- fromJSON(json_file, simplifyVector = TRUE)

            before_total_reads <- data\$summary\$before_filtering\$total_reads
            before_total_bases <- data\$summary\$before_filtering\$total_bases
            before_q20 <- data\$summary\$before_filtering\$q20_rate * 100
            before_q30 <- data\$summary\$before_filtering\$q30_rate * 100
            before_mean_length <- data\$summary\$before_filtering\$read_mean_length
            before_gc <- data\$summary\$before_filtering\$gc_content * 100

            after_total_reads <- data\$summary\$after_filtering\$total_reads
            after_total_bases <- data\$summary\$after_filtering\$total_bases
            after_q20 <- data\$summary\$after_filtering\$q20_rate * 100
            after_q30 <- data\$summary\$after_filtering\$q30_rate * 100
            after_mean_length <- data\$summary\$after_filtering\$read_mean_length
            after_gc <- data\$summary\$after_filtering\$gc_content * 100

            reads_passed <- data\$filtering_result\$passed_filter_reads
            reads_low_qual <- data\$filtering_result\$low_quality_reads
            reads_too_short <- data\$filtering_result\$too_short_reads

            reads_filtered <- reads_low_qual + reads_too_short
            if (before_total_reads > 0) {
                reads_filtered_pct <- (reads_filtered / before_total_reads) * 100
            } else {
                reads_filtered_pct <- 0
            }

            adapter_trimmed_reads <- data\$adapter_cutting\$adapter_trimmed_reads
            adapter_trimmed_bases <- data\$adapter_cutting\$adapter_trimmed_bases

            polyx_trimmed_reads <- data\$polyx_trimming\$total_polyx_trimmed_reads
            polyx_trimmed_bases <- data\$polyx_trimming\$total_polyx_trimmed_bases

            bases_removed <- before_total_bases - after_total_bases

            row <- data.frame(
                sample_id = sample_id,
                before_total_reads = before_total_reads,
                before_total_bases = before_total_bases,
                before_q20_rate = before_q20,
                before_q30_rate = before_q30,
                before_mean_length = before_mean_length,
                before_gc_content = before_gc,
                after_total_reads = after_total_reads,
                after_total_bases = after_total_bases,
                after_q20_rate = after_q20,
                after_q30_rate = after_q30,
                after_mean_length = after_mean_length,
                after_gc_content = after_gc,
                reads_passed = reads_passed,
                reads_low_quality = reads_low_qual,
                reads_too_short = reads_too_short,
                reads_filtered = reads_filtered,
                reads_filtered_percent = reads_filtered_pct,
                adapter_trimmed_reads = adapter_trimmed_reads,
                adapter_trimmed_bases = adapter_trimmed_bases,
                polyx_trimmed_reads = polyx_trimmed_reads,
                polyx_trimmed_bases = polyx_trimmed_bases,
                bases_removed = bases_removed,
                stringsAsFactors = FALSE
            )

            results <- rbind(results, row)

        }, error = function(e) {
            cat("Warning: Could not parse", json_file, ":", conditionMessage(e), "\\n")
        })
    }

    results <- results[order(results\$sample_id), ]
    rownames(results) <- NULL

    write.csv(results, "fastplong_summary.csv", row.names = FALSE)

    sink("fastplong_summary.txt")
    cat("================================================================================\\n")
    cat("                        FASTPLONG TRIMMING SUMMARY\\n")
    cat("================================================================================\\n\\n")
    cat("Summary of adapter trimming and quality filtering for long reads\\n")
    cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n\\n")

    for (i in 1:nrow(results)) {
        cat("--------------------------------------------------------------------------------\\n")
        cat(sprintf("SAMPLE: %s\\n", results\$sample_id[i]))
        cat("--------------------------------------------------------------------------------\\n")

        cat("  BEFORE FILTERING:\\n")
        cat(sprintf("    Total reads:      %s\\n",
                   format(results\$before_total_reads[i], big.mark = ",")))
        cat(sprintf("    Total bases:      %s\\n",
                   format(results\$before_total_bases[i], big.mark = ",")))
        cat(sprintf("    Mean length:      %.0f bp\\n", results\$before_mean_length[i]))
        cat(sprintf("    Q20 rate:         %.2f%%\\n", results\$before_q20_rate[i]))
        cat(sprintf("    Q30 rate:         %.2f%%\\n", results\$before_q30_rate[i]))
        cat(sprintf("    GC content:       %.2f%%\\n", results\$before_gc_content[i]))

        cat("\\n  AFTER FILTERING:\\n")
        cat(sprintf("    Total reads:      %s\\n",
                   format(results\$after_total_reads[i], big.mark = ",")))
        cat(sprintf("    Total bases:      %s\\n",
                   format(results\$after_total_bases[i], big.mark = ",")))
        cat(sprintf("    Mean length:      %.0f bp\\n", results\$after_mean_length[i]))
        cat(sprintf("    Q20 rate:         %.2f%%\\n", results\$after_q20_rate[i]))
        cat(sprintf("    Q30 rate:         %.2f%%\\n", results\$after_q30_rate[i]))
        cat(sprintf("    GC content:       %.2f%%\\n", results\$after_gc_content[i]))

        cat("\\n  FILTERING RESULTS:\\n")
        cat(sprintf("    Reads passed:     %s (%.2f%%)\\n",
                   format(results\$reads_passed[i], big.mark = ","),
                   (results\$reads_passed[i] / results\$before_total_reads[i]) * 100))
        cat(sprintf("    Low quality:      %s\\n",
                   format(results\$reads_low_quality[i], big.mark = ",")))
        cat(sprintf("    Too short:        %s\\n",
                   format(results\$reads_too_short[i], big.mark = ",")))
        cat(sprintf("    Total filtered:   %s (%.2f%%)\\n",
                   format(results\$reads_filtered[i], big.mark = ","),
                   results\$reads_filtered_percent[i]))

        cat("\\n  TRIMMING STATISTICS:\\n")
        cat(sprintf("    Adapter trimmed reads:  %s\\n",
                   format(results\$adapter_trimmed_reads[i], big.mark = ",")))
        cat(sprintf("    Adapter trimmed bases:  %s\\n",
                   format(results\$adapter_trimmed_bases[i], big.mark = ",")))
        cat(sprintf("    PolyX trimmed reads:    %s\\n",
                   format(results\$polyx_trimmed_reads[i], big.mark = ",")))
        cat(sprintf("    PolyX trimmed bases:    %s\\n",
                   format(results\$polyx_trimmed_bases[i], big.mark = ",")))
        cat(sprintf("    Total bases removed:    %s\\n",
                   format(results\$bases_removed[i], big.mark = ",")))

        q30_improvement <- results\$after_q30_rate[i] - results\$before_q30_rate[i]
        length_change <- results\$after_mean_length[i] - results\$before_mean_length[i]

        cat("\\n  QUALITY IMPROVEMENT:\\n")
        cat(sprintf("    Q30 change:       %+.2f%%\\n", q30_improvement))
        cat(sprintf("    Length change:    %+.0f bp\\n", length_change))
        cat("\\n")
    }

    cat("================================================================================\\n")
    cat("OVERALL SUMMARY:\\n")
    cat("--------------------------------------------------------------------------------\\n")
    cat(sprintf("Total samples processed: %d\\n", nrow(results)))
    cat(sprintf("Total reads before:      %s\\n",
               format(sum(results\$before_total_reads), big.mark = ",")))
    cat(sprintf("Total reads after:       %s\\n",
               format(sum(results\$after_total_reads), big.mark = ",")))
    cat(sprintf("Total bases removed:     %s\\n",
               format(sum(results\$bases_removed), big.mark = ",")))
    cat(sprintf("Mean Q30 before:         %.2f%%\\n",
               mean(results\$before_q30_rate)))
    cat(sprintf("Mean Q30 after:          %.2f%%\\n",
               mean(results\$after_q30_rate)))
    cat(sprintf("Mean Q30 improvement:    %+.2f%%\\n",
               mean(results\$after_q30_rate - results\$before_q30_rate)))
    cat(sprintf("Mean reads retained:     %.2f%%\\n",
               mean((results\$reads_passed / results\$before_total_reads) * 100)))
    cat("================================================================================\\n")
    sink()
    """
}
