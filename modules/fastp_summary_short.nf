// Process: Fastp Summary - Extract key metrics from JSON reports
process FASTP_SUMMARY_SHORT {
    container 'https://depot.galaxyproject.org/singularity/r-base:4.3.1'
    publishDir "${params.output_dir}/short_read/02_preprocessing/summaries", mode: 'copy'

    input:
    path fastp_jsons

    output:
    path "fastp_summary.csv", emit: summary
    path "fastp_summary.txt", emit: summary_txt

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

    json_files <- list.files(pattern = "*_fastp[.]json\$", full.names = TRUE)

    results <- data.frame(
        sample_id = character(),
        read_type = character(),
        before_total_reads = numeric(),
        before_total_bases = numeric(),
        before_q20_rate = numeric(),
        before_q30_rate = numeric(),
        before_gc_content = numeric(),
        after_total_reads = numeric(),
        after_total_bases = numeric(),
        after_q20_rate = numeric(),
        after_q30_rate = numeric(),
        after_gc_content = numeric(),
        reads_passed = numeric(),
        reads_filtered = numeric(),
        reads_with_adapter = numeric(),
        bases_trimmed = numeric(),
        stringsAsFactors = FALSE
    )

    for (json_file in json_files) {
        sample_id <- sub("_fastp[.]json\$", "", basename(json_file))

        tryCatch({
            data <- fromJSON(json_file, simplifyVector = TRUE)

            before_total_reads <- data\$summary\$before_filtering\$total_reads
            before_total_bases <- data\$summary\$before_filtering\$total_bases
            before_q20 <- data\$summary\$before_filtering\$q20_rate * 100
            before_q30 <- data\$summary\$before_filtering\$q30_rate * 100
            before_gc <- data\$summary\$before_filtering\$gc_content * 100

            after_total_reads <- data\$summary\$after_filtering\$total_reads
            after_total_bases <- data\$summary\$after_filtering\$total_bases
            after_q20 <- data\$summary\$after_filtering\$q20_rate * 100
            after_q30 <- data\$summary\$after_filtering\$q30_rate * 100
            after_gc <- data\$summary\$after_filtering\$gc_content * 100

            reads_passed <- data\$filtering_result\$passed_filter_reads
            reads_filtered <- data\$filtering_result\$low_quality_reads +
                             data\$filtering_result\$too_short_reads +
                             data\$filtering_result\$too_many_N_reads

            bases_trimmed <- before_total_bases - after_total_bases

            if (!is.null(data\$adapter_cutting) &&
                !is.null(data\$adapter_cutting\$adapter_trimmed_reads)) {
                reads_with_adapter <- data\$adapter_cutting\$adapter_trimmed_reads
            } else {
                reads_with_adapter <- 0
            }

            row <- data.frame(
                sample_id = sample_id,
                read_type = "PE",
                before_total_reads = before_total_reads,
                before_total_bases = before_total_bases,
                before_q20_rate = before_q20,
                before_q30_rate = before_q30,
                before_gc_content = before_gc,
                after_total_reads = after_total_reads,
                after_total_bases = after_total_bases,
                after_q20_rate = after_q20,
                after_q30_rate = after_q30,
                after_gc_content = after_gc,
                reads_passed = reads_passed,
                reads_filtered = reads_filtered,
                reads_with_adapter = reads_with_adapter,
                bases_trimmed = bases_trimmed,
                stringsAsFactors = FALSE
            )

            results <- rbind(results, row)

        }, error = function(e) {
            cat("Warning: Could not parse", json_file, ":", conditionMessage(e), "\\n")
        })
    }

    results <- results[order(results\$sample_id), ]
    rownames(results) <- NULL

    write.csv(results, "fastp_summary.csv", row.names = FALSE)

    sink("fastp_summary.txt")
    cat("================================================================================\\n")
    cat("                        FASTP TRIMMING SUMMARY\\n")
    cat("================================================================================\\n\\n")
    cat("Summary of adapter trimming and quality filtering\\n")
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
        cat(sprintf("    Q20 rate:         %.2f%%\\n", results\$before_q20_rate[i]))
        cat(sprintf("    Q30 rate:         %.2f%%\\n", results\$before_q30_rate[i]))
        cat(sprintf("    GC content:       %.2f%%\\n", results\$before_gc_content[i]))

        cat("\\n  AFTER FILTERING:\\n")
        cat(sprintf("    Total reads:      %s\\n",
                   format(results\$after_total_reads[i], big.mark = ",")))
        cat(sprintf("    Total bases:      %s\\n",
                   format(results\$after_total_bases[i], big.mark = ",")))
        cat(sprintf("    Q20 rate:         %.2f%%\\n", results\$after_q20_rate[i]))
        cat(sprintf("    Q30 rate:         %.2f%%\\n", results\$after_q30_rate[i]))
        cat(sprintf("    GC content:       %.2f%%\\n", results\$after_gc_content[i]))

        cat("\\n  FILTERING RESULTS:\\n")
        cat(sprintf("    Reads passed:     %s (%.2f%%)\\n",
                   format(results\$reads_passed[i], big.mark = ","),
                   (results\$reads_passed[i] / results\$before_total_reads[i]) * 100))
        cat(sprintf("    Reads filtered:   %s (%.2f%%)\\n",
                   format(results\$reads_filtered[i], big.mark = ","),
                   (results\$reads_filtered[i] / results\$before_total_reads[i]) * 100))
        cat(sprintf("    Reads with adapter: %s\\n",
                   format(results\$reads_with_adapter[i], big.mark = ",")))
        cat(sprintf("    Bases trimmed:    %s\\n",
                   format(results\$bases_trimmed[i], big.mark = ",")))

        q30_improvement <- results\$after_q30_rate[i] - results\$before_q30_rate[i]
        cat("\\n  QUALITY IMPROVEMENT:\\n")
        cat(sprintf("    Q30 change:       %+.2f%%\\n", q30_improvement))
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
    cat(sprintf("Total bases trimmed:     %s\\n",
               format(sum(results\$bases_trimmed), big.mark = ",")))
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
