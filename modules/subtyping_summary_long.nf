// Process: Generate comprehensive subtyping summary report
process SUBTYPING_SUMMARY_LONG {
    container 'https://depot.galaxyproject.org/singularity/r-base:4.3.1'
    publishDir "${params.output_dir}/long_read/07_subtyping", mode: 'copy'

    input:
    path subtypes
    path sample_info

    output:
    path "subtyping_summary_report.txt", emit: summary_txt
    path "subtype_distribution.csv", emit: distribution_csv
    path "subtypes_with_metadata.csv", emit: subtypes_annotated

    script:
    """
    #!/usr/bin/env Rscript

    # Read data
    subtypes <- read.delim("${subtypes}", stringsAsFactors = FALSE)
    metadata <- read.csv("${sample_info}", stringsAsFactors = FALSE)

    # Merge with metadata (base R version, no dplyr needed)
    combined <- merge(subtypes, metadata, by.x = "ID", by.y = "sample_id", all.x = TRUE)

    # Keep only relevant columns and save annotated subtypes
    subtypes_annotated <- combined[, c("ID", "consensusClass", "separationLevel", 
                                        "tissue", "condition", "aim2_group")]
    write.csv(subtypes_annotated, "subtypes_with_metadata.csv", row.names = FALSE)

    # Calculate summary statistics
    subtype_counts <- table(subtypes\$consensusClass)
    subtype_props <- prop.table(subtype_counts) * 100

    # Separation level statistics (base R version)
    sep_by_subtype <- aggregate(
        separationLevel ~ consensusClass,
        data = subtypes,
        FUN = function(x) {
            c(n = length(x),
              mean_sep = mean(x, na.rm = TRUE),
              median_sep = median(x, na.rm = TRUE),
              min_sep = min(x, na.rm = TRUE),
              max_sep = max(x, na.rm = TRUE))
        }
    )
    
    # Flatten the aggregated results
    sep_stats <- do.call(data.frame, sep_by_subtype)
    colnames(sep_stats) <- c("consensusClass", "n", "mean_sep", "median_sep", "min_sep", "max_sep")

    # Save distribution table
    write.csv(sep_stats, "subtype_distribution.csv", row.names = FALSE)

    # Generate text report
    sink("subtyping_summary_report.txt")

    cat("================================================================================\\n")
    cat("           BLADDER CANCER MOLECULAR SUBTYPING SUMMARY\\n")
    cat("                    (Consensus MIBC Classifier)\\n")
    cat("================================================================================\\n\\n")
    cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n\\n")

    cat("SAMPLES ANALYSED\\n")
    cat("--------------------------------------------------------------------------------\\n")
    cat(sprintf("Total samples: %d\\n\\n", nrow(subtypes)))

    cat("SUBTYPE DISTRIBUTION\\n")
    cat("--------------------------------------------------------------------------------\\n")
    for (i in 1:length(subtype_counts)) {
        subtype_name <- names(subtype_counts)[i]
        count <- subtype_counts[i]
        pct <- subtype_props[i]
        cat(sprintf("%-25s %3d samples (%.1f%%)\\n", subtype_name, count, pct))
    }

    cat("\\nSEPARATION LEVEL STATISTICS\\n")
    cat("--------------------------------------------------------------------------------\\n")
    cat("The separation level indicates classification confidence.\\n")
    cat("Higher values = more confident classification.\\n\\n")

    cat(sprintf("Overall mean separation:   %.3f\\n", mean(subtypes\$separationLevel, na.rm = TRUE)))
    cat(sprintf("Overall median separation: %.3f\\n\\n", median(subtypes\$separationLevel, na.rm = TRUE)))

    cat("Per-subtype separation levels:\\n\\n")
    print(sep_stats)

    cat("\\nSAMPLES WITH LOW CONFIDENCE (<0.3 separation)\\n")
    cat("--------------------------------------------------------------------------------\\n")
    low_conf <- subtypes[subtypes\$separationLevel < 0.3, ]
    if (nrow(low_conf) > 0) {
        cat(sprintf("%d samples with low confidence classification:\\n", nrow(low_conf)))
        for (i in 1:nrow(low_conf)) {
            cat(sprintf("  - %s: %s (sep = %.3f)\\n",
                       low_conf\$ID[i],
                       low_conf\$consensusClass[i],
                       low_conf\$separationLevel[i]))
        }
    } else {
        cat("All samples have adequate separation (>0.3)\\n")
    }

    cat("\\n================================================================================\\n")
    cat("DETAILED SAMPLE BREAKDOWN\\n")
    cat("================================================================================\\n\\n")
    
    # Print each sample with metadata
    for (i in 1:nrow(subtypes_annotated)) {
        sample <- subtypes_annotated[i, ]
        cat(sprintf("%-10s | %-8s | %-20s | %-15s | sep = %.3f\\n",
                   sample\$ID,
                   sample\$consensusClass,
                   sample\$condition,
                   sample\$tissue,
                   sample\$separationLevel))
    }

    cat("\\n================================================================================\\n")
    cat("SUBTYPE BY TISSUE\\n")
    cat("================================================================================\\n\\n")
    
    # Cross-tabulation of subtype by tissue
    tissue_subtype_table <- table(subtypes_annotated\$tissue, subtypes_annotated\$consensusClass)
    cat("Subtype distribution by tissue:\\n\\n")
    print(tissue_subtype_table)
    
    cat("\\n================================================================================\\n")
    cat("SUBTYPE BY CONDITION\\n")
    cat("================================================================================\\n\\n")
    
    # Cross-tabulation of subtype by condition
    condition_subtype_table <- table(subtypes_annotated\$condition, subtypes_annotated\$consensusClass)
    cat("Subtype distribution by condition:\\n\\n")
    print(condition_subtype_table)

    cat("\\n================================================================================\\n")
    cat("SUBTYPING DEFINITIONS\\n")
    cat("================================================================================\\n")
    cat("LumP     = Luminal Papillary \\n")
    cat("LumU     = Luminal Unstable \\n")
    cat("LumNS    = Luminal \\n")
    cat("Ba/Sq    = Basal/Squamous \\n")
    cat("Stroma   = Stroma-rich \\n")
    cat("NE-like  = Neuroendocrine-like \\n")
    cat("================================================================================\\n")

    sink()

    cat("Subtyping summary report generated successfully\\n")
    """
}
