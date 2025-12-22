// Process: FeatureCounts Summary - Extract assignment statistics
process FEATURECOUNTS_SUMMARY_SHORT {
    container 'https://depot.galaxyproject.org/singularity/r-base:4.3.1'
    publishDir "${params.output_dir}/short_read/04_quantification/summaries", mode: 'copy'

    input:
    path summary_file

    output:
    path "featurecounts_summary.csv", emit: summary
    path "featurecounts_summary.txt", emit: summary_txt

    script:
    """
    #!/usr/bin/env Rscript

    data <- read.table("${summary_file}", header = TRUE, sep = "\\t",
                      stringsAsFactors = FALSE, check.names = FALSE)

    status_col <- data[, 1]
    sample_cols <- colnames(data)[-1]

    sample_cols <- basename(sample_cols)
    sample_cols <- sub("[.]sorted[.]bam\$", "", sample_cols)
    sample_cols <- sub("Aligned[.]sortedByCoord[.]out[.]bam\$", "", sample_cols)

    results <- data.frame(
        sample_id = character(),
        assigned = numeric(),
        assigned_percent = numeric(),
        unassigned_unmapped = numeric(),
        unassigned_read_type = numeric(),
        unassigned_singleton = numeric(),
        unassigned_mappingquality = numeric(),
        unassigned_chimera = numeric(),
        unassigned_fragmentlength = numeric(),
        unassigned_duplicate = numeric(),
        unassigned_multimapping = numeric(),
        unassigned_secondary = numeric(),
        unassigned_nonSplit = numeric(),
        unassigned_nofeatures = numeric(),
        unassigned_overlapping_length = numeric(),
        unassigned_ambiguity = numeric(),
        total_reads = numeric(),
        stringsAsFactors = FALSE
    )

    for (i in 2:ncol(data)) {
        sample_id <- sample_cols[i-1]

        assigned <- data[data[,1] == "Assigned", i]
        unmapped <- data[data[,1] == "Unassigned_Unmapped", i]
        read_type <- data[data[,1] == "Unassigned_Read_Type", i]
        singleton <- data[data[,1] == "Unassigned_Singleton", i]
        mapping_quality <- data[data[,1] == "Unassigned_MappingQuality", i]
        chimera <- data[data[,1] == "Unassigned_Chimera", i]
        frag_length <- data[data[,1] == "Unassigned_FragmentLength", i]
        duplicate <- data[data[,1] == "Unassigned_Duplicate", i]
        multi_mapping <- data[data[,1] == "Unassigned_MultiMapping", i]
        secondary <- data[data[,1] == "Unassigned_Secondary", i]
        non_split <- data[data[,1] == "Unassigned_NonSplit", i]
        no_features <- data[data[,1] == "Unassigned_NoFeatures", i]
        overlap_length <- data[data[,1] == "Unassigned_Overlapping_Length", i]
        ambiguity <- data[data[,1] == "Unassigned_Ambiguity", i]

        safe_get <- function(x) if(length(x) == 0) 0 else x

        assigned <- safe_get(assigned)
        unmapped <- safe_get(unmapped)
        read_type <- safe_get(read_type)
        singleton <- safe_get(singleton)
        mapping_quality <- safe_get(mapping_quality)
        chimera <- safe_get(chimera)
        frag_length <- safe_get(frag_length)
        duplicate <- safe_get(duplicate)
        multi_mapping <- safe_get(multi_mapping)
        secondary <- safe_get(secondary)
        non_split <- safe_get(non_split)
        no_features <- safe_get(no_features)
        overlap_length <- safe_get(overlap_length)
        ambiguity <- safe_get(ambiguity)

        total_reads <- sum(assigned, unmapped, read_type, singleton,
                          mapping_quality, chimera, frag_length, duplicate,
                          multi_mapping, secondary, non_split, no_features,
                          overlap_length, ambiguity)

        if (total_reads > 0) {
            assigned_pct <- (assigned / total_reads) * 100
        } else {
            assigned_pct <- 0
        }

        row <- data.frame(
            sample_id = sample_id,
            assigned = assigned,
            assigned_percent = assigned_pct,
            unassigned_unmapped = unmapped,
            unassigned_read_type = read_type,
            unassigned_singleton = singleton,
            unassigned_mappingquality = mapping_quality,
            unassigned_chimera = chimera,
            unassigned_fragmentlength = frag_length,
            unassigned_duplicate = duplicate,
            unassigned_multimapping = multi_mapping,
            unassigned_secondary = secondary,
            unassigned_nonSplit = non_split,
            unassigned_nofeatures = no_features,
            unassigned_overlapping_length = overlap_length,
            unassigned_ambiguity = ambiguity,
            total_reads = total_reads,
            stringsAsFactors = FALSE
        )

        results <- rbind(results, row)
    }

    results <- results[order(results\$sample_id), ]
    rownames(results) <- NULL

    write.csv(results, "featurecounts_summary.csv", row.names = FALSE)

    sink("featurecounts_summary.txt")
    cat("================================================================================\\n")
    cat("                     FEATURECOUNTS ASSIGNMENT SUMMARY\\n")
    cat("================================================================================\\n\\n")
    cat("Gene assignment statistics from featureCounts\\n")
    cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n\\n")

    for (i in 1:nrow(results)) {
        cat("--------------------------------------------------------------------------------\\n")
        cat(sprintf("SAMPLE: %s\\n", results\$sample_id[i]))
        cat("--------------------------------------------------------------------------------\\n")

        cat(sprintf("Total reads:              %s\\n",
                   format(results\$total_reads[i], big.mark = ",")))
        cat(sprintf("\\nASSIGNED TO GENES:        %s (%.2f%%)\\n",
                   format(results\$assigned[i], big.mark = ","),
                   results\$assigned_percent[i]))

        cat("\\nUNASSIGNED READS:\\n")

        if (results\$unassigned_unmapped[i] > 0) {
            cat(sprintf("  Unmapped:               %s (%.2f%%)\\n",
                       format(results\$unassigned_unmapped[i], big.mark = ","),
                       (results\$unassigned_unmapped[i]/results\$total_reads[i])*100))
        }
        if (results\$unassigned_multimapping[i] > 0) {
            cat(sprintf("  Multi-mapping:          %s (%.2f%%)\\n",
                       format(results\$unassigned_multimapping[i], big.mark = ","),
                       (results\$unassigned_multimapping[i]/results\$total_reads[i])*100))
        }
        if (results\$unassigned_nofeatures[i] > 0) {
            cat(sprintf("  No features:            %s (%.2f%%)\\n",
                       format(results\$unassigned_nofeatures[i], big.mark = ","),
                       (results\$unassigned_nofeatures[i]/results\$total_reads[i])*100))
        }
        if (results\$unassigned_ambiguity[i] > 0) {
            cat(sprintf("  Ambiguous:              %s (%.2f%%)\\n",
                       format(results\$unassigned_ambiguity[i], big.mark = ","),
                       (results\$unassigned_ambiguity[i]/results\$total_reads[i])*100))
        }
        if (results\$unassigned_secondary[i] > 0) {
            cat(sprintf("  Secondary alignments:   %s (%.2f%%)\\n",
                       format(results\$unassigned_secondary[i], big.mark = ","),
                       (results\$unassigned_secondary[i]/results\$total_reads[i])*100))
        }
        if (results\$unassigned_mappingquality[i] > 0) {
            cat(sprintf("  Low mapping quality:    %s (%.2f%%)\\n",
                       format(results\$unassigned_mappingquality[i], big.mark = ","),
                       (results\$unassigned_mappingquality[i]/results\$total_reads[i])*100))
        }

        cat("\\n")
    }

    cat("================================================================================\\n")
    cat("OVERALL SUMMARY:\\n")
    cat("--------------------------------------------------------------------------------\\n")
    cat(sprintf("Total samples:              %d\\n", nrow(results)))
    cat(sprintf("Mean assignment rate:       %.2f%%\\n",
               mean(results\$assigned_percent, na.rm = TRUE)))
    cat(sprintf("Range of assignment:        %.2f%% - %.2f%%\\n",
               min(results\$assigned_percent, na.rm = TRUE),
               max(results\$assigned_percent, na.rm = TRUE)))

    cat("\\nMean unassigned by category:\\n")
    mean_no_features <- mean((results\$unassigned_nofeatures/results\$total_reads)*100, na.rm=TRUE)
    mean_ambiguous <- mean((results\$unassigned_ambiguity/results\$total_reads)*100, na.rm=TRUE)
    mean_multimapping <- mean((results\$unassigned_multimapping/results\$total_reads)*100, na.rm=TRUE)

    cat(sprintf("  No features:              %.2f%%\\n", mean_no_features))
    cat(sprintf("  Ambiguous:                %.2f%%\\n", mean_ambiguous))
    cat(sprintf("  Multi-mapping:            %.2f%%\\n", mean_multimapping))
    cat("================================================================================\\n")
    sink()
    """
}
