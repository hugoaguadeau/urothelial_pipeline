process ISOFORM_SWITCH_PART2 {
    container 'https://depot.galaxyproject.org/singularity/bioconductor-isoformswitchanalyzer:2.6.0--r44h3df3fcb_0'
    publishDir "${params.output_dir}/short_read/06_biological_analysis/isoform_switches/part2", mode: 'copy'

    input:
    path part1_rds
    path cpat_result
    path pfam_result

    output:
    path "switch_list_final.rds",      emit: final_rds
    path "tables/*.csv",               emit: tables
    path "plots/*.pdf",                emit: plots
    path "top_genes_plots",            emit: gene_plots, optional: true
    path "logs/*.txt",                 emit: summary 

    script:
    """
    #!/usr/bin/env Rscript

    library(IsoformSwitchAnalyzeR)
    library(dplyr)
    library(ggplot2)

    # Define absolute path for the plotting output directory
    base_dir <- getwd()
    out_dir_plots <- file.path(base_dir, "top_genes_plots")

    # Create output directories
    dir.create("tables", showWarnings = FALSE)
    dir.create("plots", showWarnings = FALSE)
    dir.create(out_dir_plots, showWarnings = FALSE)
    dir.create("logs", showWarnings = FALSE)

    cat("=== Initialising Part 2 ===\\n")

    # 1. Load Part 1 State
    switchList <- readRDS("${part1_rds}")

    # 2. MANUAL CPAT INTEGRATION
    cat("Manually integrating CPAT results...\\n")
    cpat_df <- read.table("${cpat_result}", header=TRUE, sep="\\t", stringsAsFactors=FALSE)

    col_names <- colnames(cpat_df)
    id_col_idx <- 1
    prob_col_idx <- grep("prob", col_names, ignore.case=TRUE)
    if (length(prob_col_idx) == 0) stop("Error: 'Coding_prob' column not found in CPAT.")
    prob_col_idx <- prob_col_idx[length(prob_col_idx)]

    cpat_clean <- data.frame(
        isoform_id = cpat_df[, id_col_idx],
        codingPotentialValue = cpat_df[, prob_col_idx],
        stringsAsFactors = FALSE
    )
    cpat_clean\$codingPotential <- cpat_clean\$codingPotentialValue >= 0.725

    if("codingPotentialValue" %in% colnames(switchList\$isoformFeatures)) {
        switchList\$isoformFeatures\$codingPotentialValue <- NULL
        switchList\$isoformFeatures\$codingPotential <- NULL
    }
    switchList\$isoformFeatures <- dplyr::left_join(switchList\$isoformFeatures, cpat_clean, by = "isoform_id")

    non_coding_isos <- cpat_clean\$isoform_id[ !cpat_clean\$codingPotential ]
    switchList\$isoformFeatures\$PTC[ switchList\$isoformFeatures\$isoform_id %in% non_coding_isos ] <- NA
    if (!is.null(switchList\$orfAnalysis)) {
         switchList\$orfAnalysis <- switchList\$orfAnalysis[ !(switchList\$orfAnalysis\$isoform_id %in% non_coding_isos), ]
    }

    # 3. Import Pfam Results
    cat("Importing Pfam results...\\n")
    tryCatch({
        switchList <- analyzePFAM(
            switchAnalyzeRlist = switchList,
            pathToPFAMresultFile = "${pfam_result}",
            showProgress = FALSE
        )
    }, error = function(e) { cat("Warning: analyzePFAM failed.\\n") })

    # 4. Analyse Alternative Splicing
    cat("Analysing Alternative Splicing...\\n")
    switchList <- analyzeAlternativeSplicing(switchList, quiet = TRUE)

    # 5. Predict Functional Consequences
    cat("Predicting Switch Consequences...\\n")
    cons_to_run <- c('intron_retention','coding_potential','ORF_seq_similarity','NMD_status')
    if (!is.null(switchList\$domainAnalysis)) {
        cons_to_run <- c(cons_to_run, 'domains_identified')
    }

    switchList <- analyzeSwitchConsequences(
        switchAnalyzeRlist = switchList,
        consequencesToAnalyze = cons_to_run,
        dIFcutoff = 0.1,
        showProgress = FALSE
    )

    # =========================================================================
    # REPORTING SECTION
    # =========================================================================

    cat("Generating Report Files...\\n")

    splicing_sum <- extractSplicingSummary(switchList, returnResult=TRUE, plot=FALSE)
    write.csv(splicing_sum, "tables/splicing_summary.csv", row.names=FALSE)

    pdf("plots/splicing_summary_plot.pdf", width=8, height=6)
    extractSplicingSummary(switchList, plot=TRUE, returnResult=FALSE)
    dev.off()

    splicing_enrich <- extractSplicingEnrichment(switchList, returnResult=TRUE, plot=FALSE)
    write.csv(splicing_enrich, "tables/splicing_enrichment.csv", row.names=FALSE)

    pdf("plots/splicing_enrichment_plot.pdf", width=8, height=6)
    extractSplicingEnrichment(switchList, plot=TRUE, returnResult=FALSE)
    dev.off()

    cons_sum <- extractConsequenceSummary(switchList, consequencesToAnalyze=cons_to_run, returnResult=TRUE, plot=FALSE)
    write.csv(cons_sum, "tables/consequence_summary.csv", row.names=FALSE)

    pdf("plots/consequence_summary_plot.pdf", width=10, height=8)
    extractConsequenceSummary(switchList, consequencesToAnalyze=cons_to_run, plot=TRUE, returnResult=FALSE)
    dev.off()

    cons_enrich <- extractConsequenceEnrichment(switchList, consequencesToAnalyze=cons_to_run, returnResult=TRUE, plot=FALSE)
    write.csv(cons_enrich, "tables/consequence_enrichment.csv", row.names=FALSE)

    pdf("plots/consequence_enrichment_plot.pdf", width=8, height=6)
    extractConsequenceEnrichment(switchList, consequencesToAnalyze=cons_to_run, plot=TRUE, returnResult=FALSE)
    dev.off()

    ### B. TOP CANDIDATES ###

    top_switches <- extractTopSwitches(
        switchList,
        filterForConsequences = TRUE,
        n = Inf,
        sortByQvals = TRUE
    )

    write.csv(top_switches, "tables/all_significant_switches_full.csv", row.names=FALSE)

    cols_to_keep <- c("gene_name", "gene_id", "condition_1", "condition_2",
                      "gene_switch_q_value", "isoform_switch_q_value",
                      "dIF", "switchConsequencesGene", "IF1", "IF2")

    existing_cols <- cols_to_keep[cols_to_keep %in% colnames(top_switches)]
    readable_table <- top_switches[, existing_cols]

    write.csv(readable_table, "tables/top_switches_simplified.csv", row.names=FALSE)

    ### C. VISUALISATIONS (SWITCH PLOTS) ###

    cat("Generating Top 20 Switch Plots...\\n")
    tryCatch({
        switchPlotTopSwitches(
            switchAnalyzeRlist = switchList,
            n = 20,
            filterForConsequences = TRUE,
            pathToOutput = out_dir_plots,
            fileType = "pdf",
            splitComparison = FALSE,
            splitFunctionalConsequences = FALSE
        )
    }, error = function(e) {
        cat("Warning: Could not generate top switch plots.\\n")
        print(e)
    })

    ### D. SAVE R OBJECTS ###
    saveRDS(switchList, "switch_list_final.rds")

    cat("Analysis Complete.\\n")
    
    # Create a log file because the output channel expects it
    sink("logs/final_status.txt")
    cat("Part 2 Analysis Completed Successfully.\\n")
    sink()
    """
}
