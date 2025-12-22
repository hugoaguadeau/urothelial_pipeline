process ISOFORM_SWITCH_PART1 {
    container 'https://depot.galaxyproject.org/singularity/bioconductor-isoformswitchanalyzer:2.6.0--r44h3df3fcb_0'
    publishDir "${params.output_dir}/short_read/06_biological_analysis/isoform_switches/part1", mode: 'copy'

    input:
    path transcript_counts
    path sample_info
    path gtf_file
    path genome_fasta

    output:
    path "switch_list_part1.rds",      emit: switch_list_rds
    path "isoform_switch_nt.fasta",    emit: nt_fasta
    path "isoform_switch_AA.fasta",    emit: aa_fasta
    path "part1_summary.txt",          emit: summary

    script:
    """
    #!/usr/bin/env Rscript

    # 1. SETUP LOCAL LIBRARY
    dir.create("r_libs", showWarnings = FALSE)
    .libPaths(c("r_libs", .libPaths()))

    # 2. INSTALL BSGENOME
    options(repos = c(CRAN = "https://cloud.r-project.org"))

    if (!require("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
    }

    if (!require("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
        cat("=== Installing BSgenome... ===\\n")
        BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", update=FALSE, ask=FALSE)
    }

    # 3. LOAD LIBRARIES
    library(IsoformSwitchAnalyzeR)
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(dplyr)

    cat("=== Initialising Part 1 ===\\n")

    # 4. PREPARE DATA
    counts_raw <- read.csv("${transcript_counts}", row.names=1, check.names=FALSE)
    samples_data <- read.csv("${sample_info}")

    samples_data <- samples_data[match(colnames(counts_raw), samples_data\$sample_id), ]

    design_matrix <- data.frame(
        sampleID = samples_data\$sample_id,
        condition = samples_data\$tissue
    )

    # 5. CREATE SWITCHLIST
    cat("Importing Rdata...\\n")
    switchList <- importRdata(
        isoformCountMatrix = counts_raw,
        isoformRepExpression = counts_raw,
        designMatrix = design_matrix,
        isoformExonAnnoation = "${gtf_file}",
        isoformNtFasta = NULL,
        showProgress = FALSE
    )

    # 6. RUN ANALYSIS STEPS MANUALLY

    # Step A: Pre-filter
    cat("Step A: Pre-filtering...\\n")
    switchList <- preFilter(
        switchAnalyzeRlist = switchList,
        geneExpressionCutoff = 1,
        isoformExpressionCutoff = 0,
        removeSingleIsoformGenes = TRUE,
        quiet = TRUE
    )

    # Step B: Test for Switches
    cat("Step B: Testing for Switches (DEXSeq)...\\n")
    switchList <- isoformSwitchTestDEXSeq(
        switchAnalyzeRlist = switchList,
        reduceToSwitchingGenes = TRUE,
        dIFcutoff = 0.1
    )

    # Step C: Analyse ORFs
    cat("Step C: Analysing ORFs...\\n")
    switchList <- analyzeORF(
        switchAnalyzeRlist = switchList,
        genomeObject = Hsapiens
    )

    # Step D: Extract Sequences
    cat("Step D: Extracting Sequences...\\n")
    switchList <- extractSequence(
        switchAnalyzeRlist = switchList,
        genomeObject = Hsapiens,
        dIFcutoff = 0.1,
        pathToOutput = ".",
        writeToFile = TRUE,             
        outputPrefix = "isoform_switch", 
        removeLongAAseq = FALSE,         
        alsoSplitFastaFile = FALSE       
    )

    # 7. SAVE OUTPUTS
    cat("Saving results...\\n")
    saveRDS(switchList, "switch_list_part1.rds")

    sink("part1_summary.txt")
    print(summary(switchList))
    print(extractSwitchSummary(switchList, dIFcutoff = 0.1))
    sink()
    """
}
