#!/usr/bin/env nextflow
/*
========================================================================================
    SUBWORKFLOW: PERFORM DIFFERENTIAL EXPRESSION (SHORT READS)
========================================================================================
    1. Statistical analysis (DESeq2)
    2. Visualisation (Volcano plots, PCA)
    3. PCA eigenvector analysis
    4. Gene set enrichment analysis (GSEA)
========================================================================================
*/
nextflow.enable.dsl = 2

// Import modules
include { DESEQ2_SHORT } from '../modules/deseq2_short'
include { DESEQ2_SUMMARY_SHORT } from '../modules/deseq2_summary_short'
include { ENHANCED_VOLCANO_SHORT } from '../modules/enhanced_volcano_short'
include { PCA_EIGENVECTORS_SHORT } from '../modules/pca_eigenvectors_short'
include { PCA_EIGENVECTORS_SUMMARY_SHORT } from '../modules/pca_eigenvectors_summary_short'
include { GSEA_SHORT } from '../modules/gsea_short'

workflow PERFORM_DIFFERENTIAL_EXPRESSION_SHORT {
    take:
    counts        // path: featureCounts output
    sample_info   // path: sample metadata CSV

    main:
    log.info "=== Starting Differential Expression Analysis (Short Reads) ==="

    // Step 1: Statistical analysis
    log.info "Step 1: Running DESeq2 analysis"
    DESEQ2_SHORT(counts, sample_info)
    DESEQ2_SUMMARY_SHORT(DESEQ2_SHORT.out.results)

    // Step 2: Visualisation
    log.info "Step 2: Generating volcano plots and PCA"
    ENHANCED_VOLCANO_SHORT(DESEQ2_SHORT.out.results)

    // Step 3: PCA eigenvector analysis
    log.info "Step 3: Analysing PCA drivers"
    PCA_EIGENVECTORS_SHORT(
        DESEQ2_SHORT.out.norm_counts,
        sample_info
    )
    PCA_EIGENVECTORS_SUMMARY_SHORT(
        PCA_EIGENVECTORS_SHORT.out.pca_loadings,
        PCA_EIGENVECTORS_SHORT.out.pca_var,
        PCA_EIGENVECTORS_SHORT.out.top_genes
    )

    // Step 4: Gene set enrichment analysis
    log.info "Step 4: Performing gene set enrichment analysis"
    GSEA_SHORT(
        DESEQ2_SHORT.out.ranked_genes,
        file(params.hallmark_gmt),
        file(params.reactome_gmt)
    )

    emit:
    results = DESEQ2_SHORT.out.results
    norm_counts = DESEQ2_SHORT.out.norm_counts
    
    // Summaries for reporting
    deseq2_summary = DESEQ2_SUMMARY_SHORT.out.summary
    pca_plot = DESEQ2_SHORT.out.pca_plot
    volcano_pdf = ENHANCED_VOLCANO_SHORT.out.pdf
    pca_eigenvectors_summary = PCA_EIGENVECTORS_SUMMARY_SHORT.out.summary
    gsea_hallmark = GSEA_SHORT.out.hallmark_results
    gsea_reactome = GSEA_SHORT.out.reactome_results
}