#!/usr/bin/env nextflow
/*
========================================================================================
    SUBWORKFLOW: PERFORM DIFFERENTIAL EXPRESSION (LONG READS)
========================================================================================
    1. Statistical analysis (DESeq2)
    2. Visualisation (Volcano plots, PCA)
    3. PCA eigenvector analysis
    4. Gene set enrichment analysis (GSEA)
========================================================================================
*/
nextflow.enable.dsl = 2

// Import modules
include { DESEQ2_LONG } from '../modules/deseq2_long'
include { DESEQ2_SUMMARY_LONG } from '../modules/deseq2_summary_long'
include { ENHANCED_VOLCANO_LONG } from '../modules/enhanced_volcano_long'
include { PCA_EIGENVECTORS_LONG } from '../modules/pca_eigenvectors_long'
include { PCA_EIGENVECTORS_SUMMARY_LONG } from '../modules/pca_eigenvectors_summary_long'
include { GSEA_LONG } from '../modules/gsea_long'

workflow PERFORM_DIFFERENTIAL_EXPRESSION_LONG {
    take:
    counts        // path: featureCounts output
    sample_info   // path: sample metadata CSV

    main:
    log.info "=== Starting Differential Expression Analysis (Long Reads) ==="

    // Step 1: Statistical analysis
    log.info "Step 1: Running DESeq2 analysis"
    DESEQ2_LONG(counts, sample_info)
    DESEQ2_SUMMARY_LONG(DESEQ2_LONG.out.results)

    // Step 2: Visualisation
    log.info "Step 2: Generating volcano plots and PCA"
    ENHANCED_VOLCANO_LONG(DESEQ2_LONG.out.results)

    // Step 3: PCA eigenvector analysis
    log.info "Step 3: Analysing PCA drivers"
    PCA_EIGENVECTORS_LONG(
        DESEQ2_LONG.out.norm_counts,
        sample_info
    )
    PCA_EIGENVECTORS_SUMMARY_LONG(
        PCA_EIGENVECTORS_LONG.out.pca_loadings,
        PCA_EIGENVECTORS_LONG.out.pca_var,
        PCA_EIGENVECTORS_LONG.out.top_genes
    )

    // Step 4: Gene set enrichment analysis
    log.info "Step 4: Performing gene set enrichment analysis"
    GSEA_LONG(
        DESEQ2_LONG.out.ranked_genes,
        file(params.hallmark_gmt),
        file(params.reactome_gmt)
    )

    emit:
    results = DESEQ2_LONG.out.results
    norm_counts = DESEQ2_LONG.out.norm_counts
    
    // Summaries for reporting
    deseq2_summary = DESEQ2_SUMMARY_LONG.out.summary
    pca_plot = DESEQ2_LONG.out.pca_plot
    volcano_pdf = ENHANCED_VOLCANO_LONG.out.pdf
    pca_eigenvectors_summary = PCA_EIGENVECTORS_SUMMARY_LONG.out.summary
    gsea_hallmark = GSEA_LONG.out.hallmark_results
    gsea_reactome = GSEA_LONG.out.reactome_results
}