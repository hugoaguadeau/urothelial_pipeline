#!/usr/bin/env nextflow

/*
========================================================================================
    UROTHELIAL CARCINOMA RNA-SEQ PIPELINE
========================================================================================
*/

nextflow.enable.dsl = 2

// Print pipeline header
log.info """
================================================================================
    UROTHELIAL CARCINOMA RNA-SEQ PIPELINE v2
================================================================================
    Analysis mode: ${params.run_aim1 ? 'Aim 1 ' : ''}${params.run_aim2 ? 'Aim 2' : ''}
    Short reads: ${params.run_short_read ? 'YES' : 'NO'}
    Long reads: ${params.run_long_read ? 'YES' : 'NO'}
    Output directory: ${params.output_dir}
================================================================================
"""

/*
========================================================================================
    IMPORT WORKFLOWS
========================================================================================
*/

include { AIM1_NORMAL_COMPARISON } from './workflows/aim1_normal_comparison'
include { AIM2_CANCER_COMPARISON } from './workflows/aim2_cancer_comparison'

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
    This is the dispatcher - it calls the appropriate workflow based on parameters
*/

workflow {

    // Validate that at least one aim is selected
    if (!params.run_aim1 && !params.run_aim2) {
        error "ERROR: No research aim selected. Please enable at least one aim in your config or profile."
    }

    // Validate that at least one technology is selected
    if (!params.run_short_read && !params.run_long_read) {
        error "ERROR: No sequencing technology selected. Please enable short_read and/or long_read."
    }

    // Check input file exists
    if (!file(params.input).exists()) {
        error "ERROR: Sample sheet not found: ${params.input}"
    }

    // Check reference files exist
    if (!file(params.reference_genome).exists()) {
        error "ERROR: Reference genome not found: ${params.reference_genome}"
    }
    if (!file(params.reference_gtf).exists()) {
        error "ERROR: Reference GTF not found: ${params.reference_gtf}"
    }

    // ============================================================================
    // RUN WORKFLOWS BASED ON SELECTED AIMS
    // ============================================================================

    if (params.run_aim1) {
        log.info "=== Running Aim 1: Normal Tissue Comparison (Bladder vs Ureter) ==="
        AIM1_NORMAL_COMPARISON()
    }

    if (params.run_aim2) {
        log.info "=== Running Aim 2: Cancer Tissue Comparison (Bladder Cancer vs Ureter Cancer) ==="
        AIM2_CANCER_COMPARISON()
    }
}

/*
========================================================================================
    COMPLETION HANDLERS
========================================================================================
*/

workflow.onComplete {
    log.info """
    ================================================================================
    Pipeline completed at: ${workflow.complete}
    Execution status: ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration: ${workflow.duration}
    ================================================================================
    """.stripIndent()
}

workflow.onError {
    log.error """
    ================================================================================
    Pipeline execution failed!
    Error message: ${workflow.errorMessage}
    ================================================================================
    """.stripIndent()
}

