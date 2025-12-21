#!/usr/bin/env nextflow
/*
========================================================================================
    SUBWORKFLOW: COMPARE TECHNOLOGIES
========================================================================================
    1. Pi-value correlation analysis
========================================================================================
*/
nextflow.enable.dsl = 2

// Import modules
include { PI_VALUE_COMPARISON } from '../modules/pi_value_comparison'

workflow COMPARE_TECHNOLOGIES {
    take:
    short_read_deseq2   // path: DESeq2 results (short)
    long_read_deseq2    // path: DESeq2 results (long)

    main:
    log.info "=== Starting Technology Comparison ==="

    // Step 1: Pi-value correlation
    log.info "Step 1: Calculating pi-value correlation"
    PI_VALUE_COMPARISON(short_read_deseq2, long_read_deseq2)

    emit:
    correlation_plot = PI_VALUE_COMPARISON.out.scatter_plot
}