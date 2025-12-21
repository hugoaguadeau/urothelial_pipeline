#!/usr/bin/env nextflow
/*
========================================================================================
    SUBWORKFLOW: BLADDER CANCER SUBTYPING (LONG READS)
========================================================================================
    1. Consensus classifier
    2. Subtyping summary
========================================================================================
*/
nextflow.enable.dsl = 2

// Import modules
include { CONSENSUS_CLASSIFIER_LONG } from '../modules/consensus_classifier_long'
include { SUBTYPING_SUMMARY_LONG } from '../modules/subtyping_summary_long'

workflow SUBTYPE_BLADDER_CANCER_LONG {
    take:
    tpm_counts    // path: normalised TPM counts
    sample_info   // path: sample metadata CSV

    main:
    log.info "=== Starting Bladder Cancer Molecular Subtyping ==="

    // Step 1: Consensus classifier
    log.info "Step 1: Running consensus MIBC classifier"
    CONSENSUS_CLASSIFIER_LONG(tpm_counts)

    // Step 2: Subtyping summary
    log.info "Step 2: Generating subtyping summary"
    SUBTYPING_SUMMARY_LONG(
        CONSENSUS_CLASSIFIER_LONG.out.subtypes,
        sample_info
    )

    emit:
    subtypes = CONSENSUS_CLASSIFIER_LONG.out.subtypes
    
    // Summaries for reporting
    subtyping_summary = SUBTYPING_SUMMARY_LONG.out.summary_txt
}