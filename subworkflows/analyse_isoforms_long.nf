#!/usr/bin/env nextflow
/*
========================================================================================
    SUBWORKFLOW: ANALYSE ISOFORMS (LONG READS)
========================================================================================
    1. Isoform switch detection and sequence extraction
    2. External sequence analysis (CPAT, Pfam)
    3. Integration and visualisation
========================================================================================
*/
nextflow.enable.dsl = 2

// Import modules
include { ISOFORM_SWITCH_PART1_LONG } from '../modules/isoform_switch_part1_long'
include { CPAT_LONG } from '../modules/cpat_long'
include { PFAM_LONG } from '../modules/pfam_long'
include { ISOFORM_SWITCH_PART2_LONG } from '../modules/isoform_switch_part2_long'

workflow ANALYSE_ISOFORMS_LONG {
    take:
    transcript_counts // path: StringTie/Bambu count matrix
    sample_info       // path: samples.csv
    gtf_file          // path: merged GTF
    genome_fasta      // path: genome fasta
    hexamer_file      // path: CPAT hexamer
    logit_model       // path: CPAT logit model
    pfam_db           // path: Pfam database

    main:
    log.info "=== Starting Isoform Switch Analysis (Long Reads) ==="

    // Step 1: Isoform switch detection
    log.info "Step 1: Detecting isoform switches and extracting sequences"
    ISOFORM_SWITCH_PART1_LONG(
        transcript_counts,
        sample_info,
        gtf_file,
        genome_fasta
    )

    // Step 2: External sequence analysis
    log.info "Step 2: Running external sequence analysis"
    CPAT_LONG(
        ISOFORM_SWITCH_PART1_LONG.out.nt_fasta,
        hexamer_file,
        logit_model
    )
    PFAM_LONG(
        ISOFORM_SWITCH_PART1_LONG.out.aa_fasta,
        pfam_db
    )

    // Step 3: Integration and visualisation
    log.info "Step 3: Integrating results and generating visualisations"
    ISOFORM_SWITCH_PART2_LONG(
        ISOFORM_SWITCH_PART1_LONG.out.switch_list_rds,
        CPAT_LONG.out.results,
        PFAM_LONG.out.results
    )

    emit:
    switch_list = ISOFORM_SWITCH_PART2_LONG.out.final_rds
    
    // Summaries for reporting
    plots = ISOFORM_SWITCH_PART2_LONG.out.plots
    summary = ISOFORM_SWITCH_PART2_LONG.out.summary
}