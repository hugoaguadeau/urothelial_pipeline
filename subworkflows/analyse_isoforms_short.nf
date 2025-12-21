#!/usr/bin/env nextflow
/*
========================================================================================
    SUBWORKFLOW: ANALYSE ISOFORMS (SHORT READS)
========================================================================================
    1. Isoform switch detection and sequence extraction
    2. External sequence analysis (CPAT, Pfam)
    3. Integration and visualisation
========================================================================================
*/
nextflow.enable.dsl = 2

// Import modules
include { ISOFORM_SWITCH_PART1 } from '../modules/isoform_switch_part1'
include { CPAT_SHORT } from '../modules/cpat_short'
include { PFAM_SHORT } from '../modules/pfam_short'
include { ISOFORM_SWITCH_PART2 } from '../modules/isoform_switch_part2'

workflow ANALYSE_ISOFORMS_SHORT {
    take:
    transcript_counts // path: StringTie counts
    sample_info       // path: samples.csv
    gtf_file          // path: merged GTF
    genome_fasta      // path: genome fasta
    hexamer_file      // path: CPAT hexamer
    logit_model       // path: CPAT logit model
    pfam_db           // path: Pfam database

    main:
    log.info "=== Starting Isoform Switch Analysis ==="

    // Step 1: Isoform switch detection
    log.info "Step 1: Detecting isoform switches and extracting sequences"
    ISOFORM_SWITCH_PART1(
        transcript_counts,
        sample_info,
        gtf_file,
        genome_fasta
    )

    // Step 2: External sequence analysis
    log.info "Step 2: Running external sequence analysis"
    CPAT_SHORT(
        ISOFORM_SWITCH_PART1.out.nt_fasta,
        hexamer_file,
        logit_model
    )
    PFAM_SHORT(
        ISOFORM_SWITCH_PART1.out.aa_fasta,
        pfam_db
    )

    // Step 3: Integration and visualisation
    log.info "Step 3: Integrating results and generating visualisations"
    ISOFORM_SWITCH_PART2(
        ISOFORM_SWITCH_PART1.out.switch_list_rds,
        CPAT_SHORT.out.results,
        PFAM_SHORT.out.results
    )

    emit:
    switch_list = ISOFORM_SWITCH_PART2.out.final_rds
    
    // Summaries for reporting
    plots = ISOFORM_SWITCH_PART2.out.plots
    summary = ISOFORM_SWITCH_PART2.out.summary
}