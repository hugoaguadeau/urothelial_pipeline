#!/usr/bin/env nextflow

/*
========================================================================================
    SUBWORKFLOW: PREPROCESS LONG READS
========================================================================================
    1. Raw QC
    2. Full length read identification (Pychopper)
    3. Filtering and Trimming (Fastplong)
    4. Post-trim QC
========================================================================================
*/

nextflow.enable.dsl = 2

// Import modules
include { FASTQC_LONG } from '../modules/fastqc_long'
include { FASTQC_SUMMARY_LONG } from '../modules/fastqc_summary_long'
include { SEQKIT_STATS_LONG } from '../modules/seqkit_long'
include { SEQKIT_SUMMARY_LONG } from '../modules/seqkit_summary_long'
include { PYCHOPPER_LONG } from '../modules/pychopper_long'
include { PYCHOPPER_SUMMARY_LONG } from '../modules/pychopper_summary_long'
include { FASTPLONG_LONG } from '../modules/fastplong_long'
include { FASTPLONG_SUMMARY_LONG } from '../modules/fastplong_summary_long'
include { FASTQC_TRIMMED_LONG } from '../modules/fastqc_trimmed_long'
include { FASTQC_SUMMARY_TRIMMED_LONG } from '../modules/fastqc_summary_trimmed_long'
include { SEQKIT_TRIMMED_LONG } from '../modules/seqkit_trimmed_long'
include { SEQKIT_SUMMARY_TRIMMED_LONG } from '../modules/seqkit_summary_trimmed_long'

workflow PREPROCESS_LONG_READS {

    take:
    reads    // channel: tuple(sample_id, fastq)

    main:
    log.info "=== Starting Pre-processing for Long Reads ==="

    // Step 1: Raw QC
    log.info "Step 1: Characterising raw long reads"
    FASTQC_LONG(reads)
    FASTQC_SUMMARY_LONG(FASTQC_LONG.out.zip.collect())
    SEQKIT_STATS_LONG(reads)
    SEQKIT_SUMMARY_LONG(SEQKIT_STATS_LONG.out.stats.collect())

    // Step 2: Pychopper (Full length identification)
    log.info "Step 2: Identifying full-length reads (Pychopper)"
    PYCHOPPER_LONG(reads)
    PYCHOPPER_SUMMARY_LONG(PYCHOPPER_LONG.out.stats.collect())

    // Step 3: Trimming
    log.info "Step 3: Filtering and trimming (Fastplong)"
    FASTPLONG_LONG(PYCHOPPER_LONG.out.full_length_reads)
    FASTPLONG_SUMMARY_LONG(FASTPLONG_LONG.out.json_report.collect())

    // Step 4: Post-Trim QC
    log.info "Step 4: Validating processed read quality"
    SEQKIT_TRIMMED_LONG(FASTPLONG_LONG.out.trimmed_reads)
    SEQKIT_SUMMARY_TRIMMED_LONG(SEQKIT_TRIMMED_LONG.out.stats.collect())
    
    FASTQC_TRIMMED_LONG(FASTPLONG_LONG.out.trimmed_reads)
    FASTQC_SUMMARY_TRIMMED_LONG(FASTQC_TRIMMED_LONG.out.zip.collect())

    emit:
    processed_reads = FASTPLONG_LONG.out.trimmed_reads
    
    // Summaries for reporting
    pychopper_summary = PYCHOPPER_SUMMARY_LONG.out.summary
    fastplong_summary = FASTPLONG_SUMMARY_LONG.out.summary
    raw_fastqc_summary = FASTQC_SUMMARY_LONG.out.summary
    trim_fastqc_summary = FASTQC_SUMMARY_TRIMMED_LONG.out.summary
    raw_seqkit_summary = SEQKIT_SUMMARY_LONG.out.summary
    trim_seqkit_summary = SEQKIT_SUMMARY_TRIMMED_LONG.out.summary
}
