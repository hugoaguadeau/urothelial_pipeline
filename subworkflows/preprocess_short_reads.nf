#!/usr/bin/env nextflow
/*
========================================================================================
    SUBWORKFLOW: PREPROCESS SHORT READS
========================================================================================
    1. Raw QC
    2. Adapter trimming and filtering (Fastp)
    3. Post-trim QC
    4. Comparative analysis
========================================================================================
*/
nextflow.enable.dsl = 2

// Import modules
include { FASTQC_SHORT } from '../modules/fastqc_short'
include { FASTQC_SUMMARY_SHORT } from '../modules/fastqc_summary_short'
include { SEQKIT_STATS_SHORT } from '../modules/seqkit_short'
include { SEQKIT_SUMMARY_SHORT } from '../modules/seqkit_summary_short'
include { FASTP_SHORT } from '../modules/fastp_short'
include { FASTP_SUMMARY_SHORT } from '../modules/fastp_summary_short'
include { FASTQC_TRIMMED_SHORT } from '../modules/fastqc_trimmed_short'
include { FASTQC_SUMMARY_TRIMMED_SHORT } from '../modules/fastqc_summary_trimmed_short'
include { SEQKIT_TRIMMED_SHORT } from '../modules/seqkit_trimmed_short'
include { SEQKIT_SUMMARY_TRIMMED_SHORT } from '../modules/seqkit_summary_trimmed_short'
include { QC_COMPARISON_SHORT } from '../modules/qc_comparison_short'

workflow PREPROCESS_SHORT_READS {
    take:
    reads        // channel: tuple(sample_id, [read1, read2])
    sample_info  // path: samples.csv

    main:
    log.info "=== Starting Pre-processing for Short Reads ==="

    // Step 1: Raw QC
    log.info "Step 1: Characterising raw short reads"
    FASTQC_SHORT(reads)
    FASTQC_SUMMARY_SHORT(FASTQC_SHORT.out.zip.collect())
    
    SEQKIT_STATS_SHORT(reads)
    SEQKIT_SUMMARY_SHORT(SEQKIT_STATS_SHORT.out.stats.collect())

    // Step 2: Trimming
    log.info "Step 2: Filtering and trimming (Fastp)"
    FASTP_SHORT(reads)
    FASTP_SUMMARY_SHORT(FASTP_SHORT.out.json_report.collect())

    // Step 3: Post-Trim QC
    log.info "Step 3: Validating processed read quality"
    SEQKIT_TRIMMED_SHORT(FASTP_SHORT.out.trimmed_reads)
    SEQKIT_SUMMARY_TRIMMED_SHORT(SEQKIT_TRIMMED_SHORT.out.stats.collect())
    
    FASTQC_TRIMMED_SHORT(FASTP_SHORT.out.trimmed_reads)
    FASTQC_SUMMARY_TRIMMED_SHORT(FASTQC_TRIMMED_SHORT.out.zip.collect())

    // Step 4: Comparative analysis
    log.info "Step 4: Comparing raw vs trimmed metrics"
    QC_COMPARISON_SHORT(
        FASTQC_SUMMARY_SHORT.out.summary,
        FASTQC_SUMMARY_TRIMMED_SHORT.out.summary,
        SEQKIT_SUMMARY_SHORT.out.summary,
        SEQKIT_SUMMARY_TRIMMED_SHORT.out.summary,
        FASTP_SUMMARY_SHORT.out.summary
    )

    emit:
    processed_reads = FASTP_SHORT.out.trimmed_reads
    
    // Summaries for reporting
    fastp_summary = FASTP_SUMMARY_SHORT.out.summary
    raw_fastqc_summary = FASTQC_SUMMARY_SHORT.out.summary
    trim_fastqc_summary = FASTQC_SUMMARY_TRIMMED_SHORT.out.summary
    raw_seqkit_summary = SEQKIT_SUMMARY_SHORT.out.summary
    trim_seqkit_summary = SEQKIT_SUMMARY_TRIMMED_SHORT.out.summary
    qc_comparison = QC_COMPARISON_SHORT.out.comparison_txt
}