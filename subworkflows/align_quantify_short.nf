#!/usr/bin/env nextflow
/*
========================================================================================
    SUBWORKFLOW: ALIGN AND QUANTIFY (SHORT READS)
========================================================================================
    1. Alignment (STAR)
    2. Gene quantification (featureCounts)
    3. Transcript assembly and quantification (StringTie)
========================================================================================
*/
nextflow.enable.dsl = 2

// Import modules
include { STAR_INDEX_SHORT } from '../modules/star_index_short'
include { STAR_ALIGN_SHORT } from '../modules/star_align_short'
include { STAR_ALIGNMENT_SUMMARY_SHORT } from '../modules/star_alignment_summary_short'
include { FEATURECOUNTS_SHORT } from '../modules/featurecounts_short'
include { FEATURECOUNTS_SUMMARY_SHORT } from '../modules/featurecounts_summary_short'
include { STRINGTIE_SHORT } from '../modules/stringtie_short'
include { STRINGTIE_MERGE_SHORT } from '../modules/stringtie_merge_short'
include { STRINGTIE_REQUANT_SHORT } from '../modules/stringtie_requant_short'
include { STRINGTIE_PREPDE_SHORT } from '../modules/stringtie_prepde_short'
include { STRINGTIE_SUMMARY_SHORT } from '../modules/stringtie_summary_short'

workflow ALIGN_QUANTIFY_SHORT {
    take:
    reads             // channel: tuple(sample_id, [read1, read2])
    reference_genome  // path: reference genome FASTA
    reference_gtf     // path: reference GTF annotation
    sample_info       // path: samples.csv

    main:
    log.info "=== Starting Alignment and Quantification (Short Reads) ==="

    // Step 1: Alignment
    log.info "Step 1: Building STAR index and aligning reads"
    STAR_INDEX_SHORT(reference_genome, reference_gtf)
    STAR_ALIGN_SHORT(
        reads,
        STAR_INDEX_SHORT.out.index,
        reference_gtf
    )
    STAR_ALIGNMENT_SUMMARY_SHORT(STAR_ALIGN_SHORT.out.log_final.collect())

    // Step 2: Gene quantification
    log.info "Step 2: Quantifying genes with featureCounts"
    FEATURECOUNTS_SHORT(
        STAR_ALIGN_SHORT.out.bam.map { sample_id, bam -> bam }.collect(),
        reference_gtf
    )
    FEATURECOUNTS_SUMMARY_SHORT(FEATURECOUNTS_SHORT.out.summary)

    // Step 3: Transcript quantification
    log.info "Step 3: Assembling and quantifying transcripts with StringTie"
    STRINGTIE_SHORT(
        STAR_ALIGN_SHORT.out.bam.join(STAR_ALIGN_SHORT.out.bam_index),
        reference_gtf
    )
    STRINGTIE_MERGE_SHORT(
        STRINGTIE_SHORT.out.gtf.map { id, gtf -> gtf }.collect(),
        reference_gtf
    )
    STRINGTIE_REQUANT_SHORT(
        STAR_ALIGN_SHORT.out.bam.join(STAR_ALIGN_SHORT.out.bam_index),
        STRINGTIE_MERGE_SHORT.out.merged_gtf
    )
    STRINGTIE_PREPDE_SHORT(
        STRINGTIE_REQUANT_SHORT.out.gtf.map { id, gtf -> gtf }.collect()
    )
    STRINGTIE_SUMMARY_SHORT(
        STRINGTIE_SHORT.out.gene_abundance.collect(),
        STRINGTIE_SHORT.out.transcript_abundance.collect(),
        STRINGTIE_PREPDE_SHORT.out.gene_counts,
        STRINGTIE_PREPDE_SHORT.out.transcript_counts
    )

    emit:
    gene_counts = FEATURECOUNTS_SHORT.out.counts
    transcript_counts = STRINGTIE_PREPDE_SHORT.out.transcript_counts
    merged_gtf = STRINGTIE_MERGE_SHORT.out.merged_gtf
    
    // Summaries for reporting
    star_summary = STAR_ALIGNMENT_SUMMARY_SHORT.out.summary_txt
    featurecounts_summary = FEATURECOUNTS_SUMMARY_SHORT.out.summary_txt
    stringtie_summary = STRINGTIE_SUMMARY_SHORT.out.summary_txt
}