#!/usr/bin/env nextflow
/*
========================================================================================
    SUBWORKFLOW: ALIGN AND QUANTIFY (LONG READS)
========================================================================================
    1. Alignment (Minimap2)
    2. Gene quantification (featureCounts)
    3. Transcript assembly and quantification (StringTie)
========================================================================================
*/
nextflow.enable.dsl = 2

// Import modules
include { MINIMAP2_LONG } from '../modules/minimap2_long'
include { MINIMAP2_SUMMARY_LONG } from '../modules/minimap2_summary_long'
include { FEATURECOUNTS_LONG } from '../modules/featurecounts_long'
include { FEATURECOUNTS_SUMMARY_LONG } from '../modules/featurecounts_summary_long'
include { STRINGTIE_LONG } from '../modules/stringtie_long'
include { STRINGTIE_MERGE_LONG } from '../modules/stringtie_merge_long'
include { STRINGTIE_REQUANT_LONG } from '../modules/stringtie_requant_long'
include { STRINGTIE_PREPDE_LONG } from '../modules/stringtie_prepde_long'
include { STRINGTIE_SUMMARY_LONG } from '../modules/stringtie_summary_long'

workflow ALIGN_QUANTIFY_LONG {
    take:
    reads              // channel: tuple(sample_id, fastq)
    reference_genome   // path: genome fasta
    reference_gtf      // path: gene annotation

    main:
    log.info "=== Starting Alignment and Quantification (Long Reads) ==="

    // Step 1: Alignment
    log.info "Step 1: Aligning reads with Minimap2"
    MINIMAP2_LONG(reads, reference_genome, reference_gtf)
    MINIMAP2_SUMMARY_LONG(MINIMAP2_LONG.out.bam.map { it[1] }.collect())

    // Step 2: Gene quantification
    log.info "Step 2: Quantifying genes with featureCounts"
    FEATURECOUNTS_LONG(
        MINIMAP2_LONG.out.bam.map { it[1] }.collect(),
        reference_gtf
    )
    FEATURECOUNTS_SUMMARY_LONG(FEATURECOUNTS_LONG.out.summary)

    // Step 3: Transcript quantification
    log.info "Step 3: Assembling and quantifying transcripts with StringTie"
    STRINGTIE_LONG(MINIMAP2_LONG.out.bam, reference_gtf)
    STRINGTIE_MERGE_LONG(
        STRINGTIE_LONG.out.gtf_for_merge.collect(),
        reference_gtf
    )
    STRINGTIE_REQUANT_LONG(
        MINIMAP2_LONG.out.bam,
        STRINGTIE_MERGE_LONG.out.merged_gtf
    )
    STRINGTIE_PREPDE_LONG(
        STRINGTIE_REQUANT_LONG.out.gtf_for_prepde.collect()
    )
    STRINGTIE_SUMMARY_LONG(
        STRINGTIE_PREPDE_LONG.out.gene_counts,
        STRINGTIE_PREPDE_LONG.out.transcript_counts
    )

    emit:
    bams = MINIMAP2_LONG.out.bam
    gene_counts = FEATURECOUNTS_LONG.out.counts
    transcript_counts = STRINGTIE_PREPDE_LONG.out.transcript_counts
    merged_gtf = STRINGTIE_MERGE_LONG.out.merged_gtf
    
    // Summaries for reporting
    minimap_summary = MINIMAP2_SUMMARY_LONG.out.summary_txt
    fc_summary = FEATURECOUNTS_SUMMARY_LONG.out.summary_txt
    stringtie_summary = STRINGTIE_SUMMARY_LONG.out.summary_txt
}