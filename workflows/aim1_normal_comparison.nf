#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    WORKFLOW: AIM 1 - NORMAL TISSUE COMPARISON
========================================================================================
    Compares normal bladder vs normal ureter.
    Supports both Short and Long-read technologies.
========================================================================================
*/

// --- IMPORTS: SHORT READS ---
include { PREPROCESS_SHORT_READS } from '../subworkflows/preprocess_short_reads'
include { ALIGN_QUANTIFY_SHORT } from '../subworkflows/align_quantify_short'
include { PERFORM_DIFFERENTIAL_EXPRESSION_SHORT } from '../subworkflows/perform_differential_expression_short'
include { ANALYSE_BIOLOGICAL_PATHWAYS_SHORT } from '../subworkflows/analyse_biological_pathways_short'
include { ANALYSE_ISOFORMS_SHORT } from '../subworkflows/analyse_isoforms_short'

// --- IMPORTS: LONG READS ---
include { PREPROCESS_LONG_READS } from '../subworkflows/preprocess_long_reads'
include { ALIGN_QUANTIFY_LONG } from '../subworkflows/align_quantify_long'
include { PERFORM_DIFFERENTIAL_EXPRESSION_LONG } from '../subworkflows/perform_differential_expression_long'
include { ANALYSE_BIOLOGICAL_PATHWAYS_LONG } from '../subworkflows/analyse_biological_pathways_long'
include { ANALYSE_ISOFORMS_LONG } from '../subworkflows/analyse_isoforms_long'

// --- IMPORTS: COMPARATIVE ---
include { COMPARE_TECHNOLOGIES } from '../subworkflows/compare_technologies'

workflow AIM1_NORMAL_COMPARISON {

    log.info "=== AIM 1: NORMAL TISSUE COMPARISON ==="

    // 1. Parse input channel
    samples_branched = Channel.fromPath(params.input).splitCsv(header:true)
        .filter { it.aim1_group != null && it.aim1_group != 'NA' && it.aim1_group != '' }
        .branch {
            short_read: it.seq_type == 'short'
                return tuple(it.sample_id, [file(it.fastq_1), file(it.fastq_2)])
            long_read: it.seq_type == 'long'
                return tuple(it.sample_id, file(it.fastq_1))
        }

    // 2. Load reference files from params
    reference_genome = file(params.reference_genome, checkIfExists: true)
    reference_gtf    = file(params.reference_gtf, checkIfExists: true)
    sample_info      = file(params.input, checkIfExists: true)

    // 3. Load functional annotation files (for isoform analysis)
    hexamer = file(params.cpat_hexamer)
    logit   = file(params.cpat_logit)
    pfam_db = file(params.pfam_db)

    // ========================================================================
    // SHORT READ ANALYSIS
    // ========================================================================
    if (params.run_short_read) {
        log.info "=== Starting Short Read Analysis ==="

        // Step 1: Pre-processing
        PREPROCESS_SHORT_READS(
            samples_branched.short_read,
            sample_info
        )

        // Step 2: Alignment and quantification
        ALIGN_QUANTIFY_SHORT(
            PREPROCESS_SHORT_READS.out.trimmed_reads,
            reference_genome,
            reference_gtf,
            sample_info
        )

        // Step 3: Differential expression
        PERFORM_DIFFERENTIAL_EXPRESSION_SHORT(
            ALIGN_QUANTIFY_SHORT.out.gene_counts,
            sample_info
        )

        // Step 4: Biological pathway analysis
        ANALYSE_BIOLOGICAL_PATHWAYS_SHORT(
            ALIGN_QUANTIFY_SHORT.out.gene_counts,
            sample_info,
            "aim1_group"
        )

        // Step 5: Isoform analysis
        ANALYSE_ISOFORMS_SHORT(
            ALIGN_QUANTIFY_SHORT.out.transcript_counts,
            sample_info,
            ALIGN_QUANTIFY_SHORT.out.merged_gtf,
            reference_genome,
            hexamer,
            logit,
            pfam_db
        )
    }

    // ========================================================================
    // LONG READ ANALYSIS
    // ========================================================================
    if (params.run_long_read) {
        log.info "=== Starting Long Read Analysis ==="

        // Step 1: Pre-processing
        PREPROCESS_LONG_READS(samples_branched.long_read)

        // Step 2: Alignment and quantification
        ALIGN_QUANTIFY_LONG(
            PREPROCESS_LONG_READS.out.processed_reads,
            reference_genome,
            reference_gtf
        )

        // Step 3: Differential expression
        PERFORM_DIFFERENTIAL_EXPRESSION_LONG(
            ALIGN_QUANTIFY_LONG.out.gene_counts,
            sample_info
        )

        // Step 4: Biological pathway analysis
        ANALYSE_BIOLOGICAL_PATHWAYS_LONG(
            ALIGN_QUANTIFY_LONG.out.gene_counts,
            sample_info,
            "aim1_group"
        )

        // Step 5: Isoforms
        ANALYSE_ISOFORMS_LONG(
            ALIGN_QUANTIFY_LONG.out.transcript_counts,
            sample_info,
            ALIGN_QUANTIFY_LONG.out.merged_gtf,
            reference_genome,
            hexamer,         
            logit,           
            pfam_db          
        )
    }

    // ========================================================================
    // COMPARATIVE ANALYSIS
    // ========================================================================
    if (params.run_short_read && params.run_long_read) {
        log.info "=== Starting Comparative Analysis ==="

        COMPARE_TECHNOLOGIES(
            PERFORM_DIFFERENTIAL_EXPRESSION_SHORT.out.results,
            PERFORM_DIFFERENTIAL_EXPRESSION_LONG.out.results
        )
    }
}
