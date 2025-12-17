#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    WORKFLOW: AIM 2 - CANCER TISSUE COMPARISON
========================================================================================
    Comparing bladder vs upper tract UC using long-read.
========================================================================================
*/

// --- IMPORTS ---
include { PREPROCESS_LONG_READS } from '../subworkflows/preprocess_long_reads'
include { ALIGN_QUANTIFY_LONG } from '../subworkflows/align_quantify_long'
include { PERFORM_DIFFERENTIAL_EXPRESSION_LONG } from '../subworkflows/perform_differential_expression_long'
include { ANALYSE_BIOLOGICAL_PATHWAYS_LONG } from '../subworkflows/analyse_biological_pathways_long'
include { SUBTYPE_BLADDER_CANCER_LONG } from '../subworkflows/subtype_bladder_cancer_long'
include { ANALYSE_ISOFORMS_LONG } from '../subworkflows/analyse_isoforms_long'

workflow AIM2_CANCER_COMPARISON {

    log.info "=== AIM 2: CANCER TISSUE COMPARISON ==="

    // 1. Parse input channel
    long_reads = Channel.fromPath(params.input).splitCsv(header: true)
        .filter { row ->
            row.aim2_group != null && row.aim2_group != 'NA' &&
            row.aim2_group != '' && row.seq_type == 'long'
        }
        .map { row -> tuple(row.sample_id, file(row.fastq_1)) }

    // 2. Load reference files from params
    reference_genome = file(params.reference_genome, checkIfExists: true)
    reference_gtf    = file(params.reference_gtf, checkIfExists: true)
    sample_info      = file(params.input, checkIfExists: true)

    // 3. Load functional annotation files (for isoform analysis)
    hexamer = file(params.cpat_hexamer)
    logit   = file(params.cpat_logit)
    pfam_db = file(params.pfam_db)

    // ========================================================================
    // LONG READ ANALYSIS
    // ========================================================================
    log.info "=== Starting Long Read Analysis ==="

    // Step 1: Pre-processing
    PREPROCESS_LONG_READS(long_reads)

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
        "aim2_group"
    )

    // Step 4b: Bladder cancer molecular subtyping
    SUBTYPE_BLADDER_CANCER_LONG(
        ANALYSE_BIOLOGICAL_PATHWAYS_LONG.out.tpm_counts,
        sample_info
    )

    // Step 5: Isoform analysis
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

