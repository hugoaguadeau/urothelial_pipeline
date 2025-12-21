#!/usr/bin/env nextflow

/*
========================================================================================
    SUBWORKFLOW: BIOLOGICAL ANALYSIS (LONG READS)
========================================================================================
    1. TPM normalisation
    2. HOX gene expression
    3. Luminal/basal phenotype
    4. Urothelial differentiation
    5. Signalling pathways (SHH, FGF, WNT, TGF-beta)
    6. Immune response and markers
    7. Proliferation and cell cycle
========================================================================================
*/

nextflow.enable.dsl = 2

// Import modules
include { TPM_NORMALISATION } from '../modules/tpm_normalisation'
include { HOX_GENE_ANALYSIS_LONG } from '../modules/hox_analysis_long'
include { LUMINAL_BASAL_ANALYSIS_LONG } from '../modules/luminal_basal_analysis_long'
include { UROTHELIAL_DIFFERENTIATION_ANALYSIS_LONG } from '../modules/urothelial_differentiation_long'
include { SHH_SIGNALLING_ANALYSIS_LONG } from '../modules/shh_signalling_long'
include { FGF_SIGNALLING_ANALYSIS_LONG } from '../modules/fgf_signalling_long'
include { WNT_BETACATENIN_SIGNALLING_ANALYSIS_LONG } from '../modules/wnt_betacatenin_signalling_long'
include { TGFBETA_SIGNALLING_ANALYSIS_LONG } from '../modules/tgfbeta_signalling_long'
include { INTERFERON_ANALYSIS_LONG } from '../modules/interferon_analysis_long'
include { IMMUNE_MARKERS_ANALYSIS_LONG } from '../modules/immune_markers_long'
include { PROLIFERATION_ANALYSIS_LONG } from '../modules/proliferation_analysis_long'

workflow ANALYSE_BIOLOGICAL_PATHWAYS_LONG {
    take:
    counts         // path: featureCounts output
    sample_info    // path: sample metadata CSV
    analysis_group // val: column in sample_info to filter by

    main:
    log.info "=== Starting Biological Pathway Analysis (Long Reads) ==="

    // Step 1: TPM normalisation
    log.info "Step 1: Normalising counts to TPM"
    TPM_NORMALISATION(counts, "long_read")

    // Step 2: HOX gene expression
    log.info "Step 2: Analysing HOX gene expression"
    HOX_GENE_ANALYSIS_LONG(TPM_NORMALISATION.out.tpm_counts, sample_info, analysis_group)

    // Step 3: Luminal/basal phenotype
    log.info "Step 3: Analysing luminal/basal phenotype"
    LUMINAL_BASAL_ANALYSIS_LONG(TPM_NORMALISATION.out.tpm_counts, sample_info, analysis_group)

    // Step 4: Urothelial differentiation
    log.info "Step 4: Analysing urothelial differentiation"
    UROTHELIAL_DIFFERENTIATION_ANALYSIS_LONG(TPM_NORMALISATION.out.tpm_counts, sample_info, analysis_group)

    // Step 5: SHH signalling
    log.info "Step 5: Analysing SHH signalling pathway"
    SHH_SIGNALLING_ANALYSIS_LONG(TPM_NORMALISATION.out.tpm_counts, sample_info, analysis_group)

    // Step 6: FGF signalling
    log.info "Step 6: Analysing FGF signalling pathway"
    FGF_SIGNALLING_ANALYSIS_LONG(TPM_NORMALISATION.out.tpm_counts, sample_info, analysis_group)

    // Step 7: WNT signalling
    log.info "Step 7: Analysing WNT beta-catenin signalling pathway"
    WNT_BETACATENIN_SIGNALLING_ANALYSIS_LONG(TPM_NORMALISATION.out.tpm_counts, sample_info, analysis_group)

    // Step 8: TGF-beta signalling
    log.info "Step 8: Analysing TGF-beta signalling pathway"
    TGFBETA_SIGNALLING_ANALYSIS_LONG(TPM_NORMALISATION.out.tpm_counts, sample_info, analysis_group)

    // Step 9: Interferon response
    log.info "Step 9: Analysing interferon response"
    INTERFERON_ANALYSIS_LONG(TPM_NORMALISATION.out.tpm_counts, sample_info, analysis_group)

    // Step 10: Immune markers
    log.info "Step 10: Analysing immune cell markers"
    IMMUNE_MARKERS_ANALYSIS_LONG(TPM_NORMALISATION.out.tpm_counts, sample_info, analysis_group)

    // Step 11: Proliferation
    log.info "Step 11: Analysing proliferation and cell cycle"
    PROLIFERATION_ANALYSIS_LONG(TPM_NORMALISATION.out.tpm_counts, sample_info, analysis_group)

    emit:
    tpm_counts = TPM_NORMALISATION.out.tpm_counts
    
    // Summaries for reporting
    tpm_summary = TPM_NORMALISATION.out.summary
    hox_plots = HOX_GENE_ANALYSIS_LONG.out.plots_pdf
    hox_data = HOX_GENE_ANALYSIS_LONG.out.data_csv
    hox_summary = HOX_GENE_ANALYSIS_LONG.out.summary_txt
    luminal_basal_plots = LUMINAL_BASAL_ANALYSIS_LONG.out.plots_pdf
    luminal_basal_data = LUMINAL_BASAL_ANALYSIS_LONG.out.data_csv
    luminal_basal_summary = LUMINAL_BASAL_ANALYSIS_LONG.out.summary_txt
    urothelial_plots = UROTHELIAL_DIFFERENTIATION_ANALYSIS_LONG.out.plots_pdf
    urothelial_data = UROTHELIAL_DIFFERENTIATION_ANALYSIS_LONG.out.data_csv
    urothelial_summary = UROTHELIAL_DIFFERENTIATION_ANALYSIS_LONG.out.summary_txt
    shh_plots = SHH_SIGNALLING_ANALYSIS_LONG.out.plots_pdf
    shh_data = SHH_SIGNALLING_ANALYSIS_LONG.out.data_csv
    shh_summary = SHH_SIGNALLING_ANALYSIS_LONG.out.summary_txt
    fgf_plots = FGF_SIGNALLING_ANALYSIS_LONG.out.plots_pdf
    fgf_data = FGF_SIGNALLING_ANALYSIS_LONG.out.data_csv
    fgf_summary = FGF_SIGNALLING_ANALYSIS_LONG.out.summary_txt
    wnt_plots = WNT_BETACATENIN_SIGNALLING_ANALYSIS_LONG.out.plots_pdf
    wnt_data = WNT_BETACATENIN_SIGNALLING_ANALYSIS_LONG.out.data_csv
    wnt_summary = WNT_BETACATENIN_SIGNALLING_ANALYSIS_LONG.out.summary_txt
    tgfbeta_plots = TGFBETA_SIGNALLING_ANALYSIS_LONG.out.plots_pdf
    tgfbeta_data = TGFBETA_SIGNALLING_ANALYSIS_LONG.out.data_csv
    tgfbeta_summary = TGFBETA_SIGNALLING_ANALYSIS_LONG.out.summary_txt
    interferon_plots = INTERFERON_ANALYSIS_LONG.out.plots_pdf
    interferon_data = INTERFERON_ANALYSIS_LONG.out.data_csv
    interferon_summary = INTERFERON_ANALYSIS_LONG.out.summary_txt
    immune_plots = IMMUNE_MARKERS_ANALYSIS_LONG.out.plots_pdf
    immune_data = IMMUNE_MARKERS_ANALYSIS_LONG.out.data_csv
    immune_summary = IMMUNE_MARKERS_ANALYSIS_LONG.out.summary_txt
    proliferation_plots = PROLIFERATION_ANALYSIS_LONG.out.plots_pdf
    proliferation_data = PROLIFERATION_ANALYSIS_LONG.out.data_csv
    proliferation_summary = PROLIFERATION_ANALYSIS_LONG.out.summary_txt
}