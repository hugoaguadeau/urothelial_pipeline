#!/usr/bin/env nextflow
/*
========================================================================================
    SUBWORKFLOW: BIOLOGICAL ANALYSIS (SHORT READS)
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
include { HOX_GENE_ANALYSIS } from '../modules/hox_analysis'
include { LUMINAL_BASAL_ANALYSIS } from '../modules/luminal_basal_analysis'
include { UROTHELIAL_DIFFERENTIATION_ANALYSIS } from '../modules/urothelial_differentiation'
include { SHH_SIGNALLING_ANALYSIS } from '../modules/shh_signalling'
include { FGF_SIGNALLING_ANALYSIS } from '../modules/fgf_signalling'
include { WNT_BETACATENIN_SIGNALLING_ANALYSIS } from '../modules/wnt_betacatenin_signalling'
include { TGFBETA_SIGNALLING_ANALYSIS } from '../modules/tgfbeta_signalling'
include { INTERFERON_ANALYSIS } from '../modules/interferon_analysis'
include { IMMUNE_MARKERS_ANALYSIS } from '../modules/immune_markers'
include { PROLIFERATION_ANALYSIS } from '../modules/proliferation_analysis'

workflow ANALYSE_BIOLOGICAL_PATHWAYS_SHORT {
    take:
    counts         // path: featureCounts output
    sample_info    // path: sample metadata CSV
    analysis_group // val: column in sample_info to filter by

    main:
    log.info "=== Starting Biological Pathway Analysis (Short Reads) ==="

    // Step 1: TPM normalisation
    log.info "Step 1: Normalising counts to TPM"
    TPM_NORMALISATION(counts, "short_read")

    // Step 2: HOX gene expression
    log.info "Step 2: Analysing HOX gene expression"
    HOX_GENE_ANALYSIS(TPM_NORMALISATION.out.tpm_counts, sample_info, analysis_group)

    // Step 3: Luminal/basal phenotype
    log.info "Step 3: Analysing luminal/basal phenotype"
    LUMINAL_BASAL_ANALYSIS(TPM_NORMALISATION.out.tpm_counts, sample_info, analysis_group)

    // Step 4: Urothelial differentiation
    log.info "Step 4: Analysing urothelial differentiation"
    UROTHELIAL_DIFFERENTIATION_ANALYSIS(TPM_NORMALISATION.out.tpm_counts, sample_info, analysis_group)

    // Step 5: SHH signalling
    log.info "Step 5: Analysing SHH signalling pathway"
    SHH_SIGNALLING_ANALYSIS(TPM_NORMALISATION.out.tpm_counts, sample_info, analysis_group)

    // Step 6: FGF signalling
    log.info "Step 6: Analysing FGF signalling pathway"
    FGF_SIGNALLING_ANALYSIS(TPM_NORMALISATION.out.tpm_counts, sample_info, analysis_group)

    // Step 7: WNT signalling
    log.info "Step 7: Analysing WNT beta-catenin signalling pathway"
    WNT_BETACATENIN_SIGNALLING_ANALYSIS(TPM_NORMALISATION.out.tpm_counts, sample_info, analysis_group)

    // Step 8: TGF-beta signalling
    log.info "Step 8: Analysing TGF-beta signalling pathway"
    TGFBETA_SIGNALLING_ANALYSIS(TPM_NORMALISATION.out.tpm_counts, sample_info, analysis_group)

    // Step 9: Interferon response
    log.info "Step 9: Analysing interferon response"
    INTERFERON_ANALYSIS(TPM_NORMALISATION.out.tpm_counts, sample_info, analysis_group)

    // Step 10: Immune markers
    log.info "Step 10: Analysing immune cell markers"
    IMMUNE_MARKERS_ANALYSIS(TPM_NORMALISATION.out.tpm_counts, sample_info, analysis_group)

    // Step 11: Proliferation
    log.info "Step 11: Analysing proliferation and cell cycle"
    PROLIFERATION_ANALYSIS(TPM_NORMALISATION.out.tpm_counts, sample_info, analysis_group)

    emit:
    tpm_counts = TPM_NORMALISATION.out.tpm_counts
    
    // Summaries for reporting
    tpm_summary = TPM_NORMALISATION.out.summary
    hox_plots = HOX_GENE_ANALYSIS.out.plots_pdf
    luminal_plots = LUMINAL_BASAL_ANALYSIS.out.plots_pdf
    urothelial_plots = UROTHELIAL_DIFFERENTIATION_ANALYSIS.out.plots_pdf
    shh_plots = SHH_SIGNALLING_ANALYSIS.out.plots_pdf
    fgf_plots = FGF_SIGNALLING_ANALYSIS.out.plots_pdf
    wnt_plots = WNT_BETACATENIN_SIGNALLING_ANALYSIS.out.plots_pdf
    tgfbeta_plots = TGFBETA_SIGNALLING_ANALYSIS.out.plots_pdf
    interferon_plots = INTERFERON_ANALYSIS.out.plots_pdf
    immune_plots = IMMUNE_MARKERS_ANALYSIS.out.plots_pdf
    proliferation_plots = PROLIFERATION_ANALYSIS.out.plots_pdf
}