#' Generate Full QC Report
#' 
#' Bundles all specRQC diagnostic plots into a single PDF report.
#' 
#' @param results_df Unmixed data frame from calc_residuals
#' @param M Reference matrix
#' @param output_file Path to save the PDF report
#' @param res_list Optional. List from calc_residuals with return_residuals=TRUE
#' @export
generate_qc_report <- function(results_df, M, output_file = "specRQC_Report.pdf", res_list = NULL) {
    message("Generating specRQC Summary Report...")
    
    # Check dependencies
    if (!requireNamespace("gridExtra", quietly = TRUE)) {
        stop("Package 'gridExtra' required for report generation. install.packages('gridExtra')")
    }
    
    pdf(output_file, width = 11, height = 8.5) # Landscape
    
    # Page 1: Title & Summary Stats
    grid::grid.newpage()
    summary_txt <- paste0(
        "specRQC: Spectral Unmixing Quality Control Report\n",
        "Generated on: ", Sys.time(), "\n\n",
        "Files Processed: ", length(unique(results_df$File)), "\n",
        "Total Markers: ", nrow(M), "\n",
        "Global Median RRMSE: ", round(median(results_df$Relative_RMSE) * 100, 2), "\n",
        "Global Mean RRMSE: ", round(mean(results_df$Relative_RMSE) * 100, 2), "%"
    )
    grid::grid.text(summary_txt, x = 0.5, y = 0.6, just = "center", gp = grid::gpar(fontsize = 15))
    
    # Page 2: Spectral Reference Matrix
    message("  - Adding spectra overlay...")
    p_spectra <- plot_spectra(M, output_file = NULL)
    if (!is.null(p_spectra)) print(p_spectra)
    
    # Page 3: Spectral Spread Matrix (SSM)
    message("  - Adding Spread Matrix...")
    ssm <- calculate_ssm(M)
    p_ssm <- plot_ssm(ssm, output_file = NULL)
    if (!is.null(p_ssm)) print(p_ssm)
    
    # Page 4: RRMSE Scatter Plot
    message("  - Adding RRMSE Scatter plots...")
    p_scatter <- plot_scatter_rmse(results_df, metric = "Relative_RMSE", output_file = NULL, color_limits = c(0, 5))
    if (!is.null(p_scatter)) print(p_scatter)
    
    # Page 5: Negative Population Spread (NPS)
    message("  - Adding NPS diagnostics...")
    nps_scores <- calculate_nps(results_df)
    p_nps <- plot_nps(nps_scores, output_file = NULL)
    if (!is.null(p_nps)) print(p_nps)
    
    # Page 6: Detector Residuals (if provided)
    if (!is.null(res_list)) {
        message("  - Adding detector-level residual diagnostics...")
        # res_list might be a list of lists if multiple files processed.
        # We'll just show the first one or combine them?
        # Let's show the first one as representative.
        if (!is.null(res_list$residuals)) {
            p_det <- plot_detector_residuals(res_list, top_n = 50, output_file = NULL)
            if (!is.null(p_det)) print(p_det)
        } else if (is.list(res_list) && length(res_list) > 0) {
            p_det <- plot_detector_residuals(res_list[[1]], top_n = 50, output_file = NULL)
            if (!is.null(p_det)) print(p_det)
        }
    }
    
    # Page 7+: Marker Correlations (one page each)
    message("  - Adding marker-RRMSE correlation fits...")
    p_corr <- plot_marker_correlations(results_df, metric = "Relative_RMSE", output_file = NULL, y_limits = c(0, 10))
    if (!is.null(p_corr)) print(p_corr)
    
    dev.off()
    message("Report saved to: ", output_file)
}
