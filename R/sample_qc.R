#' Generate Sample QC Plots
#' 
#' Produces RRMSE scatter plots and individual file diagnostics for unmixed experimental samples.
#' 
#' @param unmixed_list List of unmixed data (from unmix_samples).
#' @param M Reference matrix used.
#' @param output_dir Folder to save QC plots.
#' @param report_file Path to final PDF report.
#' @export
generate_sample_qc <- function(unmixed_list, 
                              M, 
                              output_dir = "samples_unmixed_qc",
                              report_file = "Experimental_Sample_Audit.pdf",
                              png_dir = "spectrQC_outputs/plots/sample_audit",
                              pd = NULL) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(png_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Combined data
    all_data <- data.table::rbindlist(lapply(unmixed_list, `[[`, "data"))
    
    # 1. Setup metadata for labels (from first sample if provided)
    
    pdf(report_file, width = 11, height = 8.5)
    
    grid::grid.newpage()
    grid::grid.text("spectrQC: Experimental Sample Audit", x = 0.5, y = 0.6, gp = grid::gpar(fontsize = 20))
    
    # Page: Unmixing Matrix
    message("  - Adding unmixing matrix...")
    W <- derive_unmixing_matrix(M, method = "OLS") 
    save_unmixing_matrix(W, "unmixing_matrix.csv") # Save as CSV as requested
    p_w <- plot_unmixing_matrix(W, pd = pd)
    ggplot2::ggsave(file.path(png_dir, "01_unmixing_matrix.png"), p_w, width = 200, height = 150, units = "mm")
    print(p_w)
    
    # Page: Global Scatters
    message("  - Adding scatter plots...")
    all_files <- unique(all_data$File)
    files_per_page <- 9
    for (p_idx in seq_len(ceiling(length(all_files) / files_per_page))) {
        start_f <- (p_idx - 1) * files_per_page + 1
        end_f <- min(p_idx * files_per_page, length(all_files))
        subset_df <- all_data[File %in% all_files[start_f:end_f]]
        p <- plot_scatter_rmse(subset_df, metric = "Relative_RMSE", output_file = file.path(png_dir, paste0("02_rrmse_scatter_pg", p_idx, ".png")))
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(width = grid::unit(180, "mm"), height = grid::unit(180, "mm")))
        grid::grid.draw(ggplot2::ggplotGrob(p))
        grid::popViewport()
    }
    
    # Page: Detector Residuals
    message("  - Adding residual deconstruction...")
    p_res <- plot_detector_residuals(unmixed_list[[1]], M = M, top_n = 50, output_file = file.path(png_dir, "03_sample_detector_residuals.png"), pd = pd)
    print(p_res)
    
    dev.off()
    message("Sample QC Report saved to: ", report_file)
}
