#' Generate a Full QC PDF Report
#'
#' Creates a multi-page report summarizing unmixing quality, including spectra,
#' detector residuals, spread matrix, RRMSE scatter plots, NPS, and correlation panels.
#'
#' @param results_df Combined unmixed data frame (typically `rbind` of sample results).
#' @param M Reference matrix used for unmixing.
#' @param output_file Output PDF file path.
#' @param res_list Optional residual object/list from `calc_residuals(..., return_residuals = TRUE)`.
#' @param png_dir Deprecated and ignored (kept for backward compatibility).
#' @param pd Optional detector metadata (`flowCore::pData(parameters(ff))`) for axis labels.
#'
#' @return Invisibly returns `NULL`; writes report to disk.
#' @export
#' @examples
#' \dontrun{
#' # results_df <- data.table::rbindlist(lapply(unmixed_list, `[[`, "data"))
#' generate_qc_report(
#'   results_df = results_df,
#'   M = M,
#'   output_file = file.path("spectrQC_outputs", "Sample_QC_Report.pdf")
#' )
#' }
generate_qc_report <- function(results_df, M, output_file = file.path("spectrQC_outputs", "Sample_QC_Report.pdf"), res_list = NULL, png_dir = NULL, pd = NULL) {
    message("Generating spectrQC Summary Report...")
    if (!is.null(png_dir)) {
        warning("png_dir is deprecated and ignored; report output is PDF-only.")
    }
    results_df <- as.data.frame(results_df, stringsAsFactors = FALSE)
    if (!("File" %in% colnames(results_df))) {
        stop("results_df must contain a 'File' column.")
    }
    out_dir <- dirname(output_file)
    if (!is.na(out_dir) && nzchar(out_dir) && out_dir != ".") {
        dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    }
    grDevices::pdf(output_file, width = 11, height = 8.5)
    on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
    
    grid::grid.newpage()
    rrmse_values <- suppressWarnings(as.numeric(results_df$Relative_RMSE))
    rrmse_values <- rrmse_values[is.finite(rrmse_values)]
    median_rrmse_txt <- if (length(rrmse_values) > 0) {
        paste0(round(stats::median(rrmse_values) * 100, 2), "%")
    } else {
        "NA (not available)"
    }
    mean_rrmse_txt <- if (length(rrmse_values) > 0) {
        paste0(round(mean(rrmse_values) * 100, 2), "%")
    } else {
        "NA (not available)"
    }
    summary_txt <- paste0(
        "spectrQC: Spectral Unmixing Quality Control Report\n",
        "Generated on: ", Sys.time(), "\n\n",
        "Files Processed: ", length(unique(results_df$File)), "\n",
        "Total Markers: ", nrow(M), "\n",
        "Global Median RRMSE: ", median_rrmse_txt, "\n",
        "Global Mean RRMSE: ", mean_rrmse_txt
    )
    grid::grid.text(summary_txt, x = 0.5, y = 0.6, just = "center", gp = grid::gpar(fontsize = 15))
    
    message("  - Adding spectra overlay...")
    p_spectra <- plot_spectra(M, pd = pd, output_file = NULL)
    if (!is.null(p_spectra)) print(p_spectra)
    
    if (!is.null(res_list)) {
        message("  - Adding detector-level residual diagnostics...")
        rep_res <- if (!is.null(res_list$residuals)) res_list else res_list[[1]]
        p_det <- plot_detector_residuals(rep_res, M = M, top_n = 50, output_file = NULL)
        if (!is.null(p_det)) print(p_det)
    }
    
    message("  - Adding Spread Matrix...")
    ssm <- calculate_ssm(M)
    p_ssm <- plot_ssm(ssm, output_file = NULL)
    if (!is.null(p_ssm)) print(p_ssm)
    
    message("  - Adding RRMSE Scatter plots...")
    all_files <- unique(results_df$File)
    files_per_page <- 9 
    n_pages_scatter <- ceiling(length(all_files) / files_per_page)
    
    for (p_idx in seq_len(n_pages_scatter)) {
        start_f <- (p_idx - 1) * files_per_page + 1
        end_f <- min(p_idx * files_per_page, length(all_files))
        subset_df <- results_df[results_df$File %in% all_files[start_f:end_f], , drop = FALSE]
        p_scatter <- plot_scatter_rmse(subset_df, metric = "Relative_RMSE", output_file = NULL, color_limits = c(0, 5))
        if (!is.null(p_scatter)) {
            grid::grid.newpage()
            grid::pushViewport(grid::viewport(width = grid::unit(180, "mm"), height = grid::unit(180, "mm")))
            grid::grid.draw(ggplot2::ggplotGrob(p_scatter))
            grid::popViewport()
        }
    }
    
    message("  - Adding NPS diagnostics...")
    nps_scores <- calculate_nps(results_df)
    p_nps <- plot_nps(nps_scores, output_file = NULL)
    if (!is.null(p_nps)) print(p_nps)
    
    message("  - Adding marker-RRMSE correlation fits...")
    exclude_cols <- c("RMSE_Score", "Relative_RMSE", "File", "FSC-A", "SSC-A", "FSC-H", "SSC-H")
    all_markers <- setdiff(colnames(results_df), exclude_cols)
    markers_per_page <- 12 
    n_pages_corr <- ceiling(length(all_markers) / markers_per_page)
    
    for (p_idx in seq_len(n_pages_corr)) {
        start_m <- (p_idx - 1) * markers_per_page + 1
        end_m <- min(p_idx * markers_per_page, length(all_markers))
        subset_markers <- all_markers[start_m:end_m]
        p_corr <- tryCatch(
            plot_marker_correlations(
                results_df,
                markers = subset_markers,
                metric = "Relative_RMSE",
                output_file = NULL,
                y_limits = c(0, 10)
            ),
            error = function(e) {
                warning("Skipping marker correlation page ", p_idx, ": ", conditionMessage(e))
                NULL
            }
        )
        if (!is.null(p_corr)) {
            grid::grid.newpage()
            grid::pushViewport(grid::viewport(width = grid::unit(250, "mm"), height = grid::unit(180, "mm")))
            grid::grid.draw(ggplot2::ggplotGrob(p_corr))
            grid::popViewport()
        }
    }
    
    message("  - Adding Recommendations page...")
    grid::grid.newpage()
    high_spread <- which(ssm > 10, arr.ind = TRUE)
    spread_msgs <- if(nrow(high_spread) > 0) sapply(1:nrow(high_spread), function(i) paste0("- ", rownames(ssm)[high_spread[i,1]], " spreads noise heavily into ", colnames(ssm)[high_spread[i,2]])) else "- No extreme noise spread detected between markers."
    has_rrmse <- any(is.finite(suppressWarnings(as.numeric(results_df$Relative_RMSE))))
    if (has_rrmse) {
        bad_files <- results_df |>
            dplyr::group_by(File) |>
            dplyr::summarise(Mean_RRMSE = mean(Relative_RMSE, na.rm = TRUE), .groups = "drop") |>
            dplyr::filter(is.finite(Mean_RRMSE), Mean_RRMSE > 0.03) |>
            dplyr::pull(File)
        rrmse_msgs <- if (length(bad_files) > 0) {
            paste0("- High overall error in: ", paste(bad_files, collapse = ", "))
        } else {
            "- Overall spectral fit error is low (<3%) for all samples."
        }
    } else {
        rrmse_msgs <- "- RRMSE metrics are unavailable in results_df (all NA). Re-run unmix_samples() with a compatible reference matrix for QC."
    }
    wrap_lines <- function(x, width = 80) paste(strwrap(x, width = width), collapse = "\n")
    rec_txt <- paste0("spectrQC: Conclusions & Recommendations\n\n", "Spectral Spread Analysis:\n", wrap_lines(paste(spread_msgs, collapse = "\n")), "\n\n", "Spectral Fit Analysis:\n", wrap_lines(rrmse_msgs), "\n\n", "General Recommendations:\n", "1. If RRMSE > 5% points cluster in specific FSC/SSC regions, they may be debris.\n", "2. Markers with high spread should not be used to resolve dim co-expressed populations.\n", "3. If detector residuals show laser-specific patterns, check instrument calibration.")
    grid::grid.text(rec_txt, x = 0.1, y = 0.9, just = c("left", "top"), gp = grid::gpar(fontsize = 11, lineheight = 1.2))
    
    grDevices::dev.off()
    message("Report saved to: ", output_file)
}
