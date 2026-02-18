#' Quick Unmix Workflow: Spectra, SCC Unmixing, and Unmixing Matrix Plot
#'
#' Convenience helper for a minimal SCC workflow:
#' 1) build reference matrix,
#' 2) plot spectra,
#' 3) unmix SCC files,
#' 4) derive and plot unmixing matrix.
#'
#' @param scc_dir Directory with SCC FCS files.
#' @param control_df Optional control mapping data.frame.
#' @param output_dir Output directory.
#' @param unmix_method Unmixing method for SCC files ("WLS", "OLS", "NNLS", "AutoSpectral").
#' @param build_qc_plots Logical. If TRUE, keep `build_reference_matrix()` gating/spectrum plot exports.
#' @param ... Additional parameters passed to `build_reference_matrix()`.
#' @return List with `M`, `W`, `unmixed_list`, `spectra_plot`, and `unmixing_plot`.
#' @export
quick_unmix <- function(
    scc_dir = "scc",
    control_df = NULL,
    output_dir = "spectrQC_outputs/quick_unmix",
    unmix_method = "WLS",
    build_qc_plots = FALSE,
    ...
) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    if (!dir.exists(scc_dir)) stop("scc_dir not found: ", scc_dir)

    build_plots_dir <- file.path(output_dir, "build_reference_plots")
    unmixed_dir <- file.path(output_dir, "scc_unmixed")
    spectra_file <- file.path(output_dir, "scc_spectra.png")
    unmixing_matrix_csv <- file.path(output_dir, "scc_unmixing_matrix.csv")
    unmixing_matrix_png <- file.path(output_dir, "scc_unmixing_matrix.png")

    M <- build_reference_matrix(
        input_folder = scc_dir,
        output_folder = build_plots_dir,
        save_qc_plots = build_qc_plots,
        control_df = control_df,
        ...
    )
    if (is.null(M) || nrow(M) == 0) stop("No valid spectra found while building reference matrix.")

    fcs_files <- list.files(scc_dir, pattern = "\\.fcs$", full.names = TRUE)
    ff_meta <- flowCore::read.FCS(fcs_files[1], transformation = FALSE, truncate_max_range = FALSE)
    pd <- flowCore::pData(flowCore::parameters(ff_meta))

    p_spectra <- plot_spectra(M, pd = pd, output_file = spectra_file)

    unmixed_list <- unmix_samples(
        sample_dir = scc_dir,
        M = M,
        method = unmix_method,
        output_dir = unmixed_dir
    )

    W <- derive_unmixing_matrix(M, method = "OLS")
    save_unmixing_matrix(W, unmixing_matrix_csv)
    p_unmix <- plot_unmixing_matrix(W, pd = pd)
    ggplot2::ggsave(unmixing_matrix_png, p_unmix, width = 200, height = 150, units = "mm")

    invisible(list(
        M = M,
        W = W,
        unmixed_list = unmixed_list,
        spectra_plot = p_spectra,
        unmixing_plot = p_unmix
    ))
}
