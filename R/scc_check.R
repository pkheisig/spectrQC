#' Generate SCC QC Report
#'
#' @param M Reference matrix.
#' @param scc_dir Directory containing SCC FCS files.
#' @param output_file Path to save the PDF report.
#' @param custom_fluorophores Named vector mapping filenames to fluorophores.
#' @param png_dir Directory to save individual PNG plots.
#' @param control_file Path to manual control mapping CSV.
#' @return Invisibly returns `NULL`; writes SCC audit PDF and PNG panels.
#' @examples
#' \dontrun{
#' generate_scc_report(
#'   M = M,
#'   scc_dir = "scc",
#'   output_file = "SCC_QC_Report.pdf",
#'   control_file = "fcs_control_file.csv"
#' )
#' }
#' @export
generate_scc_report <- function(M, scc_dir = "scc", output_file = "SCC_QC_Report.pdf", custom_fluorophores = NULL, png_dir = "spectrQC_outputs/plots/scc_audit", control_file = "fcs_control_file.csv") {
    message("Generating SCC-only QC Report...")
    out_dir <- dirname(output_file)
    if (!is.na(out_dir) && nzchar(out_dir) && out_dir != ".") {
        dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    }
    dir.create(png_dir, showWarnings = FALSE, recursive = TRUE)

    # 1. Setup metadata from first file
    fcs_files <- list.files(scc_dir, pattern = "\\.fcs$", full.names = TRUE)
    ff_meta <- flowCore::read.FCS(fcs_files[1], transformation = FALSE, truncate_max_range = FALSE)
    pd <- flowCore::pData(flowCore::parameters(ff_meta))

    # 2. Run residuals on SCCs
    res_list <- lapply(fcs_files, function(f) {
        sn <- tools::file_path_sans_ext(basename(f))
        ff <- flowCore::read.FCS(f, transformation = FALSE, truncate_max_range = FALSE)
        if (nrow(ff) > 5000) ff <- ff[sample(nrow(ff), 5000), ]
        calc_residuals(ff, M, file_name = sn, method = "WLS", return_residuals = TRUE)
    })

    all_data <- data.table::rbindlist(lapply(res_list, `[[`, "data"))

    grDevices::pdf(output_file, width = 11, height = 8.5)

    # Page 1: Title
    grid::grid.newpage()
    grid::grid.text("spectrQC: Single Color Control Audit", x = 0.5, y = 0.6, gp = grid::gpar(fontsize = 20))

    # Page 2: Reference Spectra
    message("  - Adding spectra...")
    print(plot_spectra(M, pd = pd, output_file = file.path(png_dir, "01_scc_spectra.png")))

    # Page 3: Detector Residuals (REORDERED)
    message("  - Adding residuals...")
    print(plot_detector_residuals(res_list[[1]], M = M, top_n = 50, output_file = file.path(png_dir, "02_scc_detector_residuals.png"), pd = pd))

    # Page 4: SSM
    message("  - Adding SSM...")
    ssm <- calculate_ssm(M)
    print(plot_ssm(ssm, output_file = file.path(png_dir, "03_scc_ssm.png")))

    # Page 5: SCC Diagnostics (Consistency)
    message("  - Adding diagnostics...")
    print(plot_scc_diagnostics(M, input_folder = scc_dir, custom_fluorophores = custom_fluorophores, output_folder = png_dir, control_file = control_file))

    grDevices::dev.off()
    message("SCC Report saved to: ", output_file)
}
