#' Refine SCC Matrix
#' 
#' Cleans SCCs using an RRMSE threshold, produces refined M and W, 
#' and generates a before/after comparison report.
#' 
#' @param M Initial reference matrix.
#' @param scc_dir SCC folder.
#' @param rrmse_threshold Cutoff for RRMSE cleaning.
#' @param output_dir Folder for refined matrices and unmixed SCC files.
#' @param custom_fluorophores Manual mapping.
#' @param report_file Path to comparison PDF report.
#' @return A list with refined [[M]] and [[W]].
#' @export
refine_scc_matrix <- function(M, 
                             scc_dir = "scc", 
                             rrmse_threshold = 0.05,
                             output_dir = "scc_unmixed",
                             control_file = "fcs_control_file.csv",
                             report_file = "SCC_Refinement_Report.pdf",
                             png_dir = "spectrQC_outputs/plots/scc_refinement") {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(png_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Load metadata for labels
    fcs_files <- list.files(scc_dir, pattern = "\\.fcs$", full.names = TRUE)
    if (length(fcs_files) == 0) stop("No SCC files found in ", scc_dir)
    ff_meta <- flowCore::read.FCS(fcs_files[1], transformation = FALSE, truncate_max_range = FALSE)
    pd <- flowCore::pData(flowCore::parameters(ff_meta))
    
    # Load control_df
    if (!file.exists(control_file)) stop("control_file not found: ", control_file)
    control_df <- utils::read.csv(control_file, stringsAsFactors = FALSE, check.names = FALSE)
    
    # 1. Refine M
    M_refined <- refine_reference_matrix(M, input_folder = scc_dir, rrmse_threshold = rrmse_threshold, control_df = control_df)
    
    # 2. Derive W (OLS by default for static matrix)
    W_refined <- derive_unmixing_matrix(M_refined, method = "OLS")
    
    # 3. Save matrices
    M_refined_df <- as.data.frame(M_refined, check.names = FALSE)
    M_refined_df$Marker <- rownames(M_refined)
    M_refined_df <- M_refined_df[, c("Marker", setdiff(colnames(M_refined_df), "Marker")), drop = FALSE]
    utils::write.csv(M_refined_df, file.path(output_dir, "refined_reference_matrix.csv"), row.names = FALSE, quote = TRUE)
    save_unmixing_matrix(W_refined, file.path(output_dir, "refined_unmixing_matrix.csv"))
    
    # 4. Save comparison plots
    message("  - Plotting initial vs refined spectra...")
    p1 <- plot_spectra(M, pd = pd, output_file = file.path(png_dir, "01_spectra_initial.png")) + ggplot2::ggtitle("Initial Spectra")
    p2 <- plot_spectra(M_refined, pd = pd, output_file = file.path(png_dir, "02_spectra_refined.png")) + ggplot2::ggtitle("Refined Spectra (Outliers Removed)")
    grDevices::pdf(report_file, width = 11, height = 8.5)
    grid::grid.newpage()
    grid::grid.text("spectrQC: SCC Refinement Comparison", x = 0.5, y = 0.6, gp = grid::gpar(fontsize = 20))
    print(p1)
    print(p2)
    grDevices::dev.off()
    message("  - SCC refinement report saved to: ", report_file)
    
    # 5. Unmix SCCs with refined W and save
    message("  - Saving unmixed SCC files (FCS format)...")
    
    for (f in fcs_files) {
        sn <- tools::file_path_sans_ext(basename(f))
        ff <- flowCore::read.FCS(f)
        res <- calc_residuals(ff, M_refined, method = "WLS")
        
        # Convert to matrix for FCS - standard subsetting for data.frame
        markers_to_keep <- intersect(colnames(res), rownames(M_refined))
        scatter_cols <- grep("FSC|SSC", colnames(res), value = TRUE)
        unmixed_exprs <- as.matrix(res[, c(markers_to_keep, scatter_cols)])
        
        new_ff <- flowCore::flowFrame(unmixed_exprs)
        new_ff@description <- ff@description
        flowCore::write.FCS(new_ff, file.path(output_dir, paste0(sn, "_unmixed.fcs")))
    }
    
    return(list(M = M_refined, W = W_refined))
}
