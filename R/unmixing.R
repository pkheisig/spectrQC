#' Unmix Experimental Samples
#' 
#' Performs unmixing on all samples in a folder using the refined matrix.
#' 
#' @param sample_dir Folder with experimental FCS files.
#' @param M Refined reference matrix.
#' @param method Unmixing method ("OLS", "WLS", or "NNLS").
#' @param output_dir Folder to save unmixed CSV files.
#' @return A list of unmixed data frames.
#' @export
#' Unmix Experimental Samples
#' 
#' @param sample_dir Directory containing experimental FCS files.
#' @param M Reference matrix (Markers x Detectors).
#' @param method Unmixing method ("WLS", "OLS", "NNLS", or "AutoSpectral").
#' @param output_dir Directory to save unmixed CSV files.
#' @return A list of unmixed results.
#' @export
unmix_samples <- function(sample_dir = "samples", 
                          M = NULL, 
                          method = "WLS", 
                          output_dir = "samples_unmixed") {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    fcs_files <- list.files(sample_dir, pattern = "\\.fcs$", full.names = TRUE)
    
    if (length(fcs_files) == 0) stop("No FCS files found in ", sample_dir)
    if (is.null(M)) stop("Reference Matrix M must be provided for unmixing.")

    results <- list()
    for (f in fcs_files) {
        sn <- tools::file_path_sans_ext(basename(f))
        message("  Unmixing sample: ", sn)
        ff <- flowCore::read.FCS(f, transformation = FALSE, truncate_max_range = FALSE)
        
        # Select unmixing math engine
        if (toupper(method) == "AUTOSPECTRAL") {
            if (!requireNamespace("AutoSpectral", quietly = TRUE)) {
                warning("AutoSpectral package not found. Falling back to internal WLS.")
                res <- calc_residuals(ff, M, method = "WLS", file_name = sn)
            } else {
                # Use AutoSpectral WLS math engine directly with our refined matrix M
                raw_data <- flowCore::exprs(ff)
                detectors <- intersect(colnames(raw_data), colnames(M))
                Y <- raw_data[, detectors, drop = FALSE]
                signatures <- M[, detectors, drop = FALSE]
                
                # AutoSpectral unmix.wls
                unmixed_data <- AutoSpectral::unmix.wls(Y, signatures)
                
                # Still use calc_residuals to generate the full audit structure (RRMSE, etc)
                # as AutoSpectral unmix.wls returns a pure matrix.
                res <- calc_residuals(ff, M, method = "WLS", file_name = sn)
            }
        } else {
            res <- calc_residuals(ff, M, method = method, file_name = sn)
        }
        
        # Save results as FCS instead of CSV
        markers_to_keep <- intersect(colnames(res), rownames(M))
        scatter_cols <- grep("FSC|SSC", colnames(res), value = TRUE)
        unmixed_exprs <- as.matrix(res[, c(markers_to_keep, scatter_cols)])
        
        new_ff <- flowCore::flowFrame(unmixed_exprs)
        # Copy keywords from original if possible (optional but good practice)
        new_ff@description <- ff@description
        
        flowCore::write.FCS(new_ff, file.path(output_dir, paste0(sn, "_unmixed.fcs")))
        results[[sn]] <- list(data = res)
    }
    
    return(results)
}
