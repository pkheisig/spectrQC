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
unmix_samples <- function(sample_dir = "samples", 
                          M = NULL, 
                          method = "WLS", 
                          output_dir = "samples_unmixed",
                          control_file = "fcs_control_file.csv",
                          scc_dir = "scc",
                          af_dir = "af") {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    fcs_files <- list.files(sample_dir, pattern = "\\.fcs$", full.names = TRUE)
    
    if (length(fcs_files) == 0) stop("No FCS files found in ", sample_dir)
    
    # If method is AutoSpectral, we need to ensure the control file is set up for multi-AF
    if (toupper(method) == "AUTOSPECTRAL") {
        message("  AutoSpectral mode: Updating control file and extracting signatures...")
        create_autospectral_control_file(input_folder = scc_dir, af_folder = af_dir, output_file = control_file)
        
        # We use a dummy flowFrame just to trigger the spectra extraction once
        dummy_f <- list.files(sample_dir, pattern = "\\.fcs$", full.names = TRUE)[1]
        dummy_ff <- flowCore::read.FCS(dummy_f, transformation = FALSE, truncate_max_range = FALSE)
        
        # Get the global AutoSpectral spectra
        asp_spectra <<- get_autospectral_spectra(dummy_ff, control_file = control_file, control_dir = scc_dir, af_dir = af_dir, method = method)
        M <- asp_spectra # Use these for residuals
    }

    results <- list()
    for (f in fcs_files) {
        sn <- tools::file_path_sans_ext(basename(f))
        message("  Unmixing sample: ", sn)
        ff <- flowCore::read.FCS(f, transformation = FALSE, truncate_max_range = FALSE)
        
        if (toupper(method) == "AUTOSPECTRAL") {
            # Use the pre-calculated spectra
            raw_data <- flowCore::exprs(ff)
            # Detectors are already matched in get_autospectral_spectra
            Y <- raw_data[, colnames(asp_spectra), drop = FALSE]
            unmixed_data <- AutoSpectral::unmix.wls(Y, asp_spectra)
            
            # Calculate residuals using the same spectra
            res <- calc_residuals(ff, asp_spectra, method = "WLS", file_name = sn)
        } else {
            res <- calc_residuals(ff, M, method = method, file_name = sn)
        }
        
        data.table::fwrite(res, file.path(output_dir, paste0(sn, "_unmixed.csv")))
        results[[sn]] <- list(data = res)
    }
    
    return(results)
}
