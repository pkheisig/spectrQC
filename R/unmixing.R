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
                         M, 
                         method = "WLS", 
                         output_dir = "unmixed_samples") {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    files <- list.files(sample_dir, pattern = "\\.fcs$", full.names = TRUE)
    if (length(files) == 0) stop("No samples found in ", sample_dir)
    
    results <- list()
    for (f in files) {
        name <- tools::file_path_sans_ext(basename(f))
        message("  Unmixing sample: ", name)
        ff <- flowCore::read.FCS(f)
        
        # Perform unmixing (returning residuals for QC step)
        res_list <- calc_residuals(ff, M, file_name = name, method = method, return_residuals = TRUE)
        
        # Save unmixed CSV
        data.table::fwrite(res_list$data, file.path(output_dir, paste0(name, "_unmixed.csv")))
        
        results[[name]] <- res_list
    }
    
    return(results)
}
