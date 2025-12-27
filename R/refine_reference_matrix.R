#' Refine Reference Matrix with RRMSE Cutoff
#' 
#' Cleans up single-color controls by removing events that don't fit the spectral model.
#' 
#' @param M Current reference matrix
#' @param input_folder Directory containing SCC files
#' @param rrmse_threshold Cutoff for Relative RMSE (default 0.05)
#' @param custom_fluorophores Named vector mapping filenames to fluorophore names
#' @return A refined reference matrix
#' @export
refine_reference_matrix <- function(M, 
                                   input_folder = "scc", 
                                   rrmse_threshold = 0.05,
                                   custom_fluorophores = NULL) {
    library(flowCore)
    library(data.table)
    
    fcs_files <- list.files(input_folder, pattern = "\\.fcs$", full.names = TRUE)
    if (length(fcs_files) == 0) stop("No SCC files found in ", input_folder)
    
    # We need to find which SCC matches which marker in M
    # Use the same naming logic
    sample_patterns <- get_fluorophore_patterns()
    
    get_name <- function(sn) {
        if (!is.null(custom_fluorophores) && sn %in% names(custom_fluorophores)) {
            return(custom_fluorophores[[sn]])
        }
        # Simplified match for this function
        for (type in names(sample_patterns)) {
            patterns <- sample_patterns[[type]]
            patterns <- patterns[order(-nchar(patterns))]
            for (p in patterns) {
                if (grepl(p, sn, ignore.case = TRUE)) return(p)
            }
        }
        return(sn)
    }
    
    refined_spectra <- list()
    detector_names <- colnames(M)
    
    for (f in fcs_files) {
        sn <- tools::file_path_sans_ext(basename(f))
        marker_name <- get_name(sn)
        
        if (!(marker_name %in% rownames(M))) {
            message("  Skipping ", sn, " (not in reference matrix)")
            next
        }
        
        message("Refining ", marker_name, " from ", sn, "...")
        ff <- read.FCS(f, transformation = FALSE, truncate_max_range = FALSE)
        
        # 1. Gate positive events (simplified for this refinement)
        # We reuse the peak-channel logic
        raw_data <- exprs(ff)
        # Subsample for speed if needed
        if (nrow(raw_data) > 10000) {
            raw_data <- raw_data[sample(nrow(raw_data), 10000), ]
        }
        
        # 2. Unmix against the specific signature in M
        # For an SCC, we only care about its OWN signature fit
        sig <- M[marker_name, detector_names, drop = FALSE]
        sig_norm <- sig / max(sig)
        
        Y <- raw_data[, detector_names, drop = FALSE]
        # Single-marker unmixing: abundance A = Y %*% sig^T / (sig %*% sig^T)
        A <- (Y %*% t(sig_norm)) / as.numeric(sig_norm %*% t(sig_norm))
        
        # 3. Calc residuals
        fitted <- A %*% sig_norm
        R <- Y - fitted
        rmse <- sqrt(rowMeans(R^2))
        rrmse <- rmse / pmax(rowSums(Y), 1)
        
        # 4. Filter
        keep <- rrmse <= rrmse_threshold
        n_removed <- sum(!keep)
        
        if (sum(keep) < 10) {
            warning("  Too few events passed cutoff for ", marker_name, ". Keeping all.")
            clean_data <- Y
        } else {
            message("    Removed ", n_removed, " outliers (", round(100*n_removed/length(keep), 1), "%)")
            clean_data <- Y[keep, ]
        }
        
        # 5. Extract median spectrum from cleaned data
        new_sig <- apply(clean_data, 2, median)
        refined_spectra[[marker_name]] <- new_sig / max(new_sig)
    }
    
    # Rebuild M
    M_refined <- do.call(rbind, refined_spectra)
    rownames(M_refined) <- names(refined_spectra)
    colnames(M_refined) <- detector_names
    
    # Ensure all original markers are present (if some skipped, copy from original M)
    missing <- setdiff(rownames(M), rownames(M_refined))
    if (length(missing) > 0) {
        message("Copying original signatures for missing markers: ", paste(missing, collapse = ", "))
        M_refined <- rbind(M_refined, M[missing, , drop = FALSE])
    }
    
    # Sort to match original M order
    M_refined <- M_refined[rownames(M), , drop = FALSE]
    
    return(M_refined)
}
