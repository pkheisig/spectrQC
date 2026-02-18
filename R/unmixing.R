#' Unmix Experimental Samples
#' 
#' @param sample_dir Directory containing experimental FCS files.
#' @param M Reference matrix (Markers x Detectors).
#' @param method Unmixing method ("WLS", "OLS", "NNLS", or "AutoSpectral").
#' @param output_dir Directory to save unmixed FCS files.
#' @return A named list with one element per sample. Each element contains
#'   `data` (unmixed abundances + QC metrics) and `residuals` (detector residual matrix).
#' @export
unmix_samples <- function(sample_dir = "samples", 
                          M = NULL, 
                          method = "WLS", 
                          output_dir = "samples_unmixed") {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    fcs_files <- list.files(sample_dir, pattern = "\\.fcs$", full.names = TRUE)
    
    if (length(fcs_files) == 0) stop("No FCS files found in ", sample_dir)
    if (is.null(M)) stop("Reference Matrix M must be provided for unmixing.")

    method_upper <- toupper(method)
    allowed_methods <- c("WLS", "OLS", "NNLS", "AUTOSPECTRAL")
    if (!(method_upper %in% allowed_methods)) {
        stop("method must be one of: ", paste(allowed_methods, collapse = ", "))
    }

    build_result_from_unmixed <- function(flow_frame, M_sub, abundances, file_name) {
        full_data <- flowCore::exprs(flow_frame)
        detectors <- colnames(M_sub)
        Y <- full_data[, detectors, drop = FALSE]
        fitted <- abundances %*% M_sub
        residuals <- Y - fitted

        rmse <- sqrt(rowMeans(residuals^2))
        relative_rmse <- rmse / pmax(rowSums(Y), 1)

        out <- as.data.frame(abundances)
        colnames(out) <- rownames(M_sub)
        out$RMSE_Score <- rmse
        out$Relative_RMSE <- relative_rmse

        all_cols <- colnames(full_data)
        fsc_col <- grep("^FSC[0-9]*-A$", all_cols, value = TRUE)[1]
        ssc_col <- grep("^SSC[0-9]*-A$", all_cols, value = TRUE)[1]
        if (!is.na(fsc_col)) out[["FSC-A"]] <- full_data[, fsc_col]
        if (!is.na(ssc_col)) out[["SSC-A"]] <- full_data[, ssc_col]
        out$File <- file_name

        list(data = out, residuals = residuals)
    }

    results <- list()
    for (f in fcs_files) {
        sn <- tools::file_path_sans_ext(basename(f))
        message("  Unmixing sample: ", sn)
        ff <- flowCore::read.FCS(f, transformation = FALSE, truncate_max_range = FALSE)
        
        # Select unmixing math engine
        if (method_upper == "AUTOSPECTRAL") {
            if (!requireNamespace("AutoSpectral", quietly = TRUE)) {
                warning("AutoSpectral package not found. Falling back to internal WLS.")
                res_obj <- calc_residuals(ff, M, method = "WLS", file_name = sn, return_residuals = TRUE)
            } else {
                # Use AutoSpectral WLS math engine directly with our refined matrix M
                raw_data <- flowCore::exprs(ff)
                detectors <- colnames(M)
                missing <- setdiff(detectors, colnames(raw_data))
                if (length(missing) > 0) {
                    stop("Detectors in reference matrix not found in sample '", sn, "': ", paste(missing, collapse = ", "))
                }
                Y <- raw_data[, detectors, drop = FALSE]
                signatures <- M[, detectors, drop = FALSE]
                
                # AutoSpectral unmix.wls
                unmixed_data <- AutoSpectral::unmix.wls(Y, signatures)
                abundances <- as.matrix(unmixed_data)
                colnames(abundances) <- rownames(signatures)
                res_obj <- build_result_from_unmixed(ff, signatures, abundances, sn)
            }
        } else {
            res_obj <- calc_residuals(ff, M, method = method_upper, file_name = sn, return_residuals = TRUE)
        }
        
        # Save results as FCS instead of CSV
        markers_to_keep <- intersect(colnames(res_obj$data), rownames(M))
        scatter_cols <- grep("FSC|SSC", colnames(res_obj$data), value = TRUE)
        cols_to_write <- c(markers_to_keep, scatter_cols)
        unmixed_exprs <- as.matrix(res_obj$data[, cols_to_write, drop = FALSE])
        
        new_ff <- flowCore::flowFrame(unmixed_exprs)
        # Copy keywords from original if possible (optional but good practice)
        new_ff@description <- ff@description
        
        flowCore::write.FCS(new_ff, file.path(output_dir, paste0(sn, "_unmixed.fcs")))
        results[[sn]] <- res_obj
    }
    
    return(results)
}
