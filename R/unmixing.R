#' Unmix Experimental Samples
#' 
#' @param sample_dir Directory containing experimental FCS files.
#' @param M Optional reference matrix (Markers x Detectors). Required for
#'   dynamic unmixing methods (`"WLS"`, `"OLS"`, `"NNLS"`, `"AutoSpectral"`).
#' @param W Optional static unmixing matrix (Markers x Detectors). If supplied,
#'   unmixing is performed using matrix multiplication with the transposed W.
#' @param unmixing_matrix_file Optional CSV path to a saved unmixing matrix.
#'   Used when `W` is not supplied. By default this points to the matrix produced
#'   by [autounmix_controls()].
#' @param method Unmixing method ("WLS", "OLS", "NNLS", or "AutoSpectral").
#' @param cytometer Cytometer name used when `method = "AutoSpectral"` (for example `"Aurora"`).
#' @param output_dir Directory to save unmixed FCS files.
#' @param write_fcs Logical; if `TRUE`, write unmixed FCS files to `output_dir`.
#' @return A named list with one element per sample. Each element contains
#'   `data` (unmixed abundances + QC metrics when available) and `residuals`
#'   (detector residual matrix, `NULL` if no reference matrix is available).
#' @examples
#' \dontrun{
#' unmixed <- unmix_samples(
#'   sample_dir = "samples",
#'   unmixing_matrix_file = "spectrQC_outputs/autounmix_controls/scc_unmixing_matrix.csv",
#'   output_dir = "spectrQC_outputs/unmix_samples"
#' )
#' names(unmixed)
#' }
#' @export
unmix_samples <- function(sample_dir = "samples", 
                          M = NULL, 
                          W = NULL,
                          unmixing_matrix_file = file.path("spectrQC_outputs", "autounmix_controls", "scc_unmixing_matrix.csv"),
                          method = "WLS", 
                          cytometer = "Aurora",
                          output_dir = file.path("spectrQC_outputs", "unmix_samples"),
                          write_fcs = TRUE) {
    if (!is.null(M)) {
        M <- .as_reference_matrix(M, "M")
    }
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    fcs_files <- list.files(sample_dir, pattern = "\\.fcs$", full.names = TRUE)
    
    if (length(fcs_files) == 0) stop("No FCS files found in ", sample_dir)

    read_unmixing_matrix_csv <- function(path) {
        if (!file.exists(path)) stop("unmixing_matrix_file not found: ", path)
        df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
        if (ncol(df) < 2) stop("Invalid unmixing matrix CSV: expected at least 2 columns in ", path)

        first_name <- tolower(trimws(colnames(df)[1]))
        first_col <- df[[1]]
        has_marker_col <- first_name %in% c("marker", "fluorophore", "file") || !is.numeric(first_col)
        if (has_marker_col) {
            marker_names <- trimws(as.character(first_col))
            mat_df <- df[, -1, drop = FALSE]
        } else {
            marker_names <- rownames(df)
            mat_df <- df
        }
        W_mat <- as.matrix(mat_df)
        storage.mode(W_mat) <- "numeric"
        if (is.null(marker_names) || all(marker_names == "")) {
            marker_names <- paste0("marker_", seq_len(nrow(W_mat)))
        }
        rownames(W_mat) <- marker_names
        colnames(W_mat) <- colnames(mat_df)
        W_mat
    }

    using_static_W <- FALSE
    W_use <- NULL
    if (!is.null(W)) {
        W_use <- as.matrix(W)
        using_static_W <- TRUE
    } else if (is.null(M) && !is.null(unmixing_matrix_file)) {
        W_use <- read_unmixing_matrix_csv(unmixing_matrix_file)
        using_static_W <- TRUE
    }

    if (is.null(M) && !using_static_W) {
        stop(
            "No unmixing input provided. Supply either:\n",
            " - M (reference matrix), or\n",
            " - W / unmixing_matrix_file (static unmixing matrix)."
        )
    }
    if (!is.null(M) && using_static_W) {
        message("Both M and W/unmixing_matrix_file provided. Using static unmixing matrix.")
    }

    method_upper <- toupper(method)
    allowed_methods <- c("WLS", "OLS", "NNLS", "AUTOSPECTRAL")
    if (!using_static_W && !(method_upper %in% allowed_methods)) {
        stop("method must be one of: ", paste(allowed_methods, collapse = ", "))
    }
    if (!using_static_W && method_upper == "AUTOSPECTRAL" && requireNamespace("AutoSpectral", quietly = TRUE)) {
        cytometer_candidates <- unique(c(cytometer, tolower(cytometer), toupper(cytometer)))
        ok <- FALSE
        for (cand in cytometer_candidates) {
            ok <- tryCatch({
                AutoSpectral::get.autospectral.param(cytometer = cand, figures = FALSE)
                TRUE
            }, error = function(e) FALSE)
            if (ok) break
        }
        if (!ok) stop("Unsupported cytometer for AutoSpectral: ", cytometer)
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

    build_result_from_static_unmix <- function(flow_frame, W_sub, file_name, M_for_qc = NULL) {
        full_data <- flowCore::exprs(flow_frame)
        detectors <- colnames(W_sub)
        Y <- full_data[, detectors, drop = FALSE]
        abundances <- Y %*% t(W_sub)
        colnames(abundances) <- rownames(W_sub)

        rmse <- rep(NA_real_, nrow(Y))
        relative_rmse <- rep(NA_real_, nrow(Y))
        residuals <- NULL

        if (!is.null(M_for_qc)) {
            common_markers <- intersect(rownames(W_sub), rownames(M_for_qc))
            common_detectors <- intersect(detectors, colnames(M_for_qc))
            if (length(common_markers) > 0 && length(common_detectors) > 0) {
                A_qc <- abundances[, common_markers, drop = FALSE]
                M_qc <- M_for_qc[common_markers, common_detectors, drop = FALSE]
                Y_qc <- Y[, common_detectors, drop = FALSE]
                fitted <- A_qc %*% M_qc
                residuals <- Y_qc - fitted
                rmse <- sqrt(rowMeans(residuals^2))
                relative_rmse <- rmse / pmax(rowSums(Y_qc), 1)
            }
        }

        out <- as.data.frame(abundances)
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

        if (using_static_W) {
            raw_data <- flowCore::exprs(ff)
            detectors <- colnames(W_use)
            missing <- setdiff(detectors, colnames(raw_data))
            if (length(missing) > 0) {
                stop("Detectors in unmixing matrix not found in sample '", sn, "': ", paste(missing, collapse = ", "))
            }
            res_obj <- build_result_from_static_unmix(ff, W_use, sn, M_for_qc = M)
        } else {
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
        }
        
        if (isTRUE(write_fcs)) {
            marker_source <- if (using_static_W) rownames(W_use) else rownames(M)
            markers_to_keep <- intersect(colnames(res_obj$data), marker_source)
            scatter_cols <- grep("FSC|SSC", colnames(res_obj$data), value = TRUE)
            cols_to_write <- c(markers_to_keep, scatter_cols)
            unmixed_exprs <- as.matrix(res_obj$data[, cols_to_write, drop = FALSE])
            
            new_ff <- flowCore::flowFrame(unmixed_exprs)
            # Do not copy raw metadata wholesale; old spillover/parameter keywords
            # can become inconsistent with the unmixed channels and break downstream tools.
            
            flowCore::write.FCS(new_ff, file.path(output_dir, paste0(sn, "_unmixed.fcs")))
        }
        results[[sn]] <- res_obj
    }
    
    return(results)
}
