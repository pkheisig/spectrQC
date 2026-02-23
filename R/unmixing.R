#' Unmix Experimental Samples
#' 
#' @param sample_dir Directory containing experimental FCS files.
#' @param M Optional reference matrix (Markers x Detectors). Required for
#'   dynamic unmixing methods (`"WLS"`, `"OLS"`, `"NNLS"`, `"AutoSpectral"`).
#' @param W Optional static unmixing matrix (Markers x Detectors). If supplied,
#'   unmixing is performed using matrix multiplication with the transposed W.
#' @param unmixing_matrix_file Optional CSV path to a saved unmixing matrix.
#'   Used when `W` is not supplied. By default this points to the matrix produced
#'   by [autounmix_controls()]. In static mode, if `M` is not supplied,
#'   spectrQC attempts to auto-load a paired reference matrix from the same
#'   directory (for QC metrics).
#' @param method Unmixing method ("WLS", "OLS", "NNLS", or "AutoSpectral").
#' @param cytometer Cytometer name used when `method = "AutoSpectral"` (for example `"Aurora"`).
#' @param output_dir Directory to save unmixed FCS files.
#' @param write_fcs Logical; if `TRUE`, write unmixed FCS files to `output_dir`.
#' @return A named list with one element per sample. Each element contains
#'   `data` (unmixed abundances + QC metrics when available) and `residuals`
#'   (detector residual matrix, `NULL` if no compatible reference matrix is available).
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
    normalize_key <- function(x) {
        out <- toupper(trimws(as.character(x)))
        out <- gsub("[^A-Z0-9]", "", out)
        out[is.na(out)] <- ""
        out
    }

    map_names <- function(source_names, target_names) {
        source_keys <- normalize_key(source_names)
        target_keys <- normalize_key(target_names)
        idx <- match(source_keys, target_keys)
        mapped <- target_names[idx]
        mapped[is.na(idx)] <- NA_character_
        mapped
    }

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

    read_reference_matrix_csv <- function(path) {
        if (!file.exists(path)) stop("reference matrix file not found: ", path)
        df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
        .as_reference_matrix(df, paste0("reference matrix file ('", path, "')"))
    }

    using_static_W <- FALSE
    W_use <- NULL
    auto_loaded_qc_matrix <- FALSE
    if (!is.null(W)) {
        W_use <- as.matrix(W)
        using_static_W <- TRUE
    } else if (is.null(M) && !is.null(unmixing_matrix_file)) {
        W_use <- read_unmixing_matrix_csv(unmixing_matrix_file)
        using_static_W <- TRUE
    }

    # For static unmixing from CSV, auto-load the paired reference matrix when available
    # so QC metrics (RMSE/RRMSE) can be computed without extra user inputs.
    if (is.null(M) && is.null(W) && using_static_W && !is.null(unmixing_matrix_file)) {
        matrix_dir <- dirname(unmixing_matrix_file)
        matrix_file <- basename(unmixing_matrix_file)
        ref_candidates <- unique(c(
            file.path(matrix_dir, "scc_reference_matrix.csv"),
            file.path(matrix_dir, sub("unmixing_matrix", "reference_matrix", matrix_file)),
            file.path(matrix_dir, sub("unmixing", "reference", matrix_file))
        ))
        ref_candidates <- ref_candidates[file.exists(ref_candidates)]
        if (length(ref_candidates) > 0) {
            M <- read_reference_matrix_csv(ref_candidates[[1]])
            auto_loaded_qc_matrix <- TRUE
            message("Auto-loaded reference matrix for QC metrics: ", ref_candidates[[1]])
        } else {
            warning(
                "Static unmixing is using only W/unmixing_matrix_file. ",
                "No paired reference matrix was found for QC metrics, so RMSE/RRMSE will be NA."
            )
        }
    }

    if (is.null(M) && !using_static_W) {
        stop(
            "No unmixing input provided. Supply either:\n",
            " - M (reference matrix), or\n",
            " - W / unmixing_matrix_file (static unmixing matrix)."
        )
    }
    if (!auto_loaded_qc_matrix && !is.null(M) && using_static_W) {
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
        fsc_a_col <- grep("^FSC[0-9]*-A$", all_cols, value = TRUE)[1]
        ssc_a_col <- grep("^SSC[0-9]*-A$", all_cols, value = TRUE)[1]
        fsc_h_col <- grep("^FSC[0-9]*-H$", all_cols, value = TRUE)[1]
        ssc_h_col <- grep("^SSC[0-9]*-H$", all_cols, value = TRUE)[1]
        if (!is.na(fsc_a_col)) out[["FSC-A"]] <- full_data[, fsc_a_col]
        if (!is.na(ssc_a_col)) out[["SSC-A"]] <- full_data[, ssc_a_col]
        if (!is.na(fsc_h_col)) out[["FSC-H"]] <- full_data[, fsc_h_col]
        if (!is.na(ssc_h_col)) out[["SSC-H"]] <- full_data[, ssc_h_col]
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
            marker_map <- data.frame(
                w_name = rownames(W_sub),
                m_name = map_names(rownames(W_sub), rownames(M_for_qc)),
                stringsAsFactors = FALSE
            )
            marker_map <- marker_map[!is.na(marker_map$m_name) & marker_map$m_name != "", , drop = FALSE]
            marker_map <- marker_map[!duplicated(marker_map$m_name), , drop = FALSE]

            detector_map <- data.frame(
                w_name = detectors,
                m_name = map_names(detectors, colnames(M_for_qc)),
                stringsAsFactors = FALSE
            )
            detector_map <- detector_map[!is.na(detector_map$m_name) & detector_map$m_name != "", , drop = FALSE]
            detector_map <- detector_map[!duplicated(detector_map$m_name), , drop = FALSE]

            if (nrow(marker_map) > 0 && nrow(detector_map) > 0) {
                A_qc <- abundances[, marker_map$w_name, drop = FALSE]
                colnames(A_qc) <- marker_map$m_name
                M_qc <- M_for_qc[marker_map$m_name, detector_map$m_name, drop = FALSE]
                Y_qc <- Y[, detector_map$w_name, drop = FALSE]
                colnames(Y_qc) <- detector_map$m_name
                fitted <- A_qc %*% M_qc
                residuals <- Y_qc - fitted
                rmse <- sqrt(rowMeans(residuals^2, na.rm = TRUE))
                rmse[rowSums(is.finite(residuals)) == 0] <- NA_real_
                total_intensity <- rowSums(Y_qc, na.rm = TRUE)
                relative_rmse <- rmse / pmax(total_intensity, 1)
                relative_rmse[!is.finite(relative_rmse)] <- NA_real_
            } else {
                warning(
                    "QC metrics unavailable for sample '", file_name,
                    "': static unmixing matrix and reference matrix have no compatible marker/detector names."
                )
            }
        } else {
            warning(
                "QC metrics unavailable for sample '", file_name,
                "': no reference matrix (M) supplied for static unmixing."
            )
        }

        out <- as.data.frame(abundances)
        out$RMSE_Score <- rmse
        out$Relative_RMSE <- relative_rmse
        all_cols <- colnames(full_data)
        fsc_a_col <- grep("^FSC[0-9]*-A$", all_cols, value = TRUE)[1]
        ssc_a_col <- grep("^SSC[0-9]*-A$", all_cols, value = TRUE)[1]
        fsc_h_col <- grep("^FSC[0-9]*-H$", all_cols, value = TRUE)[1]
        ssc_h_col <- grep("^SSC[0-9]*-H$", all_cols, value = TRUE)[1]
        if (!is.na(fsc_a_col)) out[["FSC-A"]] <- full_data[, fsc_a_col]
        if (!is.na(ssc_a_col)) out[["SSC-A"]] <- full_data[, ssc_a_col]
        if (!is.na(fsc_h_col)) out[["FSC-H"]] <- full_data[, fsc_h_col]
        if (!is.na(ssc_h_col)) out[["SSC-H"]] <- full_data[, ssc_h_col]
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
