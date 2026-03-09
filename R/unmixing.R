#' Unmix Experimental Samples
#' 
#' @param sample_dir Directory containing experimental FCS files.
#' @param M Optional reference matrix (Markers x Detectors). Required for
#'   dynamic unmixing methods (`"WLS"`, `"OLS"`, `"NNLS"`).
#' @param W Optional static unmixing matrix (Markers x Detectors). If supplied,
#'   unmixing is performed using matrix multiplication with the transposed W.
#' @param unmixing_matrix_file Optional CSV path to a saved unmixing matrix.
#'   Used when `W` is not supplied. By default this points to the matrix produced
#'   by [autounmix_controls()].
#' @param method Unmixing method ("WLS", "OLS", or "NNLS").
#' @param cytometer Reserved for compatibility with older workflows.
#' @param output_dir Directory to save unmixed FCS files.
#' @param write_fcs Logical; if `TRUE`, write unmixed FCS files to `output_dir`.
#' @return A named list with one element per sample. Each element contains
#'   `data` (unmixed abundances plus retained acquisition parameters) and
#'   `residuals` (detector residual matrix when available, otherwise `NULL`).
#' @examples
#' \dontrun{
#' unmixed <- unmix_samples(
#'   sample_dir = "samples",
#'   unmixing_matrix_file = "spectreasy_outputs/autounmix_controls/scc_unmixing_matrix.csv",
#'   output_dir = "spectreasy_outputs/unmix_samples"
#' )
#' names(unmixed)
#' }
#' @export
unmix_samples <- function(sample_dir = "samples", 
                          M = NULL, 
                          W = NULL,
                          unmixing_matrix_file = file.path("spectreasy_outputs", "autounmix_controls", "scc_unmixing_matrix.csv"),
                          method = "WLS", 
                          cytometer = "Aurora",
                          output_dir = file.path("spectreasy_outputs", "unmix_samples"),
                          write_fcs = TRUE) {
    if (!is.null(M)) {
        M <- .as_reference_matrix(M, "M")
    }
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    fcs_files <- list.files(sample_dir, pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE)
    
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

    resolve_secondary_label_map <- function(primary_names) {
        primary_names <- trimws(as.character(primary_names))
        labels <- stats::setNames(primary_names, primary_names)
        if (length(primary_names) == 0) return(labels)

        opt_control <- getOption("spectreasy.control_file", "")
        sample_parent <- dirname(normalizePath(sample_dir, mustWork = FALSE))
        candidates <- unique(c(
            as.character(opt_control),
            .resolve_control_file_path("fcs_mapping.csv"),
            file.path(sample_parent, "fcs_mapping.csv")
        ))
        candidates <- candidates[!is.na(candidates) & nzchar(trimws(candidates))]
        existing <- candidates[file.exists(candidates)]
        if (length(existing) == 0) return(labels)

        control_df <- tryCatch(
            utils::read.csv(existing[1], stringsAsFactors = FALSE, check.names = FALSE),
            error = function(e) NULL
        )
        if (is.null(control_df) || nrow(control_df) == 0) return(labels)

        lower_names <- tolower(colnames(control_df))
        fluor_idx <- match("fluorophore", lower_names)
        marker_idx <- match("marker", lower_names)
        if (is.na(fluor_idx) || is.na(marker_idx)) return(labels)

        fluor_vals <- trimws(as.character(control_df[[fluor_idx]]))
        marker_vals <- trimws(as.character(control_df[[marker_idx]]))
        keep <- !is.na(fluor_vals) & fluor_vals != "" & !is.na(marker_vals) & marker_vals != ""
        if (!any(keep)) return(labels)

        key <- tolower(fluor_vals[keep])
        val <- marker_vals[keep]
        map_ci <- stats::setNames(val, key)
        map_ci <- map_ci[!duplicated(names(map_ci))]

        for (nm in names(labels)) {
            k <- tolower(trimws(nm))
            if (k %in% names(map_ci) && nzchar(map_ci[[k]])) {
                labels[[nm]] <- map_ci[[k]]
            }
        }
        labels
    }

    apply_feature_secondary_labels <- function(target_ff, source_ff, marker_cols, secondary_label_map) {
        pd_new <- flowCore::pData(flowCore::parameters(target_ff))
        if (!all(c("name", "desc") %in% colnames(pd_new))) {
            return(target_ff)
        }

        param_names <- as.character(pd_new$name)
        desc_vals <- param_names

        # Preserve source descriptions for passthrough acquisition parameters when available.
        pd_src <- tryCatch(flowCore::pData(flowCore::parameters(source_ff)), error = function(e) NULL)
        if (!is.null(pd_src) && all(c("name", "desc") %in% colnames(pd_src))) {
            src_name <- as.character(pd_src$name)
            src_desc <- as.character(pd_src$desc)
            src_desc[is.na(src_desc)] <- ""
            src_map <- stats::setNames(src_desc, src_name)
            matched <- param_names %in% names(src_map)
            desc_vals[matched] <- ifelse(
                nzchar(src_map[param_names[matched]]),
                src_map[param_names[matched]],
                param_names[matched]
            )
        }

        marker_cols <- as.character(marker_cols)
        for (mk in marker_cols) {
            idx <- which(param_names == mk)
            if (length(idx) == 0) next
            lbl <- if (mk %in% names(secondary_label_map)) secondary_label_map[[mk]] else mk
            if (is.na(lbl) || !nzchar(trimws(lbl))) lbl <- mk
            desc_vals[idx] <- as.character(lbl)
        }

        pd_new$desc <- desc_vals
        flowCore::parameters(target_ff) <- methods::new("AnnotatedDataFrame", data = pd_new)
        target_ff
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
    allowed_methods <- c("WLS", "OLS", "NNLS")
    if (!using_static_W && !(method_upper %in% allowed_methods)) {
        stop("method must be one of: ", paste(allowed_methods, collapse = ", "))
    }

    build_result_from_unmixed <- function(flow_frame, M_sub, abundances, file_name) {
        full_data <- flowCore::exprs(flow_frame)
        detectors <- colnames(M_sub)
        Y <- full_data[, detectors, drop = FALSE]
        fitted <- abundances %*% M_sub
        residuals <- Y - fitted

        out <- as.data.frame(abundances)
        colnames(out) <- rownames(M_sub)

        out <- .append_passthrough_parameters(out, full_data, detector_names = detectors)
        out$File <- file_name

        list(data = out, residuals = residuals)
    }

    build_result_from_static_unmix <- function(flow_frame, W_sub, file_name) {
        full_data <- flowCore::exprs(flow_frame)
        detectors <- colnames(W_sub)
        Y <- full_data[, detectors, drop = FALSE]
        abundances <- Y %*% t(W_sub)
        colnames(abundances) <- rownames(W_sub)
        residuals <- NULL

        out <- as.data.frame(abundances)
        out <- .append_passthrough_parameters(out, full_data, detector_names = detectors)
        out$File <- file_name

        list(data = out, residuals = residuals)
    }

    results <- list()
    marker_source_all <- if (using_static_W) rownames(W_use) else rownames(M)
    secondary_label_map <- resolve_secondary_label_map(marker_source_all)

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
            res_obj <- build_result_from_static_unmix(ff, W_use, sn)
        } else {
            res_obj <- calc_residuals(ff, M, method = method_upper, file_name = sn, return_residuals = TRUE)
        }
        
        if (isTRUE(write_fcs)) {
            marker_source <- if (using_static_W) rownames(W_use) else rownames(M)
            markers_to_keep <- intersect(colnames(res_obj$data), marker_source)
            passthrough_cols <- .get_passthrough_parameter_names(colnames(res_obj$data))
            cols_to_write <- unique(c(markers_to_keep, passthrough_cols))
            unmixed_exprs <- as.matrix(res_obj$data[, cols_to_write, drop = FALSE])
            
            new_ff <- flowCore::flowFrame(unmixed_exprs)
            new_ff <- apply_feature_secondary_labels(
                target_ff = new_ff,
                source_ff = ff,
                marker_cols = markers_to_keep,
                secondary_label_map = secondary_label_map
            )
            # Do not copy raw metadata wholesale; old spillover/parameter keywords
            # can become inconsistent with the unmixed channels and break downstream tools.
            
            flowCore::write.FCS(new_ff, file.path(output_dir, paste0(sn, "_unmixed.fcs")))
        }
        results[[sn]] <- res_obj
    }
    
    return(results)
}
