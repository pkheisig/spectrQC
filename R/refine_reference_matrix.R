#' Refine a Reference Matrix by Filtering High-Residual Events
#'
#' Reprocesses SCC files against an initial matrix `M`, removes high-RRMSE events,
#' and recomputes cleaner marker spectra.
#'
#' @param M Initial reference matrix (markers x detectors).
#' @param input_folder SCC directory with `.fcs` files.
#' @param rrmse_threshold Maximum RRMSE retained during event cleaning.
#' @param control_df Optional control mapping data.frame or CSV path.
#'
#' @return Refined reference matrix with same detector columns and marker order as `M`.
#' @export
#' @examples
#' \dontrun{
#' M_refined <- refine_reference_matrix(
#'   M = M,
#'   input_folder = "scc",
#'   rrmse_threshold = 0.05,
#'   control_df = "fcs_control_file.csv"
#' )
#' }
refine_reference_matrix <- function(M, 
                                   input_folder = "scc", 
                                   rrmse_threshold = 0.05,
                                   control_df = NULL) {
    if (is.character(control_df) && length(control_df) == 1 && !is.na(control_df)) {
        if (!file.exists(control_df)) stop("control_df file not found: ", control_df)
        control_df <- utils::read.csv(control_df, stringsAsFactors = FALSE, check.names = FALSE)
    } else if (!is.null(control_df) && !is.data.frame(control_df)) {
        stop("control_df must be either a data.frame or a single CSV path.")
    }
    if (is.data.frame(control_df)) {
        if (!("filename" %in% colnames(control_df))) {
            stop("control_df is missing required column: filename")
        }
        if (!("fluorophore" %in% colnames(control_df))) control_df$fluorophore <- ""
        if (!("channel" %in% colnames(control_df))) control_df$channel <- ""
        if (!("universal.negative" %in% colnames(control_df))) control_df$universal.negative <- ""

        control_df$filename <- trimws(as.character(control_df$filename))
        control_df$fluorophore <- trimws(as.character(control_df$fluorophore))
        control_df$channel <- trimws(as.character(control_df$channel))
        control_df$universal.negative <- trimws(as.character(control_df$universal.negative))
    }
    
    fcs_files <- list.files(input_folder, pattern = "\\.fcs$", full.names = TRUE)
    if (length(fcs_files) == 0) stop("No SCC files found in ", input_folder)
    
    get_control_rows <- function(df, filename) {
        if (is.null(df) || !("filename" %in% colnames(df))) {
            return(data.frame())
        }
        fn <- as.character(df$filename)
        df[fn == filename, ]
    }

    # We need to find which SCC matches which marker in M
    sample_patterns <- get_fluorophore_patterns()
    
    get_name <- function(sn_ext) {
        sn <- tools::file_path_sans_ext(sn_ext)
        if (!is.null(control_df)) {
            row_info <- get_control_rows(control_df, sn_ext)
            if (nrow(row_info) == 0) row_info <- get_control_rows(control_df, sn)
            if (nrow(row_info) > 0 && !is.na(row_info$fluorophore[1])) return(row_info$fluorophore[1])
        }
        # Fallback
        for (type in names(sample_patterns)) {
            patterns <- sample_patterns[[type]]
            patterns <- patterns[order(-nchar(patterns))]
            for (p in patterns) {
                if (grepl(p, sn, ignore.case = TRUE)) return(p)
            }
        }
        return(sn)
    }
    
    refined_spectra_raw <- list()
    detector_names <- colnames(M)
    
    # 1. First Pass: Find AF baseline
    af_data_raw <- NULL
    for (f in fcs_files) {
        sn_ext <- basename(f)
        if (get_name(sn_ext) == "AF") {
            message("  Establishing AF baseline for refinement from ", sn_ext)
            ff <- flowCore::read.FCS(f, transformation = FALSE, truncate_max_range = FALSE)
            raw_data <- flowCore::exprs(ff)
            if (nrow(raw_data) > 5000) raw_data <- raw_data[sample(nrow(raw_data), 5000), ]
            af_data_raw <- apply(raw_data[, detector_names, drop = FALSE], 2, median)
            break
        }
    }
    
    # 2. Second Pass: Refine all markers
    for (f in fcs_files) {
        sn_ext <- basename(f)
        marker_name <- get_name(sn_ext)
        
        if (!(marker_name %in% rownames(M))) {
            next
        }
        
        if (marker_name == "AF") next
        
        # Get info from control_df
        row_info <- get_control_rows(control_df, sn_ext)
        if (nrow(row_info) == 0) row_info <- get_control_rows(control_df, tools::file_path_sans_ext(sn_ext))
        
        message("Refining ", marker_name, " (", sn_ext, ")")
        ff <- flowCore::read.FCS(f, transformation = FALSE, truncate_max_range = FALSE)
        raw_data <- flowCore::exprs(ff)
        
        # 1. Unmix against the specific signature in M
        sig <- M[marker_name, detector_names, drop = FALSE]
        sig_norm <- sig / max(sig)
        
        Y <- raw_data[, detector_names, drop = FALSE]
        # Single-marker unmixing
        A <- (Y %*% t(sig_norm)) / as.numeric(sig_norm %*% t(sig_norm))
        
        # 2. Calc residuals and RRMSE
        fitted <- A %*% sig_norm
        R <- Y - fitted
        rmse <- sqrt(rowMeans(R^2))
        rrmse <- rmse / pmax(rowSums(Y), 1)
        
        # 3. Filter
        keep <- rrmse <= rrmse_threshold
        if (sum(keep) < 10) {
            warning("  Too few events passed cutoff for ", marker_name, ". Keeping all.")
            clean_Y <- Y
        } else {
            clean_Y <- Y[keep, ]
        }
        
        # 4. Extract signature with proper negative handling
        peak_ch <- if (nrow(row_info) > 0 && !is.na(row_info$channel[1]) && row_info$channel[1] != "") {
            row_info$channel[1]
        } else {
            detector_names[which.max(sig_norm)]
        }
        
        vals <- clean_Y[, peak_ch]
        
        neg_idx <- vals <= quantile(vals, 0.15)
        pos_idx <- vals >= quantile(vals, 0.95)
        if (sum(pos_idx) < 50) pos_idx <- vals >= quantile(vals, 0.90)
        
        pos_spectrum <- apply(clean_Y[pos_idx, , drop=FALSE], 2, median)
        neg_spectrum <- apply(clean_Y[neg_idx, , drop=FALSE], 2, median)
        
        # Universal negative logic
        val <- if (nrow(row_info) > 0 && "universal.negative" %in% colnames(row_info)) {
            trimws(as.character(row_info$universal.negative[1]))
        } else {
            ""
        }
        use_universal <- toupper(val) %in% c("TRUE", "AF")
        
        if (use_universal && !is.null(af_data_raw)) {
            sig_pure <- pmax(pos_spectrum - af_data_raw, 0)
        } else {
            sig_pure <- pmax(pos_spectrum - neg_spectrum, 0)
        }
        
        refined_spectra_raw[[marker_name]] <- if(max(sig_pure)>0) sig_pure else pos_spectrum
    }
    
    # Add AF
    if (!is.null(af_data_raw)) {
        refined_spectra_raw[["AF"]] <- af_data_raw
    }
    
    # Rebuild M with normalization
    M_refined_list <- lapply(refined_spectra_raw, function(s) {
        if (max(s) == 0) return(s)
        s / max(s)
    })
    M_refined <- do.call(rbind, M_refined_list)
    rownames(M_refined) <- names(M_refined_list)
    colnames(M_refined) <- detector_names
    
    # Ensure all original markers are present
    missing <- setdiff(rownames(M), rownames(M_refined))
    if (length(missing) > 0) {
        M_refined <- rbind(M_refined, M[missing, , drop = FALSE])
    }
    
    # Sort to match original M order
    M_refined <- M_refined[rownames(M), , drop = FALSE]
    
    return(M_refined)
}
