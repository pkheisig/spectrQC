' Create AutoSpectral Control File
#' 
#' Generates the CSV file required by the AutoSpectral package.
#' 
#' @param input_folder Directory containing single-stained control FCS files.
#' @param output_file Path where the CSV will be saved (default: "fcs_control_file.csv").
#' @param custom_fluorophores Optional named vector to map filenames to fluorophore names.
#' @return A data frame containing the control file information.
#' @export
create_autospectral_control_file <- function(input_folder = "scc", 
                                          af_folder = "af",
                                          output_file = "fcs_control_file.csv",
                                          custom_fluorophores = NULL) {
    scc_files <- list.files(input_folder, pattern = "\\.fcs$", full.names = FALSE)
    af_files <- if(dir.exists(af_folder)) list.files(af_folder, pattern = "\\.fcs$", full.names = FALSE) else character(0)
    
    if (length(scc_files) == 0) stop("No FCS files found in ", input_folder)
    
    # 1. Identify patterns
    patterns_list <- get_fluorophore_patterns()
    specific_dyes <- c(patterns_list$beads, patterns_list$cells)
    generic_terms <- c("LIVE/DEAD", "Fixable Viability", "Viability", "FVD")
    specific_dyes <- setdiff(specific_dyes, generic_terms)
    specific_dyes <- specific_dyes[order(-nchar(specific_dyes))]
    fluor_patterns <- c(specific_dyes, generic_terms)

    # 2. Build entries for SCC files
    scc_rows <- list()
    for (fn in scc_files) {
        fluor <- "Unknown"
        if (!is.null(custom_fluorophores) && fn %in% names(custom_fluorophores)) {
            fluor <- custom_fluorophores[fn]
        } else {
            for (p in fluor_patterns) {
                if (grepl(p, fn, ignore.case = TRUE)) { fluor <- p; break }
            }
        }
        
        # Check for internal AF
        if (grepl("Unstained|US_UT", fn, ignore.case = TRUE)) fluor <- "AF_Internal"
        
        scc_rows[[fn]] <- data.frame(
            filename = fn,
            fluorophore = fluor,
            marker = "",
            channel = "",
            control.type = "cells",
            universal.negative = "",
            large.gate = "",
            is.viability = "",
            stringsAsFactors = FALSE
        )
    }
    
    # 3. Build entries for AF folder files
    af_rows <- list()
    for (i in seq_along(af_files)) {
        fn <- af_files[i]
        tag <- if(i == 1) "AF" else paste0("AF_", i)
        af_rows[[fn]] <- data.frame(
            filename = fn,
            fluorophore = tag,
            marker = "Autofluorescence",
            channel = "",
            control.type = "cells",
            universal.negative = fn, # Self-referencing for negative control
            large.gate = "TRUE",
            is.viability = "",
            stringsAsFactors = FALSE
        )
    }
    
    # 4. Combine
    # If we have files in 'af' folder, they take precedence over SCC-based detections
    df <- do.call(rbind, c(af_rows, scc_rows))
    
    # Strict uniqueness check: one row per filename
    # We keep the first occurrence (which would be the AF-folder one if duplicated)
    df <- df[!duplicated(df$filename), ]
    
    # Final Fixes
    primary_af_file <- if(length(af_files) > 0) af_files[1] else scc_files[grep("Unstained|US_UT", scc_files, ignore.case = TRUE)][1]
    if (is.na(primary_af_file)) stop("No Unstained/AF file found.")
    
    # Ensure all markers point to the primary AF file as their universal negative
    df$universal.negative <- primary_af_file
    
    # EXCEPTION: ALL AF-tagged files should have an empty string for universal.negative
    # This prevents AutoSpectral from trying to gate a negative for a negative
    df$universal.negative[grepl("^AF", df$fluorophore)] <- ""
    
    # Ensure the primary AF file itself is tagged as 'AF' in the fluorophore column
    df$fluorophore[df$filename == primary_af_file] <- "AF"

    # Identify channels
    for (i in seq_len(nrow(df))) {
        fn <- df$filename[i]
        path <- if(file.exists(file.path(af_folder, fn))) file.path(af_folder, fn) else file.path(input_folder, fn)
        tryCatch({
            ff <- flowCore::read.FCS(path, transformation = FALSE, truncate_max_range = FALSE)
            pd <- flowCore::pData(flowCore::parameters(ff))
            fl_pd <- get_sorted_detectors(pd)
            if (length(fl_pd$names) > 0) {
                medians <- apply(flowCore::exprs(ff)[, fl_pd$names, drop = FALSE], 2, median)
                df$channel[i] <- names(medians)[which.max(medians)]
            }
        }, error = function(e) NULL)
    }

    # USE STANDARD write.csv for compatibility with AutoSpectral's reader
    write.csv(df, output_file, row.names = FALSE, quote = TRUE)
    return(df)
}

#' Unmix using AutoSpectral
#' 
#' A wrapper function to use AutoSpectral's unmixing methods.
#' 
#' @param flow_frame A flowFrame object.
#' @param control_file Path to the AutoSpectral control file.
#' @param control_dir Directory containing control FCS files.
#' @param method Unmixing method ("WLS", "OLS", or "Poisson").
#' @param cytometer Cytometer type (e.g., "Aurora", "ID7000").
#' @export
#' Get Spectra via AutoSpectral (Robust Multi-AF)
#' 
#' Extracts SCC signatures using internal logic and AF signatures from the af/ folder,
#' then combines them for unmixing using the AutoSpectral WLS engine.
#' @export
get_autospectral_spectra <- function(flow_frame, 
                                   control_file = "fcs_control_file.csv", 
                                   control_dir = "scc",
                                   af_dir = "af",
                                   method = "WLS",
                                   cytometer = "Aurora") {
    if (!requireNamespace("AutoSpectral", quietly = TRUE)) {
        stop("Package 'AutoSpectral' required. Install from GitHub.")
    }
    
    # 1. Get detector info
    pd <- flowCore::pData(flowCore::parameters(flow_frame))
    det_info <- get_sorted_detectors(pd)
    detector_names <- det_info$names

    # 2. Extract SCC signatures
    message("  - Extracting reference signatures from single-color controls...")
    M_scc <- build_reference_matrix(input_folder = control_dir, control_df = data.table::fread(control_file))
    M_scc <- M_scc[, detector_names, drop = FALSE]

    # 3. Extract Multi-AF signatures
    if (dir.exists(af_dir)) {
        message("  - Extracting additional AF signatures from '", af_dir, "' folder...")
        M_af <- extract_af_signatures(af_dir, detector_names)
        M_expanded <- rbind(M_scc, M_af)
    } else {
        M_expanded <- M_scc
    }
    
    return(M_expanded[, detector_names, drop = FALSE])
}

#' Internal helper to extract gated AF signatures
extract_af_signatures <- function(af_dir, detector_names) {
    af_files <- list.files(af_dir, pattern = "\\.fcs$", full.names = TRUE)
    if (length(af_files) == 0) return(matrix(0, nrow = 0, ncol = length(detector_names)))
    
    spectra <- list()
    for (f in af_files) {
        ff <- flowCore::read.FCS(f, transformation = FALSE, truncate_max_range = FALSE)
        raw_data <- flowCore::exprs(ff)
        pd <- flowCore::pData(flowCore::parameters(ff))
        fsc <- pd$name[grepl("^FSC", pd$name) & grepl("-A$", pd$name)][1]
        ssc <- pd$name[grepl("^SSC", pd$name) & grepl("-A$", pd$name)][1]
        
        # 20-80% quantile gate on FSC/SSC to isolate main population
        idx <- which(raw_data[, fsc] > quantile(raw_data[, fsc], 0.2) & 
                     raw_data[, fsc] < quantile(raw_data[, fsc], 0.8) &
                     raw_data[, ssc] > quantile(raw_data[, ssc], 0.2) &
                     raw_data[, ssc] < quantile(raw_data[, ssc], 0.8))
        
        if (length(idx) < 10) idx <- seq_len(nrow(raw_data)) 
        
        sig <- apply(raw_data[idx, detector_names, drop = FALSE], 2, median)
        sig_norm <- sig / max(sig)
        
        tag <- if(length(spectra) == 0) "AF" else paste0("AF_", length(spectra) + 1)
        spectra[[tag]] <- sig_norm
    }
    return(do.call(rbind, spectra))
}

