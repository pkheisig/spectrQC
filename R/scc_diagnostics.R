#' Plot SCC Diagnostics
#' 
#' Generates RRMSE distribution plots for each single-color control.
#' This helps the user decide on an appropriate RRMSE cutoff for refinement.
#' 
#' @param M Initial reference matrix
#' @param input_folder SCC folder
#' @param output_folder Folder to save diagnostic plots
#' @param custom_fluorophores Custom mapping
#' @export
plot_scc_diagnostics <- function(M, 
                                input_folder = "scc", 
                                output_folder = "spectrQC_outputs/plots/scc_diagnostics",
                                control_file = "fcs_control_file.csv",
                                custom_fluorophores = NULL) {
    library(flowCore)
    library(ggplot2)
    library(data.table)
    
    if (!is.null(output_folder)) {
        dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)
    }
    
    fcs_files <- list.files(input_folder, pattern = "\\.fcs$", full.names = TRUE)
    sample_patterns <- get_fluorophore_patterns()
    
    # Load control file if available for precise mapping
    file_map <- list()
    if (file.exists(control_file)) {
        tryCatch({
            cdf <- data.table::fread(control_file)
            if ("fluorophore" %in% colnames(cdf) && "filename" %in% colnames(cdf)) {
                # Map Fluorophore -> Filename
                # Handle duplicates by taking the first one
                valid_cdf <- cdf[fluorophore != "" & !is.na(fluorophore)]
                file_map <- split(valid_cdf$filename, valid_cdf$fluorophore)
            }
        }, error = function(e) message("Warning: Could not read control file: ", e$message))
    }

    get_name <- function(sn) {
        if (!is.null(custom_fluorophores) && sn %in% names(custom_fluorophores)) {
            return(custom_fluorophores[[sn]])
        }
        for (type in names(sample_patterns)) {
            patterns <- sample_patterns[[type]]
            patterns <- patterns[order(-nchar(patterns))]
            for (p in patterns) {
                if (grepl(p, sn, ignore.case = TRUE)) return(p)
            }
        }
        return(sn)
    }
    
    all_rrmse <- list()
    detector_names <- colnames(M)
    marker_names <- rownames(M)
    
    for (marker_name in marker_names) {
        # Find file matching marker_name
        match_f <- NULL
        
        # 1. Try Control File Map
        if (marker_name %in% names(file_map)) {
            target_fns <- file_map[[marker_name]]
            # Look for these filenames in fcs_files
            for (fn in target_fns) {
                hits <- fcs_files[basename(fcs_files) == fn]
                if (length(hits) > 0) {
                    match_f <- hits[1]
                    break
                }
            }
        }
        
        # 2. Fallback to heuristic
        if (is.null(match_f)) {
            for (f in fcs_files) {
                if (get_name(tools::file_path_sans_ext(basename(f))) == marker_name) {
                    match_f <- f
                    break
                }
            }
        }
        
        if (is.null(match_f)) {
            message("  Warning: No SCC file found for marker ", marker_name)
            next
        }
        
        message("  SCC Diagnostic: ", marker_name)
        ff <- read.FCS(match_f, transformation = FALSE, truncate_max_range = FALSE)
        raw_data <- exprs(ff)
        if (nrow(raw_data) > 5000) raw_data <- raw_data[sample(nrow(raw_data), 5000), ]
        
        sig <- M[marker_name, detector_names, drop = FALSE]
        sig_norm <- sig / max(sig)
        
        Y <- raw_data[, detector_names, drop = FALSE]
        A <- (Y %*% t(sig_norm)) / as.numeric(sig_norm %*% t(sig_norm))
        
        fitted <- A %*% sig_norm
        R <- Y - fitted
        rmse <- sqrt(rowMeans(R^2))
        rrmse <- (rmse / pmax(rowSums(Y), 1)) * 100 # Percentage
        
        all_rrmse[[marker_name]] <- data.table(Marker = marker_name, RRMSE = rrmse)
    }
    
    res <- rbindlist(all_rrmse)
    
    p <- ggplot(res, aes(x = RRMSE, fill = Marker)) +
        geom_density(alpha = 0.5) +
        scale_y_log10() + # Log scale to see tiny outlier peaks
        geom_vline(xintercept = 5, linetype = "dashed", color = "red") +
        facet_wrap(~Marker, scales = "free_y") +
        labs(title = "SCC Spectral Consistency (Initial Pass)",
             subtitle = "Dashed red line = 5% cutoff. Y-axis is LOG SCALE to highlight rare outliers.",
             x = "RRMSE (%)", y = "Density (Log)") +
        theme_minimal() +
        theme(legend.position = "none")
    
    if (!is.null(output_folder)) {
        ggsave(file.path(output_folder, "scc_rrmse_dist.png"), p, width = 250, height = 180, units = "mm")
    }
    
    # Also save a summary table
    summary <- res[, .(
        Median_RRMSE = median(RRMSE),
        Pct_Above_5 = sum(RRMSE > 5) / .N * 100
    ), by = Marker]
    
    if (!is.null(output_folder)) {
        fwrite(summary, file.path(output_folder, "scc_qc_summary.csv"))
    }
    
    return(p)
}
