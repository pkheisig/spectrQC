#' Plot Detector-Level Residuals
#' 
#' Identifies which detectors contribute most to the unmixing error for high-RRMSE cells.
#' Overlays the reference signatures to help identify the source of the mismatch.
#' 
#' @param res_list List returned by calc_residuals with return_residuals=TRUE
#' @param M Reference matrix
#' @param top_n Number of top high-error cells to analyze
#' @param output_file Path to save the plot
#' @param pd Optional pData for descriptive labels
#' @export
plot_detector_residuals <- function(res_list, M, top_n = 50, output_file = "detector_residuals.png", width = 250, height = 120, pd = NULL) {
    data <- res_list$data
    residuals <- res_list$residuals
    
    if (is.null(residuals) || nrow(residuals) == 0) {
        warning("No residuals provided to plot_detector_residuals. Skipping plot.")
        return(NULL)
    }

    # Identify high-error cells
    idx <- order(data$Relative_RMSE, decreasing = TRUE)[1:min(top_n, nrow(data))]
    R_sub <- residuals[idx, , drop = FALSE]
    
    # Convert to long format for plotting
    long <- as.data.frame(R_sub)
    long$Cell <- seq_len(nrow(R_sub))
    long <- tidyr::pivot_longer(long, cols = -Cell, names_to = "Detector", values_to = "Residual")
    
    # Sort and label detectors
    det_names <- colnames(M)
    if (!is.null(pd)) {
        det_info <- get_sorted_detectors(pd)
        common <- intersect(det_info$names, det_names)
        levels_sorted <- common
        labels_sorted <- det_info$labels[match(common, det_info$names)]
    } else {
        nums <- as.numeric(gsub("[^0-9]", "", det_names))
        levels_sorted <- det_names[order(nums)]
        labels_sorted <- levels_sorted
    }
    
    long$Detector <- factor(long$Detector, levels = levels_sorted, labels = labels_sorted)
    
    # Prepare M for overlay
    # Scale M to the max absolute residual for visual comparison
    max_res <- max(abs(long$Residual), na.rm = TRUE)
    M_overlay <- as.data.frame(M[, levels_sorted, drop = FALSE])
    M_overlay$Fluorophore <- rownames(M)
    M_long <- tidyr::pivot_longer(M_overlay, cols = -Fluorophore, names_to = "Detector", values_to = "Signature")
    M_long$Signature <- M_long$Signature * max_res
    M_long$Detector <- factor(M_long$Detector, levels = levels_sorted, labels = labels_sorted)
    
    p <- ggplot2::ggplot() +
        # Spectra overlay (background)
        ggplot2::geom_line(data = M_long, ggplot2::aes(Detector, Signature, group = Fluorophore, color = Fluorophore), 
                           alpha = 0.5, linewidth = 0.5) +
        # Residual boxplots (foreground)
        ggplot2::geom_boxplot(data = long, ggplot2::aes(Detector, Residual), 
                              outlier.size = 0.5, fill = "steelblue", alpha = 0.6) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
        ggplot2::labs(title = paste("Residual Contributions for Top", top_n, "High-Error Cells"),
                      subtitle = "Background: Reference Spectra scaled to max residual. Foreground: Residual distribution.",
                      x = "Detector", y = "Residual Value / Scaled Signature") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 6))
    
    if (!is.null(output_file)) {
        ggplot2::ggsave(output_file, p, width = width, height = height, units = "mm", dpi = 300)
    }
    return(p)
}

#' Calculate Negative Population Spread (NPS)
#' 
#' Quantifies the spreading of negative populations after unmixing. 
#' High spread indicates poor unmixing or high spillover noise.
#' 
#' @param data Unmixed data frame
#' @param markers Vector of marker names to analyze
#' @return A data frame with NPS values (MAD) per marker per file
#' @export
calculate_nps <- function(data, markers = NULL) {
    if (is.null(markers)) {
        exclude <- c("RMSE_Score", "Relative_RMSE", "FSC-A", "SSC-A", "File")
        markers <- setdiff(colnames(data), exclude)
    }
    
    # For each marker, isolate the negative population (intensity < threshold)
    # We use a robust estimate (Median Absolute Deviation) of the spread
    nps_results <- data |> 
        dplyr::group_by(File) |> 
        dplyr::summarize(dplyr::across(dplyr::all_of(markers), function(x) {
            # Heuristic: negative population is around 0. 
            # We take values between -2SD and +2SD to avoid true positives
            # But simpler: just take the MAD of all values < quantile(0.2)
            neg_subset <- x[x < quantile(x, 0.2)]
            stats::mad(neg_subset, na.rm = TRUE)
        })) |> 
        tidyr::pivot_longer(cols = -File, names_to = "Marker", values_to = "NPS")
    
    return(nps_results)
}

#' Plot Negative Population Spread
#' @export
plot_nps <- function(nps_results, output_file = "nps_plot.png", width = 200) {
    p <- ggplot2::ggplot(nps_results, ggplot2::aes(Marker, NPS, fill = File)) +
        ggplot2::geom_bar(stat = "identity", position = "dodge") +
        ggplot2::labs(title = "Negative Population Spread (Unmixing Noise Floor)",
                      subtitle = "High MAD indicates unmixing-induced spreading error.\nThis is often caused by spectral overlap with bright markers in the sample.",
                      y = "Spread (MAD)", x = "Unmixed Marker") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    
    if (!is.null(output_file)) {
        ggplot2::ggsave(output_file, p, width = width, height = 120, units = "mm", dpi = 300)
    }
    return(p)
}
