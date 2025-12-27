#' Plot Detector-Level Residuals
#' 
#' Identifies which detectors contribute most to the unmixing error for high-RRMSE cells.
#' 
#' @param res_list List returned by calc_residuals with return_residuals=TRUE
#' @param top_n Number of top high-error cells to analyze
#' @param output_file Path to save the plot
#' @export
plot_detector_residuals <- function(res_list, top_n = 50, output_file = "detector_residuals.png") {
    data <- res_list$data
    residuals <- res_list$residuals
    
    # Identify high-error cells
    idx <- order(data$Relative_RMSE, decreasing = TRUE)[1:top_n]
    R_sub <- residuals[idx, , drop = FALSE]
    
    # Convert to long format for plotting
    long <- as.data.frame(R_sub)
    long$Cell <- seq_len(nrow(R_sub))
    long <- tidyr::pivot_longer(long, cols = -Cell, names_to = "Detector", values_to = "Residual")
    
    # Calculate median absolute residual per detector
    detector_order <- long |> 
        dplyr::group_by(Detector) |> 
        dplyr::summarize(Median_Abs_Res = median(abs(Residual))) |> 
        dplyr::arrange(desc(Median_Abs_Res))
    
    long$Detector <- factor(long$Detector, levels = detector_order$Detector)
    
    p <- ggplot2::ggplot(long, ggplot2::aes(Detector, Residual)) +
        ggplot2::geom_boxplot(outlier.size = 0.5, fill = "steelblue", alpha = 0.7) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        ggplot2::labs(title = paste("Residual Contributions for Top", top_n, "High-Error Cells"),
                      subtitle = "Positive residuals = under-explained signal; Negative = over-explained",
                      x = "Detector", y = "Residual Value") +
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
plot_nps <- function(nps_results, output_file = "nps_plot.png") {
    p <- ggplot2::ggplot(nps_results, ggplot2::aes(Marker, NPS, fill = File)) +
        ggplot2::geom_bar(stat = "identity", position = "dodge") +
        ggplot2::labs(title = "Negative Population Spread (Unmixing Noise Floor)",
                      y = "Spread (MAD)", x = "Unmixed Marker") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    
    if (!is.null(output_file)) {
        ggplot2::ggsave(output_file, p, width = width, height = 120, units = "mm", dpi = 300)
    }
    return(p)
}
