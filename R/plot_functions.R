plot_scatter_rmse <- function(data,
                              metric = "RMSE_Score",
                              output_file = "scatter_rmse.png",
                              width = 180, # 1.5x larger default
                              height = 180, 
                              unit = "mm",
                              dpi = 600,
                              max_cells_per_file = 2000,
                              color_limits = NULL) {
    # Check required columns
    if (!all(c("FSC-A", "SSC-A", metric, "File") %in% colnames(data))) {
        stop(paste("Data must contain FSC-A, SSC-A, File and", metric, "columns"))
    }

    # Subsample per file for speed
    data <- data |>
        dplyr::group_by(File) |>
        dplyr::slice_sample(n = max_cells_per_file) |>
        dplyr::ungroup()
        
    plot_data <- data.table::copy(data.table::as.data.table(data))
    
    # If Relative_RMSE, convert to percentage for better readability
    legend_name <- metric
    if (metric == "Relative_RMSE") {
        plot_data[, Relative_RMSE := Relative_RMSE * 100]
        legend_name <- "RRMSE (%)"
        if (is.null(color_limits)) color_limits <- c(0, 5) # Default 0-5% range for RRMSE
    }

    # Trim outliers for the plot coordinates
    # IMPORTANT: we arrange BY metric so high values are plotted LAST (on top)
    plot_data <- plot_data |>
        dplyr::filter(
            `FSC-A` >= quantile(`FSC-A`, 0.005, na.rm=TRUE) & `FSC-A` <= quantile(`FSC-A`, 0.995, na.rm=TRUE),
            `SSC-A` >= quantile(`SSC-A`, 0.005, na.rm=TRUE) & `SSC-A` <= quantile(`SSC-A`, 0.995, na.rm=TRUE)
        ) |>
        dplyr::arrange(.data[[metric]])

    # Calculate number of columns for facet grid
    n_files <- length(unique(plot_data$File))
    ncols <- ceiling(sqrt(n_files))

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(`FSC-A`, `SSC-A`, color = .data[[metric]])) +
        ggplot2::geom_point(size = 0.1, alpha = 0.8) + # Larger and more opaque
        ggplot2::scale_color_gradientn(
            colors = c("gray98", "gray90", "orange", "red", "black"),
            name = legend_name,
            limits = color_limits,
            oob = scales::squish
        ) +
        ggplot2::facet_wrap(~File, scales = "fixed", ncol = ncols) +
        ggplot2::labs(x = "FSC-A", y = "SSC-A") +
        ggplot2::theme_minimal(base_size = 8) +
        ggplot2::theme(
            aspect.ratio = 1, # Force square facets
            panel.grid = ggplot2::element_blank(),
            panel.background = ggplot2::element_rect(fill = "white", color = "gray95"),
            strip.text = ggplot2::element_text(size = 5),
            axis.text = ggplot2::element_text(size = 5),
            axis.title = ggplot2::element_text(size = 6),
            legend.title = ggplot2::element_text(size = 6),
            legend.text = ggplot2::element_text(size = 5)
        )

    if (!is.null(output_file)) {
        ggplot2::ggsave(output_file, p, width = width, height = height, units = unit, dpi = dpi)
    }
    return(p)
}



#' Plot Marker-RRMSE Correlations
#' 
#' Plots the unmixed intensity of each marker against the RRMSE.
#' The red line is a GAM fit representing the trend of unmixing error vs signal.
#' Ideally, the line should be flat or decreasing. A sharp upward trend indicates 
#' detector non-linearity or spectral mismatch for that specific fluorophore.
#' 
#' @param data Unmixed results data.table
#' @param metric Metric to plot ("RMSE_Score" or "Relative_RMSE")
#' @param markers Vector of marker names
#' @param output_file Path to save plot
#' @param show_smooth Logical. If TRUE, adds a red trend line (GAM).
#' @param y_limits Optional y-axis limits
#' @export
plot_marker_correlations <- function(data, 
                                     metric = "RMSE_Score",
                                     markers = NULL,
                                     output_file = "marker_correlations.png",
                                     width = 150,
                                     height = 100,
                                     unit = "mm",
                                     dpi = 600,
                                     max_cells = 5000,
                                     y_limits = NULL,
                                     show_smooth = TRUE) {
    if (is.null(markers)) {
        exclude_cols <- c("RMSE_Score", "Relative_RMSE", "File", "FSC-A", "SSC-A", "FSC-H", "SSC-H")
        markers <- setdiff(colnames(data), exclude_cols)
    }

    # Subsample for speed
    if (nrow(data) > max_cells) {
        data <- data[sample(nrow(data), max_cells), ]
    }
    
    plot_data <- data.table::copy(data.table::as.data.table(data))
    y_name <- metric
    if (metric == "Relative_RMSE") {
        plot_data[, Relative_RMSE := Relative_RMSE * 100]
        y_name <- "RRMSE (%)"
        if (is.null(y_limits)) y_limits <- c(0, 10) # Default 0-10% for correlations
    }

    long <- tidyr::pivot_longer(plot_data,
        cols = dplyr::all_of(markers),
        names_to = "Marker", values_to = "Intensity"
    )

    # Trim top/bottom 0.5% per marker
    long <- long |>
        dplyr::group_by(Marker) |>
        dplyr::filter(
            Intensity >= quantile(Intensity, 0.005, na.rm=TRUE) &
                Intensity <= quantile(Intensity, 0.995, na.rm=TRUE)
        ) |>
        dplyr::ungroup()

    p <- ggplot2::ggplot(long, ggplot2::aes(Intensity, .data[[metric]])) +
        ggplot2::geom_point(size = 0.1, alpha = 0.1) +
        {if(metric == "Relative_RMSE") ggplot2::geom_hline(yintercept = 5, color = "darkred", linetype = "dashed", linewidth = 0.8, alpha = 0.7)} +
        {if(show_smooth) ggplot2::geom_smooth(method = "gam", color = "red", linewidth = 0.5, formula = y ~ s(x, bs = "cs"))} +
        ggplot2::facet_wrap(~Marker, scales = "free_x", ncol = ceiling(length(markers) / 3)) +
        ggplot2::labs(x = "Unmixed Abundance", y = y_name) +
        ggplot2::theme_minimal(base_size = 8) +
        ggplot2::theme(
            aspect.ratio = 1, # Make square
            panel.grid.minor = ggplot2::element_blank(),
            strip.text = ggplot2::element_text(size = 5),
            axis.text = ggplot2::element_text(size = 5),
            axis.title = ggplot2::element_text(size = 6)
        )
    
    if (!is.null(y_limits)) {
        p <- p + ggplot2::coord_cartesian(ylim = y_limits)
    }

    if (!is.null(output_file)) {
        ggplot2::ggsave(output_file, p, width = width, height = height, units = unit, dpi = dpi)
    }
    return(p)
}


#' Plot Spectral Overlays
#' 
#' @param ref_matrix Reference matrix (Markers x Detectors)
#' @param pd Optional pData for descriptive labels
#' @param output_file Path to save plot
#' @export
plot_spectra <- function(ref_matrix,
                         pd = NULL,
                         output_file = "spectra_overlay.png",
                         width = 250,
                         height = 100,
                         unit = "mm",
                         dpi = 600,
                         theme_custom = NULL) {
    detectors <- colnames(ref_matrix)
    
    # 1. Get Sorted Detectors and Labels
    if (!is.null(pd)) {
        det_info <- get_sorted_detectors(pd)
        # Filter ref_matrix to match
        common <- intersect(det_info$names, detectors)
        ref_matrix <- ref_matrix[, common, drop = FALSE]
        detectors <- common
        # Get labels only for common
        labels <- det_info$labels[match(common, det_info$names)]
    } else {
        # Fallback to numerical sort
        nums <- as.numeric(gsub("[^0-9]", "", detectors))
        ord <- order(nums)
        ref_matrix <- ref_matrix[, ord, drop = FALSE]
        detectors <- colnames(ref_matrix)
        labels <- detectors
    }

    long <- data.frame(
        Fluorophore = rep(rownames(ref_matrix), ncol(ref_matrix)),
        Detector = rep(detectors, each = nrow(ref_matrix)),
        Intensity = as.vector(ref_matrix)
    )
    
    long$Detector <- factor(long$Detector, levels = detectors, labels = labels)

    p <- ggplot2::ggplot(long, ggplot2::aes(Detector, Intensity, color = Fluorophore, group = Fluorophore)) +
        ggplot2::geom_line(linewidth = 0.7) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5)) +
        ggplot2::labs(x = "Detector", y = "Normalized Intensity")

    if (!is.null(output_file)) {
        ggplot2::ggsave(output_file, p, width = width, height = height, unit = unit, dpi = dpi)
    }
    return(p)
}
