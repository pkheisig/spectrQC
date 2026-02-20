#' Plot FSC/SSC Scatter Colored by Unmixing Error
#'
#' Produces per-file FSC/SSC scatter panels colored by `RMSE_Score` or
#' `Relative_RMSE`, useful to detect whether high-error cells cluster in
#' specific morphological regions.
#'
#' @param data Unmixed data frame containing at least `FSC-A`, `SSC-A`, `File`,
#'   and the chosen metric column.
#' @param metric Error metric column name (`"RMSE_Score"` or `"Relative_RMSE"`).
#' @param output_file Optional output image path; set `NULL` to skip writing.
#' @param width Plot width for `ggsave()`.
#' @param height Plot height for `ggsave()`.
#' @param unit Plot size unit (for example `"mm"` or `"in"`).
#' @param dpi Output resolution.
#' @param max_cells_per_file Subsampling cap per file for rendering speed.
#' @param color_limits Optional fixed limits for the color scale.
#'
#' @return A `ggplot` object.
#' @export
#' @examples
#' \dontrun{
#' p <- plot_scatter_rmse(
#'   data = results_df,
#'   metric = "Relative_RMSE",
#'   output_file = "rrmse_scatter.png"
#' )
#' print(p)
#' }
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
        plot_data$Relative_RMSE <- plot_data$Relative_RMSE * 100
        legend_name <- "RRMSE (%)"
        if (is.null(color_limits)) color_limits <- c(0, 5) # Default 0-5% range for RRMSE
    }

    # Trim outliers for the plot coordinates
    # IMPORTANT: we arrange BY metric so high values are plotted LAST (on top)
    plot_data <- plot_data |>
        dplyr::filter(
            `FSC-A` >= quantile(`FSC-A`, 0.005, na.rm=TRUE) & `FSC-A` <= quantile(`FSC-A`, 0.995, na.rm=TRUE),
            `SSC-A` >= quantile(`SSC-A`, 0.005, na.rm=TRUE) & `SSC-A` <= quantile(`SSC-A`, 0.995, na.rm=TRUE)
        )
    plot_data <- plot_data[order(plot_data[[metric]]), , drop = FALSE]

    # Calculate number of columns for facet grid
    n_files <- length(unique(plot_data$File))
    ncols <- ceiling(sqrt(n_files))

    p <- ggplot2::ggplot(plot_data, ggplot2::aes_string(x = "FSC-A", y = "SSC-A", color = metric)) +
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
#' @param width Plot width.
#' @param height Plot height.
#' @param unit Plot size unit.
#' @param dpi Output resolution.
#' @param max_cells Maximum events sampled before plotting.
#' @param show_smooth Logical. If TRUE, adds a red trend line (GAM).
#' @param y_limits Optional y-axis limits
#' @return A `ggplot` object.
#' @export
#' @examples
#' \dontrun{
#' p <- plot_marker_correlations(
#'   data = results_df,
#'   metric = "Relative_RMSE",
#'   output_file = "marker_rrmse_correlations.png"
#' )
#' print(p)
#' }
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
        plot_data$Relative_RMSE <- plot_data$Relative_RMSE * 100
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

    p <- ggplot2::ggplot(long, ggplot2::aes_string(x = "Intensity", y = metric)) +
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
#' @param width Plot width.
#' @param height Plot height.
#' @param unit Plot size unit.
#' @param dpi Output resolution.
#' @param theme_custom Optional ggplot theme object (reserved).
#' @return A `ggplot` object.
#' @export
#' @examples
#' m <- diag(3)
#' rownames(m) <- c("FITC", "PE", "APC")
#' colnames(m) <- c("B2-A", "YG1-A", "R1-A")
#' p <- plot_spectra(ref_matrix = m, output_file = NULL)
#' print(p)
plot_spectra <- function(ref_matrix,
                         pd = NULL,
                         output_file = "spectra_overlay.png",
                         width = 250,
                         height = 100,
                         unit = "mm",
                         dpi = 600,
                         theme_custom = NULL) {
    ref_matrix <- .as_reference_matrix(ref_matrix, "ref_matrix")
    detectors <- colnames(ref_matrix)
    
    # 1. Get Sorted Detectors and Labels
    if (!is.null(pd)) {
        det_info <- get_sorted_detectors(pd)
        # Filter ref_matrix to match
        common <- intersect(det_info$names, detectors)
        
        if (length(common) == 0) {
            warning("No matching detectors found between reference matrix and provided metadata. Ignoring metadata.")
            message("Ref Matrix cols (first 5): ", paste(utils::head(detectors, 5), collapse=", "))
            message("Metadata names (first 5): ", paste(utils::head(det_info$names, 5), collapse=", "))

            # Fallback to numerical sort
            nums <- as.numeric(gsub("[^0-9]", "", detectors))
            ord <- order(nums)
            ref_matrix <- ref_matrix[, ord, drop = FALSE]
            detectors <- colnames(ref_matrix)
            labels <- detectors
        } else {
            ref_matrix <- ref_matrix[, common, drop = FALSE]
            detectors <- common
            # Get labels only for common
            labels <- det_info$labels[match(common, det_info$names)]
        }
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


#' Plot NxN Unmixing Scatter Matrix
#'
#' Creates a lower-triangle marker-vs-marker scatter matrix from unmixed SCC output.
#' Each populated panel uses events from only one single-stain file (row stain),
#' while the column stain is permuted across the lower triangle.
#'
#' @param unmixed_list Output list from `unmix_samples()`.
#' @param sample_to_marker Named character vector mapping sample IDs (basename
#'   without extension) to fluorophore/marker names.
#' @param markers Optional character vector of marker columns to include.
#' @param output_file Path to PNG output.
#' @param max_points_per_sample Maximum events sampled per sample for plotting.
#' @param transform One of `"asinh"` or `"none"`.
#' @param asinh_cofactor Cofactor used when `transform = "asinh"`.
#' @param panel_size_mm Size per matrix panel in millimeters.
#' @return Invisibly returns `output_file`.
#' @export
#' @examples
#' \dontrun{
#' plot_unmixing_scatter_matrix(
#'   unmixed_list = unmixed_list,
#'   markers = rownames(M),
#'   output_file = "scc_unmixing_scatter_matrix.png"
#' )
#' }
plot_unmixing_scatter_matrix <- function(
    unmixed_list,
    sample_to_marker = NULL,
    markers = NULL,
    output_file = "unmixing_scatter_matrix.png",
    max_points_per_sample = 3000,
    transform = c("none", "asinh"),
    asinh_cofactor = 150,
    panel_size_mm = 30
) {
    transform <- match.arg(transform)
    if (!is.list(unmixed_list) || length(unmixed_list) == 0) {
        stop("unmixed_list must be a non-empty list from unmix_samples().")
    }

    extract_data <- function(x) {
        if (is.list(x) && "data" %in% names(x) && is.data.frame(x$data)) return(x$data)
        if (is.data.frame(x)) return(x)
        NULL
    }

    data_list <- lapply(unmixed_list, extract_data)
    valid <- !vapply(data_list, is.null, logical(1))
    data_list <- data_list[valid]
    if (length(data_list) == 0) stop("No valid data frames found in unmixed_list.")

    sample_ids <- names(data_list)
    if (is.null(sample_ids) || any(sample_ids == "")) {
        sample_ids <- paste0("sample_", seq_along(data_list))
    }
    names(data_list) <- sample_ids

    normalize_id <- function(x) tools::file_path_sans_ext(basename(as.character(x)))
    sample_keys <- normalize_id(sample_ids)

    if (is.null(sample_to_marker)) {
        sample_to_marker <- stats::setNames(sample_ids, sample_ids)
    }
    marker_names_raw <- names(sample_to_marker)
    sample_to_marker <- as.character(sample_to_marker)
    names(sample_to_marker) <- marker_names_raw
    names(sample_to_marker) <- normalize_id(names(sample_to_marker))
    sample_to_marker <- sample_to_marker[!is.na(names(sample_to_marker)) & names(sample_to_marker) != ""]
    sample_to_marker <- sample_to_marker[!duplicated(names(sample_to_marker))]

    sample_stains <- sample_to_marker[sample_keys]
    names(sample_stains) <- sample_ids

    if (is.null(markers)) {
        markers <- unique(sample_stains[!is.na(sample_stains) & sample_stains != ""])
    } else {
        markers <- unique(as.character(markers))
    }
    markers <- markers[!is.na(markers) & markers != ""]
    if (length(markers) < 2) stop("Need at least two marker names for scatter matrix.")

    panel_data <- list()
    panel_limits <- list()
    k <- 1
    lim_k <- 1

    compute_limits <- function(v) {
        v <- v[is.finite(v)]
        if (length(v) == 0) return(c(-1, 1))
        q <- stats::quantile(v, probs = c(0.005, 0.995), na.rm = TRUE, names = FALSE, type = 7)
        lo <- q[1]
        hi <- q[2]
        if (!is.finite(lo) || !is.finite(hi)) {
            lo <- min(v, na.rm = TRUE)
            hi <- max(v, na.rm = TRUE)
        }
        lo <- min(lo, 0)
        hi <- max(hi, 0)
        if (lo == hi) {
            pad <- max(1e-6, abs(lo) * 0.1)
            lo <- lo - pad
            hi <- hi + pad
        } else {
            pad <- (hi - lo) * 0.05
            lo <- lo - pad
            hi <- hi + pad
        }
        c(lo, hi)
    }

    for (i in seq_along(data_list)) {
        d <- data_list[[i]]
        row_stain <- sample_stains[[i]]
        if (is.na(row_stain) || row_stain == "" || !(row_stain %in% markers)) next
        if (!(row_stain %in% colnames(d))) next

        marker_idx <- match(row_stain, markers)
        if (is.na(marker_idx) || marker_idx <= 1) next
        x_markers <- markers[seq_len(marker_idx - 1)]
        x_markers <- x_markers[x_markers %in% colnames(d)]
        if (length(x_markers) == 0) next

        for (xm in x_markers) {
            d_pair <- d[, c(xm, row_stain), drop = FALSE]
            d_pair <- d_pair[stats::complete.cases(d_pair), , drop = FALSE]
            if (nrow(d_pair) == 0) next
            if (nrow(d_pair) > max_points_per_sample) {
                d_pair <- d_pair[sample.int(nrow(d_pair), max_points_per_sample), , drop = FALSE]
            }

            x_vals <- d_pair[[xm]]
            y_vals <- d_pair[[row_stain]]
            if (transform == "asinh") {
                x_vals <- asinh(x_vals / asinh_cofactor)
                y_vals <- asinh(y_vals / asinh_cofactor)
            }

            x_lim <- compute_limits(x_vals)
            y_lim <- compute_limits(y_vals)
            x_plot <- pmax(pmin(x_vals, x_lim[2]), x_lim[1])
            y_plot <- pmax(pmin(y_vals, y_lim[2]), y_lim[1])

            panel_data[[k]] <- data.frame(
                x = x_plot,
                y = y_plot,
                panel_col = xm,
                panel_row = row_stain,
                stringsAsFactors = FALSE
            )
            panel_limits[[lim_k]] <- data.frame(
                panel_col = xm,
                panel_row = row_stain,
                x_low = x_lim[1],
                x_high = x_lim[2],
                y_low = y_lim[1],
                y_high = y_lim[2],
                stringsAsFactors = FALSE
            )
            k <- k + 1
            lim_k <- lim_k + 1
        }
    }
    if (length(panel_data) == 0) stop("No panel data available for scatter matrix plot.")
    plot_df <- do.call(rbind, panel_data)
    lim_df <- do.call(rbind, panel_limits)

    plot_df$panel_col <- factor(plot_df$panel_col, levels = markers)
    plot_df$panel_row <- factor(plot_df$panel_row, levels = markers)
    lim_df$panel_col <- factor(lim_df$panel_col, levels = markers)
    lim_df$panel_row <- factor(lim_df$panel_row, levels = markers)

    panel_used <- unique(lim_df[, c("panel_row", "panel_col"), drop = FALSE])
    panel_grid <- expand.grid(
        panel_row = factor(markers, levels = markers),
        panel_col = factor(markers, levels = markers),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )

    p <- ggplot2::ggplot() +
        ggplot2::geom_blank(data = panel_grid, ggplot2::aes(x = 0, y = 0)) +
        ggplot2::geom_blank(data = lim_df, ggplot2::aes(x = x_low, y = y_low)) +
        ggplot2::geom_blank(data = lim_df, ggplot2::aes(x = x_high, y = y_high)) +
        ggplot2::geom_hline(
            data = panel_used,
            ggplot2::aes(yintercept = 0),
            inherit.aes = FALSE,
            color = "grey45",
            linewidth = 0.1
        ) +
        ggplot2::geom_vline(
            data = panel_used,
            ggplot2::aes(xintercept = 0),
            inherit.aes = FALSE,
            color = "grey45",
            linewidth = 0.1
        ) +
        ggplot2::geom_point(
            data = plot_df,
            ggplot2::aes(x = x, y = y),
            color = "black",
            alpha = 0.6,
            size = 0.25,
            stroke = 0
        ) +
        ggplot2::facet_grid(panel_row ~ panel_col, drop = FALSE, switch = "y", scales = "free") +
        ggplot2::labs(title = "Unmixing Scatter Matrix") +
        ggplot2::theme_bw(base_size = 7) +
        ggplot2::theme(
            legend.position = "none",
            panel.grid = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(color = "grey60", linewidth = 0.25),
            strip.background = ggplot2::element_rect(fill = "grey95", color = "grey80"),
            strip.text = ggplot2::element_text(size = 6, face = "bold"),
            axis.title = ggplot2::element_blank(),
            axis.text = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(),
            axis.line = ggplot2::element_blank(),
            panel.spacing = grid::unit(0.3, "mm"),
            plot.title = ggplot2::element_text(size = 10, face = "bold", hjust = 0.5)
        )

    n <- length(markers)
    panel_size_mm <- as.numeric(panel_size_mm)[1]
    if (!is.finite(panel_size_mm) || panel_size_mm <= 0) panel_size_mm <- 30
    width_mm <- max(220, panel_size_mm * n + 35)
    height_mm <- max(220, panel_size_mm * n + 30)
    ggplot2::ggsave(output_file, p, width = width_mm, height = height_mm, units = "mm", dpi = 350)

    invisible(output_file)
}
