build_reference_matrix <- function(
  input_folder = "scc",
  output_folder = "autogate_plots",
  custom_fluorophores = NULL,
  default_sample_type = "beads",
  gating_opts = NULL,
  histogram_pct_beads = 0.98,
  histogram_direction_beads = "both",
  histogram_pct_cells = 0.35,
  histogram_direction_cells = "both",
  outlier_percentile = 0.02,
  debris_percentile = 0.02,
  bead_gate_scale = 1.3,
  histogram_min_x_log = 2,
  max_clusters = 6,
  min_cluster_proportion = 0.03,
  gate_contour_beads = 0.999999999,
  gate_contour_cells = 0.95,
  subsample_n = 5000
) {
    # If gating_opts provided, use those values
    if (!is.null(gating_opts)) {
        histogram_pct_beads <- gating_opts$histogram_pct_beads
        histogram_direction_beads <- gating_opts$histogram_direction_beads
        histogram_pct_cells <- gating_opts$histogram_pct_cells
        histogram_direction_cells <- gating_opts$histogram_direction_cells
    }

    library(mclust)
    library(data.table)
    library(sp)

    sample_patterns <- get_fluorophore_patterns()

    fcs_dir <- normalizePath(input_folder, mustWork = TRUE)
    fcs_files <- list.files(fcs_dir, pattern = "\\.fcs$", full.names = TRUE)

    out_path <- normalizePath(output_folder, mustWork = FALSE)
    dir.create(paste0(out_path, "/fsc_ssc"), showWarnings = FALSE, recursive = TRUE)
    dir.create(paste0(out_path, "/histogram"), showWarnings = FALSE, recursive = TRUE)
    dir.create(paste0(out_path, "/spectrum"), showWarnings = FALSE, recursive = TRUE)

    get_ellipse <- function(mean, sigma, level = 0.95, n = 100, scale = 1.0) {
        chi2_val <- qchisq(level, df = 2)
        eig <- eigen(sigma)
        a <- sqrt(eig$values[1] * chi2_val) * scale
        b <- sqrt(eig$values[2] * chi2_val) * scale
        angle <- atan2(eig$vectors[2, 1], eig$vectors[1, 1])
        theta <- seq(0, 2 * pi, length.out = n)
        ellipse_x <- a * cos(theta)
        ellipse_y <- b * sin(theta)
        rot_x <- ellipse_x * cos(angle) - ellipse_y * sin(angle)
        rot_y <- ellipse_x * sin(angle) + ellipse_y * cos(angle)
        data.table::data.table(x = rot_x + mean[1], y = rot_y + mean[2])
    }

    fit_gmm_populations <- function(data, max_k = 5, min_prop = 0.05) {
        fit <- mclust::Mclust(data, G = 1:max_k, verbose = FALSE)
        if (is.null(fit)) {
            return(NULL)
        }
        n_clusters <- fit$G
        proportions <- fit$parameters$pro
        means <- fit$parameters$mean
        
        # Robustly extract covariance matrices for each cluster
        sigmas <- lapply(1:n_clusters, function(k) {
            # mclust stores variances differently depending on the model
            if (fit$modelName %in% c("EII", "VII")) {
                # Spherical models: scalar variance * Identity
                return(diag(fit$parameters$variance$sigmasq[k], nrow = 2))
            } else {
                # Other models: sigma is a d x d x G array
                return(fit$parameters$variance$sigma[, , k])
            }
        })
        
        main_pops <- which(proportions >= min_prop)
        list(
            fit = fit, n_clusters = n_clusters, model_name = fit$modelName, BIC = fit$BIC,
            proportions = proportions, means = means, sigmas = sigmas,
            main_populations = main_pops, classification = fit$classification
        )
    }

    get_sample_type <- function(filename, patterns, default) {
        # Sort patterns by length (longest first) for each type
        for (type in names(patterns)) {
            type_patterns <- patterns[[type]]
            # Sort by length descending so longer patterns match first
            type_patterns <- type_patterns[order(-nchar(type_patterns))]
            for (pattern in type_patterns) {
                # Use fixed string matching (case insensitive)
                if (grepl(pattern, filename, fixed = FALSE, ignore.case = TRUE)) {
                    return(list(type = type, pattern = pattern))
                }
            }
        }
        list(type = default, pattern = "default")
    }

    select_bead_population <- function(gmm_result) {
        main_pops <- gmm_result$main_populations
        if (length(main_pops) == 0) {
            return(NULL)
        }
        proportions <- gmm_result$proportions[main_pops]
        best_pop <- main_pops[which.max(proportions)]
        list(selected = best_pop, reason = "largest population", proportion = gmm_result$proportions[best_pop])
    }

    select_cell_populations <- function(gmm_result, debris_threshold, detection_limit) {
        main_pops <- gmm_result$main_populations
        if (length(main_pops) == 0) {
            return(NULL)
        }
        valid_pops <- c()
        for (k in main_pops) {
            mean_fsc <- gmm_result$means[1, k]
            mean_ssc <- gmm_result$means[2, k]
            sigma <- gmm_result$sigmas[[k]]
            eig <- eigen(sigma)$values
            elongation <- max(eig) / min(eig)
            ell_extent <- sqrt(max(eig) * qchisq(0.95, df = 2))
            max_extent_fsc <- mean_fsc + ell_extent
            max_extent_ssc <- mean_ssc + ell_extent
            if (mean_fsc >= debris_threshold && mean_fsc <= detection_limit * 0.8 &&
                mean_ssc <= detection_limit * 0.8 && max_extent_fsc <= detection_limit &&
                max_extent_ssc <= detection_limit && elongation <= 50) {
                valid_pops <- c(valid_pops, k)
            }
        }
        list(selected = valid_pops, proportions = gmm_result$proportions[valid_pops])
    }

    create_merged_gate <- function(gmm_result, populations, level = 0.95, n_points = 100,
                                   clip_max_x = Inf, clip_max_y = Inf, scale = 1.0) {
        if (length(populations) == 0) {
            return(NULL)
        }
        if (length(populations) == 1) {
            ell <- get_ellipse(gmm_result$means[, populations], gmm_result$sigmas[[populations]], level, n_points, scale)
        } else {
            all_points <- data.table::rbindlist(lapply(populations, function(k) {
                get_ellipse(gmm_result$means[, k], gmm_result$sigmas[[k]], level, n_points, scale)
            }))
            hull_idx <- chull(all_points$x, all_points$y)
            ell <- all_points[c(hull_idx, hull_idx[1]), ]
        }
        ell[, x := pmin(x, clip_max_x)]
        ell[, y := pmin(y, clip_max_y)]
        ell[, x := pmax(x, 0)]
        ell[, y := pmax(y, 0)]
        ell
    }

    results_list <- list()

    for (fcs_file in fcs_files) {
        sn <- tools::file_path_sans_ext(basename(fcs_file))

        sample_info <- get_sample_type(sn, sample_patterns, default_sample_type)
        message("Processing: ", sn, " → type=", sample_info$type, ", pattern=", sample_info$pattern)

        fluor_name <- if (!is.null(custom_fluorophores) && sn %in% names(custom_fluorophores)) {
            custom_fluorophores[[sn]]
        } else {
            # Use expanded pattern list to get a clean fluorophore name
            found_name <- sample_info$pattern
            # Fallback to the pattern found by get_sample_type
            found_name
        }

        ff <- flowCore::read.FCS(fcs_file, transformation = FALSE, truncate_max_range = FALSE)
        pd <- flowCore::pData(flowCore::parameters(ff))

        fsc <- pd$name[grepl("^FSC", pd$name) & grepl("-A$", pd$name)][1]
        ssc <- pd$name[grepl("^SSC", pd$name) & grepl("-A$", pd$name)][1]
        fl_idx <- which(grepl("^FL[0-9]+-A$", pd$name))
        if (length(fl_idx) == 0) fl_idx <- which(grepl("^FL[0-9]+-Comp$", pd$name))
        fl_channels <- pd$name[fl_idx]

        data_raw <- flowCore::exprs(ff)[, c(fsc, ssc)]
        n_total <- nrow(data_raw)

        fsc_max <- quantile(data_raw[, 1], 1 - outlier_percentile, na.rm = TRUE)
        ssc_max <- quantile(data_raw[, 2], 1 - outlier_percentile, na.rm = TRUE)
        plot_max <- max(fsc_max, ssc_max) * 1.05

        valid_idx <- which(data_raw[, 1] < fsc_max & data_raw[, 2] < ssc_max &
            data_raw[, 1] > 0 & data_raw[, 2] > 0)
        data_filtered <- data_raw[valid_idx, ]

        if (sample_info$type %in% c("cells", "unstained")) {
            debris_threshold <- quantile(data_filtered[, 1], debris_percentile, na.rm = TRUE)
            debris_idx <- which(data_filtered[, 1] >= debris_threshold)
            data_filtered <- data_filtered[debris_idx, ]
        } else {
            debris_threshold <- 0
        }

        if (!is.null(subsample_n) && nrow(data_filtered) > subsample_n) {
            sample_idx <- sample(nrow(data_filtered), subsample_n)
            data_fit <- data_filtered[sample_idx, ]
        } else {
            data_fit <- data_filtered
        }

        gmm_result <- fit_gmm_populations(data_fit, max_k = max_clusters, min_prop = min_cluster_proportion)
        if (is.null(gmm_result)) {
            message("  SKIP: GMM failed")
            next
        }

        if (sample_info$type == "beads") {
            selection <- select_bead_population(gmm_result)
            selected_pops <- selection$selected
            gate_level <- gate_contour_beads
        } else {
            selection <- select_cell_populations(gmm_result, debris_threshold, min(fsc_max, ssc_max))
            selected_pops <- selection$selected
            gate_level <- gate_contour_cells
            if (length(selected_pops) == 0) {
                selection <- select_bead_population(gmm_result)
                selected_pops <- selection$selected
            }
        }
        if (length(selected_pops) == 0) {
            message("  SKIP: No populations selected")
            next
        }

        if (sample_info$type %in% c("cells", "unstained")) {
            final_gate <- create_merged_gate(gmm_result, selected_pops, gate_level,
                clip_max_x = fsc_max * 0.95, clip_max_y = ssc_max * 0.95
            )
        } else {
            final_gate <- create_merged_gate(gmm_result, selected_pops, gate_level, scale = bead_gate_scale)
        }

        all_data <- flowCore::exprs(ff)
        in_gate <- sp::point.in.polygon(all_data[, fsc], all_data[, ssc], final_gate$x, final_gate$y) > 0
        gated_data <- all_data[in_gate, ]
        n_gated <- nrow(gated_data)

        if (n_gated < 100 || length(fl_channels) == 0) {
            message("  SKIP: n_gated=", n_gated, " fl_channels=", length(fl_channels))
            next
        }

        q99 <- apply(gated_data[, fl_channels, drop = FALSE], 2, quantile, probs = 0.99)
        peak_channel <- fl_channels[which.max(q99)]
        peak_vals <- gated_data[, peak_channel]
        vals_log <- log10(pmax(peak_vals, 1))

        if (sample_info$type %in% c("unstained", "cells")) {
            vals_above <- vals_log
        } else {
            d_full <- density(vals_log, n = 512)
            mid_idx <- which(d_full$x > histogram_min_x_log & d_full$x < quantile(vals_log, 0.98))
            if (length(mid_idx) > 0) {
                local_min_idx <- mid_idx[which.min(d_full$y[mid_idx])]
                min_x <- d_full$x[local_min_idx]
            } else {
                min_x <- histogram_min_x_log
            }
            vals_above <- vals_log[vals_log >= min_x]
        }

        if (sample_info$type %in% c("unstained", "cells")) {
            pct <- histogram_pct_cells
            direction <- histogram_direction_cells
        } else {
            pct <- histogram_pct_beads
            direction <- histogram_direction_beads
        }

        if (direction == "right") {
            lower_q <- 0.5
            upper_q <- 0.5 + pct
        } else if (direction == "left") {
            lower_q <- 0.5 - pct
            upper_q <- 0.5
        } else {
            lower_q <- 0.5 - pct / 2
            upper_q <- 0.5 + pct / 2
        }
        lower_q <- max(0, lower_q)
        upper_q <- min(1, upper_q)

        gate_min <- 10^quantile(vals_above, lower_q)
        gate_max <- 10^quantile(vals_above, upper_q)

        in_hist_gate <- peak_vals >= gate_min & peak_vals <= gate_max
        final_gated_data <- gated_data[in_hist_gate, ]

        if (nrow(final_gated_data) < 10) next

        channels_desc <- pd$desc[match(fl_channels, pd$name)]
        channels_desc <- sub("-(A|H|W)$", "", sub("-Comp$", "", channels_desc))
        laser_nm <- as.integer(sub("^([0-9]+)nm.*$", "\\1", channels_desc))
        center_nm <- as.integer(sub("^([0-9]+).*$", "\\1", sub("^.*-\\s*", "", channels_desc)))
        ord <- order(laser_nm, center_nm)
        fl_channels_ord <- fl_channels[ord]

        mat_fl <- final_gated_data[, fl_channels_ord, drop = FALSE]
        spectrum <- apply(mat_fl, 2, median)
        spectrum_norm <- spectrum / max(spectrum)
        channels_desc_ord <- channels_desc[ord]
        laser_nm_ord <- laser_nm[ord]

        # FSC/SSC plot
        fsc_desc <- pd$desc[pd$name == fsc]
        ssc_desc <- pd$desc[pd$name == ssc]
        dt_plot <- data.table(x = data_raw[, 1], y = data_raw[, 2])
        x_lim <- c(0, plot_max)
        y_lim <- c(0, plot_max)
        x_breaks <- seq(0, plot_max, length.out = 201)
        y_breaks <- seq(0, plot_max, length.out = 201)
        dt_plot[, x_bin := findInterval(x, x_breaks, rightmost.closed = TRUE)]
        dt_plot[, y_bin := findInterval(y, y_breaks, rightmost.closed = TRUE)]
        dt2d <- dt_plot[x_bin >= 1 & x_bin <= 200 & y_bin >= 1 & y_bin <= 200, .(count = .N), by = .(x_bin, y_bin)]
        dt2d[, `:=`(x = x_breaks[x_bin], y = y_breaks[y_bin], fill = log10(count + 1))]

        p1 <- ggplot(dt2d, aes(x, y, fill = fill)) +
            geom_tile(width = diff(x_breaks)[1], height = diff(y_breaks)[1]) +
            scale_fill_gradientn(colors = c("#0000FF", "#00FFFF", "#00FF00", "#FFFF00", "#FF0000"), guide = "none") +
            geom_path(data = final_gate, aes(x, y), inherit.aes = FALSE, color = "red", linewidth = 1) +
            labs(
                title = paste0(sn, " - FSC/SSC"), subtitle = paste0(round(100 * n_gated / n_total, 1), "% gated"),
                x = fsc_desc, y = ssc_desc
            ) +
            theme_minimal(base_size = 11) +
            theme(
                legend.position = "none", panel.grid = element_blank(),
                panel.background = element_rect(fill = "white", color = NA),
                plot.background = element_rect(fill = "white", color = NA)
            ) +
            coord_cartesian(xlim = x_lim, ylim = y_lim)
        ggsave(paste0(out_path, "/fsc_ssc/", sn, "_fsc_ssc.png"), p1, width = 5, height = 5, dpi = 300)

        # Histogram plot
        peak_desc <- pd$desc[pd$name == peak_channel]
        gate_min_log <- log10(gate_min)
        gate_max_log <- log10(gate_max)
        n_hist_gated <- sum(in_hist_gate)

        p2 <- ggplot(data.table(x = vals_log), aes(x)) +
            geom_density(fill = "grey80", color = "grey40") +
            geom_vline(xintercept = gate_min_log, color = "red", linewidth = 1) +
            geom_vline(xintercept = gate_max_log, color = "red", linewidth = 1) +
            annotate("rect", xmin = gate_min_log, xmax = gate_max_log, ymin = -Inf, ymax = Inf, alpha = 0.15, fill = "red") +
            labs(
                title = paste0(sn, " - ", peak_desc),
                subtitle = paste0(round(100 * n_hist_gated / n_gated, 1), "% of FSC/SSC gated"),
                x = paste0("log10(", peak_desc, ")")
            ) +
            theme_minimal() +
            theme(legend.position = "none")
        ggsave(paste0(out_path, "/histogram/", sn, "_histogram.png"), p2, width = 6.5, height = 4, dpi = 300)

        # Spectrum plot
        log_mat <- log10(pmax(mat_fl, 1e-3))
        min_y <- floor(min(log_mat))
        max_y <- max(3, ceiling(max(log_mat)))
        breaks <- seq(min_y, max_y, length.out = 1001)
        bin_mid <- (breaks[-1] + breaks[-length(breaks)]) / 2
        bin_height <- breaks[2] - breaks[1]

        counts_mat <- vapply(seq_len(ncol(log_mat)), function(j) {
            as.numeric(hist(log_mat[, j], breaks = breaks, plot = FALSE)$counts)
        }, numeric(length(bin_mid)))
        rownames(counts_mat) <- as.character(seq_along(bin_mid))
        colnames(counts_mat) <- as.character(seq_len(ncol(log_mat)))

        dt <- as.data.table(as.table(counts_mat))
        setnames(dt, c("bin_idx", "ch_idx", "count"))
        dt[, bin_idx := as.integer(as.character(bin_idx))]
        dt[, ch_idx := as.integer(as.character(ch_idx))]
        dt[, y_orig := bin_mid[bin_idx]]
        dt[, fill := log10(count + 1)]
        dt <- dt[count >= 2]
        y_power <- 1.1
        dt[, y := y_orig^y_power]
        fill_lo <- min(dt$fill)
        fill_hi <- quantile(dt$fill, 0.87)
        vlines <- which(diff(laser_nm_ord) != 0) + 0.5
        y_breaks_orig <- 0:ceiling(max_y)
        y_breaks_trans <- y_breaks_orig^y_power
        y_labels <- sapply(y_breaks_orig, function(x) bquote(10^.(x)))

        p3 <- ggplot(dt, aes(ch_idx, y, fill = fill)) +
            geom_tile(width = 0.7, height = bin_height * 3) +
            scale_fill_gradientn(
                colors = colorRampPalette(c("#0000FF", "#00FFFF", "#00FF00", "#FFFF00", "#FF0000"))(100),
                limits = c(fill_lo, fill_hi), oob = scales::squish
            ) +
            {
                if (length(vlines) > 0) geom_vline(xintercept = vlines, color = "grey30", linewidth = 0.5)
            } +
            scale_x_continuous(breaks = seq_along(fl_channels_ord), labels = channels_desc_ord) +
            scale_y_continuous(limits = c(0, (max_y + 0.5)^y_power), breaks = y_breaks_trans, labels = y_labels) +
            coord_cartesian(expand = FALSE) +
            labs(
                title = paste0(sn, " - Spectrum (", format(nrow(final_gated_data), big.mark = ","), " events)"),
                x = NULL, y = "Intensity"
            ) +
            theme(
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
                axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12),
                axis.line = element_line(linewidth = 0.5), plot.title = element_text(size = 15),
                legend.position = "none", panel.grid.major = element_blank(),
                panel.background = element_rect(fill = "white")
            )
        ggsave(paste0(out_path, "/spectrum/", sn, "_spectrum.png"), p3, width = 300, height = 120, units = "mm", dpi = 600)

        results_list[[sn]] <- data.table::data.table(
            sample = sn,
            fluorophore = fluor_name,
            type = sample_info$type,
            n_total = n_total,
            n_final = nrow(final_gated_data),
            spectrum = list(spectrum_norm) # Store spectrum
        )
    }

    if (length(results_list) == 0) {
        warning("No samples were successfully processed")
        return(invisible(NULL))
    }

    results_dt <- data.table::rbindlist(results_list)
    data.table::fwrite(
        results_dt[, .(sample, fluorophore, type, n_total, n_final)],
        paste0(out_path, "/gating_summary.csv")
    )

    # Build reference matrix from collected spectra
    spectra_list <- results_dt$spectrum
    names(spectra_list) <- results_dt$fluorophore

    # Get detector names from first spectrum
    detector_names <- names(spectra_list[[1]])

    M <- do.call(rbind, lapply(spectra_list, function(s) {
        s_norm <- pmax(s, 0)
        s_norm / max(s_norm, na.rm = TRUE)
    }))
    rownames(M) <- names(spectra_list)
    colnames(M) <- detector_names

    # Save reference matrix
    M_df <- as.data.frame(M)
    M_df$file <- rownames(M_df)
    M_df <- M_df[, c("file", colnames(M))]
    data.table::fwrite(M_df, "reference_matrix.csv")

    message("Reference matrix saved to reference_matrix.csv")
    return(M)
}
