build_reference_matrix <- function(
  input_folder = "scc",
  output_folder = "gating_and_spectrum_plots",
  control_df = NULL,
  include_multi_af = FALSE,
  af_dir = "af",
  default_sample_type = "beads",
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
    library(mclust)
    library(data.table)
    library(sp)
    library(ggplot2)

    sample_patterns <- get_fluorophore_patterns()
    fcs_files <- list.files(input_folder, pattern = "\\.fcs$", full.names = TRUE)
    out_path <- normalizePath(output_folder, mustWork = FALSE)
    dir.create(file.path(out_path, "fsc_ssc"), showWarnings = FALSE, recursive = TRUE)
    dir.create(file.path(out_path, "histogram"), showWarnings = FALSE, recursive = TRUE)
    dir.create(file.path(out_path, "spectrum"), showWarnings = FALSE, recursive = TRUE)

    # 1. Establish Detector Sorting and Labels
    ff_meta <- flowCore::read.FCS(fcs_files[1], transformation = FALSE, truncate_max_range = FALSE)
    det_info <- get_sorted_detectors(flowCore::pData(flowCore::parameters(ff_meta)))
    detector_names <- det_info$names
    detector_labels <- det_info$labels

    message("Found ", length(detector_names), " spectral detectors. Sorting by laser...")

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
        sigmas <- lapply(1:fit$G, function(k) {
            if (fit$modelName %in% c("EII", "VII")) {
                diag(fit$parameters$variance$sigmasq[k], nrow = 2)
            } else {
                fit$parameters$variance$sigma[, , k]
            }
        })
        list(fit = fit, proportions = fit$parameters$pro, means = fit$parameters$mean, sigmas = sigmas, main_populations = which(fit$parameters$pro >= min_prop))
    }

    get_sample_type <- function(filename, patterns, default) {
        for (type in names(patterns)) {
            pats <- patterns[[type]][order(-nchar(patterns[[type]]))]
            for (p in pats) {
                if (grepl(p, filename, fixed = FALSE, ignore.case = TRUE)) {
                    return(list(type = type, pattern = p))
                }
            }
        }
        list(type = default, pattern = "default")
    }

    select_bead_population <- function(gmm_result) {
        if (length(gmm_result$main_populations) == 0) {
            return(NULL)
        }
        best <- gmm_result$main_populations[which.max(gmm_result$proportions[gmm_result$main_populations])]
        list(selected = best)
    }

    select_cell_populations <- function(gmm_result, debris_threshold, detection_limit) {
        valid <- c()
        for (k in gmm_result$main_populations) {
            if (gmm_result$means[1, k] >= debris_threshold && gmm_result$means[1, k] <= detection_limit * 0.8) valid <- c(valid, k)
        }
        list(selected = valid)
    }

    create_merged_gate <- function(gmm_result, populations, level, scale = 1.0, clip_x = Inf, clip_y = Inf) {
        if (length(populations) == 0) {
            return(NULL)
        }
        if (length(populations) == 1) {
            ell <- get_ellipse(gmm_result$means[, populations], gmm_result$sigmas[[populations]], level, scale = scale)
        } else {
            all_pts <- data.table::rbindlist(lapply(populations, function(k) get_ellipse(gmm_result$means[, k], gmm_result$sigmas[[k]], level, scale = scale)))
            ell <- all_pts[chull(x, y), ]
        }
        ell[, x := pmax(0, pmin(x, clip_x))]
        ell[, y := pmax(0, pmin(y, clip_y))]
        return(ell)
    }

    results_list <- list()

    # 1.1 Add files from AF directory if requested
    fcs_files_all <- fcs_files
    if (include_multi_af && dir.exists(af_dir)) {
        af_files <- list.files(af_dir, pattern = "\\.fcs$", full.names = TRUE)
        message("Found ", length(af_files), " extra AF files in '", af_dir, "'")
        # Pre-pend AF files so they are processed
        fcs_files_all <- c(af_files, fcs_files)
    }

    af_data_raw <- NULL
    af_fn <- NULL
    if (!is.null(control_df)) {
        af_rows <- control_df[fluorophore == "AF"]
        if (nrow(af_rows) > 0) af_fn <- tools::file_path_sans_ext(basename(af_rows$filename[1]))
    }
    if (is.null(af_fn)) {
        af_idx_tmp <- grep("Unstained|US_UT", fcs_files, ignore.case = TRUE)
        if (length(af_idx_tmp) > 0) af_fn <- tools::file_path_sans_ext(basename(fcs_files[af_idx_tmp[1]]))
    }
    if (!is.null(af_fn)) {
        af_path <- fcs_files[grep(af_fn, fcs_files, fixed = TRUE)]
        if (length(af_path) > 0) {
            ff_af <- flowCore::read.FCS(af_path[1], transformation = FALSE, truncate_max_range = FALSE)
            af_data_raw <- apply(flowCore::exprs(ff_af)[, detector_names, drop = FALSE], 2, median)
        }
    }

    for (fcs_file in fcs_files_all) {
        sn_ext <- basename(fcs_file)
        sn <- tools::file_path_sans_ext(sn_ext)

        # Determine if this is an extra AF file
        is_extra_af <- FALSE
        if (include_multi_af && grepl(normalizePath(af_dir), normalizePath(fcs_file), fixed = TRUE)) {
            is_extra_af <- TRUE
        }

        row_info <- if (!is.null(control_df)) control_df[filename == sn_ext | filename == sn] else data.table::data.table()
        sample_info <- get_sample_type(sn, sample_patterns, default_sample_type)

        fluor_name <- if (nrow(row_info) > 0 && !is.na(row_info$fluorophore[1])) row_info$fluorophore[1] else sample_info$pattern

        if (is_extra_af) {
            # Assign a unique AF name if not in control_df
            if (nrow(row_info) == 0) {
                fluor_name <- paste0("AF_", sn)
            }
            sample_info$type <- "cells" # Force cells for AF folder
        }

        if (fluor_name == "AF" && !is_extra_af) next

        message("Processing SCC: ", fluor_name, " (", sn, ")")
        ff <- flowCore::read.FCS(fcs_file, transformation = FALSE, truncate_max_range = FALSE)
        pd <- flowCore::pData(flowCore::parameters(ff))
        raw_data <- flowCore::exprs(ff)
        fsc <- pd$name[grepl("^FSC", pd$name) & grepl("-A$", pd$name)][1]
        ssc <- pd$name[grepl("^SSC", pd$name) & grepl("-A$", pd$name)][1]

        data_raw_scatter <- raw_data[, c(fsc, ssc)]
        fsc_max <- quantile(data_raw_scatter[, 1], 1 - outlier_percentile, na.rm = TRUE)
        ssc_max <- quantile(data_raw_scatter[, 2], 1 - outlier_percentile, na.rm = TRUE)
        valid_idx <- which(data_raw_scatter[, 1] < fsc_max & data_raw_scatter[, 2] < ssc_max & data_raw_scatter[, 1] > 0 & data_raw_scatter[, 2] > 0)
        data_filtered <- data_raw_scatter[valid_idx, ]
        if (sample_info$type %in% c("cells", "unstained")) {
            debris_threshold <- quantile(data_filtered[, 1], debris_percentile, na.rm = TRUE)
            data_filtered <- data_filtered[data_filtered[, 1] >= debris_threshold, ]
        } else {
            debris_threshold <- 0
        }

        if (!is.null(subsample_n) && nrow(data_filtered) > subsample_n) {
            data_fit <- data_filtered[sample(nrow(data_filtered), subsample_n), ]
        } else {
            data_fit <- data_filtered
        }

        gmm_result <- fit_gmm_populations(data_fit, max_k = max_clusters, min_prop = min_cluster_proportion)
        if (is.null(gmm_result)) next
        if (sample_info$type == "beads") {
            selected_pops <- select_bead_population(gmm_result)$selected
            gate_level <- gate_contour_beads
        } else {
            selected_pops <- select_cell_populations(gmm_result, debris_threshold, min(fsc_max, ssc_max))$selected
            gate_level <- gate_contour_cells
            if (length(selected_pops) == 0) selected_pops <- select_bead_population(gmm_result)$selected
        }
        if (length(selected_pops) == 0) next
        final_gate <- create_merged_gate(gmm_result, selected_pops, gate_level, scale = if (sample_info$type == "beads") bead_gate_scale else 1.0, clip_x = fsc_max * 0.95, clip_y = ssc_max * 0.95)

        gated_data <- raw_data[sp::point.in.polygon(raw_data[, fsc], raw_data[, ssc], final_gate$x, final_gate$y) > 0, ]
        if (nrow(gated_data) < 100) next

        peak_channel <- if (nrow(row_info) > 0 && !is.na(row_info$channel[1]) && row_info$channel[1] != "") row_info$channel[1] else detector_names[which.max(apply(gated_data[, detector_names, drop = FALSE], 2, function(x) quantile(x, 0.999)))]
        peak_vals <- gated_data[, peak_channel]
        vals_log <- log10(pmax(peak_vals, 1))
        if (sample_info$type %in% c("unstained", "cells")) {
            vals_above <- vals_log
        } else {
            d_full <- density(vals_log, n = 512)
            mid_idx <- which(d_full$x > histogram_min_x_log & d_full$x < quantile(vals_log, 0.98))
            min_x <- if (length(mid_idx) > 0) d_full$x[mid_idx[which.min(d_full$y[mid_idx])]] else histogram_min_x_log
            vals_above <- vals_log[vals_log >= min_x]
        }
        pct <- if (sample_info$type %in% c("unstained", "cells")) histogram_pct_cells else histogram_pct_beads
        dir <- if (sample_info$type %in% c("unstained", "cells")) histogram_direction_cells else histogram_direction_beads
        if (dir == "right") {
            lq <- 0.5
            uq <- 0.5 + pct
        } else if (dir == "left") {
            lq <- 0.5 - pct
            uq <- 0.5
        } else {
            lq <- 0.5 - pct / 2
            uq <- 0.5 + pct / 2
        }
        gate_min <- 10^quantile(vals_above, max(0, lq))
        gate_max <- 10^quantile(vals_above, min(1, uq))
        final_gated_data <- gated_data[peak_vals >= gate_min & peak_vals <= gate_max, ]
        neg_gated_data <- gated_data[peak_vals <= 10^quantile(vals_log, 0.15), ]
        if (nrow(final_gated_data) < 10) next

        pos_spectrum_raw <- apply(final_gated_data[, detector_names, drop = FALSE], 2, median)
        neg_spectrum_raw <- apply(gated_data[peak_vals <= 10^quantile(vals_log, 0.15), detector_names, drop = FALSE], 2, median)
        use_univ <- if (nrow(row_info) > 0 && !is.na(row_info$universal.negative[1])) (row_info$universal.negative[1] %in% c("TRUE", TRUE, "AF")) else FALSE
        final_neg <- if (use_univ && !is.null(af_data_raw)) af_data_raw else neg_spectrum_raw
        sig_pure <- pmax(pos_spectrum_raw - final_neg, 0)
        if (max(sig_pure) <= 0) sig_pure <- pmax(pos_spectrum_raw, 0)
        spectrum_norm <- sig_pure / max(sig_pure)

        # FSC/SSC plot
        fsc_desc <- pd$desc[pd$name == fsc]
        ssc_desc <- pd$desc[pd$name == ssc]
        dt_plot_gating <- data.table::data.table(x = raw_data[, fsc], y = raw_data[, ssc])
        x_breaks <- seq(0, max(fsc_max, ssc_max) * 1.05, length.out = 201)
        y_breaks <- x_breaks
        dt2d <- dt_plot_gating[findInterval(x, x_breaks) >= 1 & findInterval(x, x_breaks) <= 200 & findInterval(y, x_breaks) >= 1 & findInterval(y, x_breaks) <= 200, .(count = .N), by = .(x_bin = findInterval(x, x_breaks), y_bin = findInterval(y, x_breaks))]
        dt2d[, `:=`(x = x_breaks[x_bin], y = x_breaks[y_bin], fill = log10(count + 1))]
        p1 <- ggplot2::ggplot(dt2d, ggplot2::aes(x, y, fill = fill)) +
            ggplot2::geom_tile(width = diff(x_breaks)[1], height = diff(y_breaks)[1]) +
            ggplot2::scale_fill_gradientn(colors = c("#0000FF", "#00FFFF", "#00FF00", "#FFFF00", "#FF0000"), guide = "none") +
            ggplot2::geom_path(data = final_gate, ggplot2::aes(x, y), inherit.aes = FALSE, color = "red", linewidth = 1) +
            ggplot2::labs(title = paste0(sn, " - FSC/SSC"), subtitle = paste0(round(100 * nrow(gated_data) / nrow(raw_data), 1), "% gated"), x = pd$desc[pd$name == fsc], y = pd$desc[pd$name == ssc]) +
            ggplot2::theme_minimal() +
            ggplot2::theme(legend.position = "none", panel.grid = element_blank(), panel.background = element_rect(fill = "white", color = NA)) +
            ggplot2::coord_cartesian(xlim = c(0, max(fsc_max, ssc_max) * 1.05), ylim = c(0, max(fsc_max, ssc_max) * 1.05))
        ggplot2::ggsave(file.path(out_path, "fsc_ssc", paste0(sn, "_fsc_ssc.png")), p1, width = 5, height = 5, dpi = 300)

        p2 <- ggplot2::ggplot(data.table::data.table(x = vals_log), ggplot2::aes(x)) +
            ggplot2::geom_density(fill = "grey80", color = "grey40") +
            ggplot2::geom_vline(xintercept = log10(gate_min), color = "red", linewidth = 1) +
            ggplot2::geom_vline(xintercept = log10(gate_max), color = "red", linewidth = 1) +
            ggplot2::annotate("rect", xmin = log10(gate_min), xmax = log10(gate_max), ymin = -Inf, ymax = Inf, alpha = 0.15, fill = "red") +
            ggplot2::labs(title = paste0(sn, " - ", peak_channel), subtitle = paste0(round(100 * nrow(final_gated_data) / nrow(gated_data), 1), "% gated"), x = paste0("log10(", peak_channel, ")")) +
            ggplot2::theme_minimal() +
            ggplot2::theme(legend.position = "none")
        ggplot2::ggsave(file.path(out_path, "histogram", paste0(sn, "_histogram.png")), p2, width = 6.5, height = 4, dpi = 300)

        # Spectrum plot logic - based on working code from autogate_contour.R
        log_mat <- log10(pmax(final_gated_data[, detector_names, drop = FALSE], 1e-3))
        min_y <- floor(min(log_mat))
        max_y <- ceiling(max(log_mat))
        breaks <- seq(min_y, max_y, length.out = 151)
        bin_mid <- (breaks[-1] + breaks[-length(breaks)]) / 2
        bin_height <- breaks[2] - breaks[1]
        counts_mat <- vapply(seq_len(ncol(log_mat)), function(j) as.numeric(hist(log_mat[, j], breaks = breaks, plot = FALSE)$counts), numeric(length(bin_mid)))
        rownames(counts_mat) <- as.character(seq_along(bin_mid))
        colnames(counts_mat) <- as.character(seq_len(ncol(log_mat)))
        dt_c <- data.table::as.data.table(as.table(counts_mat))
        data.table::setnames(dt_c, c("bin_idx", "ch_idx", "count"))
        dt_c[, `:=`(bin_idx = as.integer(as.character(bin_idx)), ch_idx = as.integer(as.character(ch_idx)))]
        dt_c[, y_orig := bin_mid[bin_idx]]
        dt_c[, fill := log10(count + 1)]
        min_bin_count <- 3
        dt_c <- dt_c[count >= min_bin_count]
        y_power <- 1.5
        dt_c[, y := y_orig^y_power]

        if (nrow(dt_c) == 0) {
            message("  Warning: dt_c is empty for ", sn, ". Max count: ", max(counts_mat))
        } else {
            message("  Plotting ", nrow(dt_c), " bins for ", sn)
        }

        vlines <- which(diff(det_info$laser_nm) != 0) + 0.5
        fill_lo <- min(dt_c$fill, na.rm = TRUE)
        fill_hi <- quantile(dt_c$fill, 0.96, na.rm = TRUE)
        y_breaks_orig <- 0:ceiling(max_y)
        y_breaks_trans <- y_breaks_orig^y_power
        y_labels <- sapply(y_breaks_orig, function(x) bquote(10^.(x)))

        p3 <- ggplot2::ggplot(dt_c, ggplot2::aes(ch_idx, y, fill = fill)) +
            ggplot2::geom_tile(width = 0.7, height = bin_height * 3) +
            ggplot2::scale_fill_gradientn(
                colors = c("#0000FF", "#00FFFF", "#00FF00", "#FFFF00", "#FF0000"),
                limits = c(fill_lo, fill_hi), oob = scales::squish
            ) +
            ggplot2::scale_x_continuous(breaks = seq_along(detector_names), labels = detector_labels) +
            ggplot2::scale_y_continuous(limits = c(0, (max_y + 0.5)^y_power), breaks = y_breaks_trans, labels = y_labels) +
            ggplot2::coord_cartesian(expand = FALSE) +
            ggplot2::labs(title = paste0(sn, " - Spectrum"), x = NULL, y = "Intensity") +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6), legend.position = "none", panel.background = ggplot2::element_rect(fill = "white", color = NA))
        ggplot2::ggsave(file.path(out_path, "spectrum", paste0(sn, "_spectrum.png")), p3, width = 300, height = 120, units = "mm", dpi = 600)

        results_list[[sn]] <- data.table::data.table(sample = sn, fluorophore = fluor_name, type = sample_info$type, n_total = nrow(raw_data), n_final = nrow(final_gated_data), spectrum = list(spectrum_norm))
    }
    results_dt <- data.table::rbindlist(results_list)
    spectra_list <- results_dt$spectrum
    names(spectra_list) <- results_dt$fluorophore
    if (!is.null(af_data_raw)) spectra_list[["AF"]] <- af_data_raw / max(af_data_raw)
    M <- do.call(rbind, spectra_list)
    colnames(M) <- detector_names
    M_df <- as.data.frame(M)
    M_df$file <- rownames(M_df)
    data.table::fwrite(M_df[, c("file", colnames(M))], "reference_matrix.csv")
    message("Reference matrix (", nrow(M), " markers) saved.")
    return(M)
}
