#!/usr/bin/env Rscript

# =============================================================================
#
#   AUTOMATED SPECTRAL FLOW CYTOMETRY GATING
#   ─────────────────────────────────────────
#   
#   This script automatically gates FCS files in two steps:
#     1) FSC/SSC gate to select your population of interest
#     2) Histogram gate on the brightest channel to select positive events
#   
#   HOW TO USE:
#     1. Place your FCS files in the "scc" folder
#     2. Adjust the settings below (START HERE section)
#     3. Run the script: Rscript autogate_contour.R
#     4. Find your plots in the output folder
#
# =============================================================================


# *****************************************************************************
# *                                                                           *
# *   ▼▼▼  START HERE - SETTINGS YOU'LL WANT TO CHANGE  ▼▼▼                   *
# *                                                                           *
# *****************************************************************************

# ─────────────────────────────────────────────────────────────────────────────
# INPUT / OUTPUT FOLDERS
# ─────────────────────────────────────────────────────────────────────────────
# Where are your FCS files? (relative to this script's location)
input_folder <- "scc"

# Where should output plots be saved?
output_folder <- "autogate_plots"


# ─────────────────────────────────────────────────────────────────────────────
# SAMPLE TYPE DETECTION
# ─────────────────────────────────────────────────────────────────────────────
# The script needs to know what type each sample is (beads, cells, or unstained)
# because they require different gating strategies.
#
# Add keywords from your filenames to each category below.
# Example: If your bead files contain "FITC" in the name, add "FITC" to beads.
#
# TIP: The script checks these patterns in order - put more specific patterns first!

sample_patterns <- list(
  
  # UNSTAINED CONTROLS (autofluorescence reference)
  unstained = c("US_UT", "Unstained", "unstained", "blank", "Blank", "AF only"),
  
  # SINGLE-STAINED BEADS (compensation controls)
  # This list includes ~60 common fluorophores - add more as needed!
  beads = c(
    # Alexa Fluor dyes
    "Alexa 350", "Alexa 405", "Alexa 430", "Alexa 488", "Alexa 514", 
    "Alexa 532", "Alexa 546", "Alexa 555", "Alexa 568", "Alexa 594", 
    "Alexa 610", "Alexa 633", "Alexa 647", "Alexa 660", "Alexa 680", 
    "Alexa 700", "Alexa 750", "Alexa 790",
    # FITC and basic dyes
    "FITC", "GFP", "YFP", "CFP", "mCherry", "DsRed", "tdTomato",
    # PE and PE conjugates
    "PE", "PE-Cy5", "PE-Cy5.5", "PE-Cy7", "PE-CF594", "PE-Texas Red", 
    "PE-Dazzle 594", "PE-Fire 640", "PE-Fire 700", "PE-Fire 810",
    "PE-Vio 615", "PE-Vio 770",
    # APC and APC conjugates
    "APC", "APC-Cy7", "APC-H7", "APC-R700", "APC-Fire 750", "APC-Fire 810",
    "APC-eFluor 780", "APC-Vio 770",
    # PerCP conjugates
    "PerCP", "PerCP-Cy5.5", "PerCP-eFluor 710", "PerCP-Vio 700",
    # BV (Brilliant Violet) dyes
    "BV421", "BV480", "BV510", "BV570", "BV605", "BV650", 
    "BV711", "BV750", "BV785",
    # BUV (Brilliant UV) dyes
    "BUV395", "BUV496", "BUV563", "BUV615", "BUV661", "BUV737", "BUV805",
    # Pacific dyes
    "Pacific Blue", "Pacific Orange",
    # Cy dyes
    "Cy3", "Cy5", "Cy5.5", "Cy7",
    # eFluor dyes (for beads)
    "eFluor 450", "eFluor 660",
    # Other common dyes
    "V450", "V500", "V500-C", "DAPI", "Hoechst", "PI", "7-AAD",
    "Indo-1", "Fluo-3", "Fluo-4", "Calcium Green",
    # SuperBright dyes
    "SB436", "SB600", "SB645", "SB702", "SB780",
    # Spark dyes
    "Spark Blue 550", "Spark NIR 685", "Spark YG 570", "Spark YG 581"
  ),
  
  # SINGLE-STAINED CELLS (e.g., viability markers)
  cells = c(
    "eFluor 506", "eFluor 520", "eFluor 780", "eFluor 455UV",
    "Zombie Aqua", "Zombie NIR", "Zombie UV", "Zombie Violet", 
    "Zombie Green", "Zombie Red", "Zombie Yellow",
    "LIVE/DEAD", "Fixable Viability", "FVD", "FVS",
    "Calcein", "SYTOX", "Viability"
  )
)

# If a file doesn't match any pattern above, treat it as:
# Options: "beads", "cells", or "unstained"
default_sample_type <- "beads"


# ─────────────────────────────────────────────────────────────────────────────
# HISTOGRAM GATING (Second Gate)
# ─────────────────────────────────────────────────────────────────────────────
# After FSC/SSC gating, events are gated on their brightest fluorescence channel.
# These settings control how much of the fluorescent peak to capture.

# ┌─────────────────────────────────────────────────────────────────────────┐
# │ FOR BEADS:                                                              │
# └─────────────────────────────────────────────────────────────────────────┘

# What percentage of events to capture? (0.0 to 1.0)
# Examples: 0.9 = 90%, 0.5 = 50%, 1.0 = 100%
histogram_pct_beads <- 0.5

# Which side of the peak to capture? (always starts from the CENTER of the peak)
#   "right" = capture brighter events (from center going right)
#   "left"  = capture dimmer events (from center going left)
#   "both"  = capture both sides equally (symmetric around center)
#
# Visual guide (| marks the center/median):
#                         CENTER
#                           │
#         ┌─────────────────│─────────────────┐
#   "left" │ ░░░░░░░████████│░░░░░░░░░░░░░░░░ │  ← X% leftward from center
#  "right" │ ░░░░░░░░░░░░░░░│████████░░░░░░░░ │  ← X% rightward from center  
#   "both" │ ░░░░░░░░░░█████│█████░░░░░░░░░░░ │  ← X/2% each side of center
#         └─────────────────│─────────────────┘
#
# Example with pct = 0.4:
#   "right" → captures events from 50th to 90th percentile (40% rightward)
#   "left"  → captures events from 10th to 50th percentile (40% leftward)
#   "both"  → captures events from 30th to 70th percentile (20% each side)
histogram_direction_beads <- "right"

# ┌─────────────────────────────────────────────────────────────────────────┐
# │ FOR CELLS (including unstained):                                        │
# └─────────────────────────────────────────────────────────────────────────┘

histogram_pct_cells <- 0.3
histogram_direction_cells <- "left"

# ┌─────────────────────────────────────────────────────────────────────────┐
# │ PER-FLUOROPHORE OVERRIDES (optional):                                   │
# └─────────────────────────────────────────────────────────────────────────┘
# Override the default bead/cell settings for specific fluorophores.
# These take priority over the general settings above.
#
# Format: "fluorophore_pattern" = c(pct, "direction")
# The pattern is matched against the FILENAME (case-insensitive).
#
# Examples:
#   "PE" = c(0.3, "left")        → PE files: 30% leftward from center
#   "PerCP-Cy5.5" = c(0.4, "right") → PerCP-Cy5.5: 40% rightward
#   "Alexa 647" = c(0.5, "both")    → Alexa 647: 50% symmetric

histogram_overrides <- list(
  # Add your overrides here, e.g.:
  "PE" = c(0.5, "left")
  # "PerCP-Cy5.5" = c(0.4, "right")
)


# *****************************************************************************
# *                                                                           *
# *   ▲▲▲  END OF COMMON SETTINGS  ▲▲▲                                        *
# *                                                                           *
# *   Most users can stop here! Advanced settings are below.                  *
# *                                                                           *
# *****************************************************************************


# ─────────────────────────────────────────────────────────────────────────────
# ADVANCED SETTINGS (usually don't need to change)
# ─────────────────────────────────────────────────────────────────────────────

# FSC/SSC Gate Settings
# ---------------------
outlier_percentile <- 0.02       # Remove top X% of extreme FSC/SSC values (0.02 = 2%)
debris_percentile <- 0.02        # For cells: exclude bottom X% as debris
bead_gate_scale <- 1.3           # Expand bead gate by this factor (1.3 = 30% larger)

# Histogram Detection
# -------------------
histogram_min_x_log <- 2         # Ignore fluorescence below 10^2 when finding peaks

# Clustering Algorithm (GMM)
# --------------------------
max_clusters <- 6                # Maximum populations to detect
min_cluster_proportion <- 0.03   # Ignore populations smaller than 3%
gate_contour_beads <- 0.9999999999  # Contour level for beads (higher = wider)
gate_contour_cells <- 0.98          # Contour level for cells
subsample_n <- 5000              # Use max N events for faster fitting


# =============================================================================
# INTERNAL VARIABLES (do not modify below this line)
# =============================================================================

# Map user-friendly variable names to internal names
out_dir <- output_folder
histogram_peak_pct <- histogram_pct_beads
histogram_direction <- histogram_direction_beads
histogram_peak_pct_cells <- histogram_pct_cells

# =============================================================================
# LOAD PACKAGES
# =============================================================================

suppressPackageStartupMessages({
  library(flowCore)
  library(data.table)
  library(ggplot2)
  library(mclust)
  library(scales)
})

# =============================================================================
# INPUT VALIDATION
# =============================================================================

# Check histogram settings are valid
stopifnot(
  "histogram_pct_beads must be between 0 and 1" = 
    histogram_pct_beads >= 0 && histogram_pct_beads <= 1,
  "histogram_pct_cells must be between 0 and 1" = 
    histogram_pct_cells >= 0 && histogram_pct_cells <= 1,
  "histogram_direction_beads must be 'both', 'left', or 'right'" = 
    histogram_direction_beads %in% c("both", "left", "right"),
  "histogram_direction_cells must be 'both', 'left', or 'right'" = 
    histogram_direction_cells %in% c("both", "left", "right"),
  "default_sample_type must be 'beads', 'cells', or 'unstained'" =
    default_sample_type %in% c("beads", "cells", "unstained")
)

# Check input folder exists
dir <- getwd()
fcs_dir <- paste0(dir, "/", input_folder)
if (!dir.exists(fcs_dir)) {
  stop("\n",
       "ERROR: Input folder not found!\n",
       "       Expected: ", fcs_dir, "\n",
       "       Please create this folder and add your FCS files.\n")
}

# Check for FCS files
fcs_files <- list.files(fcs_dir, pattern = "\\.fcs$", full.names = TRUE)
if (length(fcs_files) == 0) {
  stop("\n",
       "ERROR: No FCS files found!\n",
       "       Looked in: ", fcs_dir, "\n",
       "       Please add .fcs files to this folder.\n")
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

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
  data.table(x = rot_x + mean[1], y = rot_y + mean[2])
}

fit_gmm_populations <- function(data, max_k = 5, min_prop = 0.05) {
  fit <- Mclust(data, G = 1:max_k, verbose = FALSE)
  if (is.null(fit)) return(NULL)
  
  n_clusters <- fit$G
  proportions <- fit$parameters$pro
  means <- fit$parameters$mean
  
  if (fit$modelName %in% c("EII", "VII")) {
    sigmas <- lapply(1:n_clusters, function(k) diag(2) * fit$parameters$variance$sigmasq[k])
  } else if (fit$modelName %in% c("EEI", "VEI", "EVI", "VVI")) {
    sigmas <- lapply(1:n_clusters, function(k) diag(fit$parameters$variance$sigma[, , k]))
  } else {
    sigmas <- lapply(1:n_clusters, function(k) fit$parameters$variance$sigma[, , k])
  }
  
  main_pops <- which(proportions >= min_prop)
  
  list(fit = fit, n_clusters = n_clusters, model_name = fit$modelName, BIC = fit$BIC,
       proportions = proportions, means = means, sigmas = sigmas,
       main_populations = main_pops, classification = fit$classification)
}

get_sample_type <- function(filename, patterns, default) {
  for (type in names(patterns)) {
    for (pattern in patterns[[type]]) {
      if (grepl(pattern, filename, ignore.case = TRUE)) {
        return(list(type = type, pattern = pattern))
      }
    }
  }
  list(type = default, pattern = "default")
}

select_bead_population <- function(gmm_result) {
  main_pops <- gmm_result$main_populations
  if (length(main_pops) == 0) return(NULL)
  proportions <- gmm_result$proportions[main_pops]
  best_pop <- main_pops[which.max(proportions)]
  list(selected = best_pop, reason = "largest population", 
       proportion = gmm_result$proportions[best_pop])
}

select_cell_populations <- function(gmm_result, debris_threshold, detection_limit) {
  main_pops <- gmm_result$main_populations
  if (length(main_pops) == 0) return(NULL)
  
  valid_pops <- c()
  for (k in main_pops) {
    mean_fsc <- gmm_result$means[1, k]
    mean_ssc <- gmm_result$means[2, k]
    sigma <- gmm_result$sigmas[[k]]
    eig <- eigen(sigma)$values
    elongation <- max(eig) / min(eig)
    
    # More aggressive exclusion: 80% of detection limit
    # Also check if ellipse would extend past limit
    ell_extent <- sqrt(max(eig) * qchisq(0.95, df = 2))
    max_extent_fsc <- mean_fsc + ell_extent
    max_extent_ssc <- mean_ssc + ell_extent
    
    if (mean_fsc >= debris_threshold && 
        mean_fsc <= detection_limit * 0.8 && 
        mean_ssc <= detection_limit * 0.8 && 
        max_extent_fsc <= detection_limit &&
        max_extent_ssc <= detection_limit &&
        elongation <= 50) {
      valid_pops <- c(valid_pops, k)
    }
  }
  list(selected = valid_pops, proportions = gmm_result$proportions[valid_pops])
}

create_merged_gate <- function(gmm_result, populations, level = 0.95, n_points = 100, 
                                clip_max_x = Inf, clip_max_y = Inf, scale = 1.0) {
  if (length(populations) == 0) return(NULL)
  
  if (length(populations) == 1) {
    ell <- get_ellipse(gmm_result$means[, populations], gmm_result$sigmas[[populations]], 
                       level, n_points, scale)
  } else {
    all_points <- rbindlist(lapply(populations, function(k) {
      get_ellipse(gmm_result$means[, k], gmm_result$sigmas[[k]], level, n_points, scale)
    }))
    hull_idx <- chull(all_points$x, all_points$y)
    ell <- all_points[c(hull_idx, hull_idx[1]), ]
  }
  
  # Clip to detection limits
  ell[, x := pmin(x, clip_max_x)]
  ell[, y := pmin(y, clip_max_y)]
  ell[, x := pmax(x, 0)]
  ell[, y := pmax(y, 0)]
  
  ell
}

# =============================================================================
# MAIN SCRIPT
# =============================================================================

out_path <- paste0(dir, "/", out_dir)

dir.create(paste0(out_path, "/fsc_ssc"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(out_path, "/histogram"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(out_path, "/spectrum"), showWarnings = FALSE, recursive = TRUE)

message("\n", paste(rep("=", 70), collapse = ""))
message("AUTOMATED GATING: FSC/SSC + Histogram + Spectrum")
message(paste(rep("=", 70), collapse = ""))
message("\nInput folder:  ", fcs_dir)
message("Output folder: ", out_path)
message("Files found:   ", length(fcs_files), " FCS files\n")

results_list <- list()

for (fcs_file in fcs_files) {
  sn <- tools::file_path_sans_ext(basename(fcs_file))
  message("\n", paste(rep("-", 60), collapse = ""))
  message("Sample: ", sn)
  message(paste(rep("-", 60), collapse = ""))
  
  sample_info <- get_sample_type(sn, sample_patterns, default_sample_type)
  message("Type: ", toupper(sample_info$type), " (matched: ", sample_info$pattern, ")")
  
  ff <- read.FCS(fcs_file, transformation = FALSE, truncate_max_range = FALSE)
  pd <- pData(parameters(ff))
  
  fsc <- pd$name[grepl("^FSC", pd$name) & grepl("-A$", pd$name)][1]
  ssc <- pd$name[grepl("^SSC", pd$name) & grepl("-A$", pd$name)][1]
  fsc_desc <- pd$desc[pd$name == fsc]
  ssc_desc <- pd$desc[pd$name == ssc]
  
  fl_idx <- which(grepl("^FL[0-9]+-A$", pd$name))
  if (length(fl_idx) == 0) fl_idx <- which(grepl("^FL[0-9]+-Comp$", pd$name))
  fl_channels <- pd$name[fl_idx]
  
  data_raw <- exprs(ff)[, c(fsc, ssc)]
  n_total <- nrow(data_raw)
  message("Total events: ", format(n_total, big.mark = ","))
  
  # ==========================================================================
  # PRE-FILTERING
  # ==========================================================================
  
  fsc_max <- quantile(data_raw[, 1], 1 - outlier_percentile, na.rm = TRUE)
  ssc_max <- quantile(data_raw[, 2], 1 - outlier_percentile, na.rm = TRUE)
  plot_max <- max(fsc_max, ssc_max) * 1.05
  
  valid_idx <- which(data_raw[, 1] < fsc_max & data_raw[, 2] < ssc_max & 
                     data_raw[, 1] > 0 & data_raw[, 2] > 0)
  data_filtered <- data_raw[valid_idx, ]
  n_filtered <- nrow(data_filtered)
  message("After filtering top ", outlier_percentile * 100, "% outliers: ", format(n_filtered, big.mark = ","))
  
  if (sample_info$type %in% c("cells", "unstained")) {
    debris_threshold <- quantile(data_filtered[, 1], debris_percentile, na.rm = TRUE)
    debris_idx <- which(data_filtered[, 1] >= debris_threshold)
    data_filtered <- data_filtered[debris_idx, ]
    message("After removing debris: ", format(nrow(data_filtered), big.mark = ","))
  } else {
    debris_threshold <- 0
  }
  
  if (!is.null(subsample_n) && nrow(data_filtered) > subsample_n) {
    sample_idx <- sample(nrow(data_filtered), subsample_n)
    data_fit <- data_filtered[sample_idx, ]
  } else {
    data_fit <- data_filtered
  }
  
  # ==========================================================================
  # FIT GMM
  # ==========================================================================
  
  message("Fitting GMM...")
  gmm_result <- fit_gmm_populations(data_fit, max_k = max_clusters, min_prop = min_cluster_proportion)
  
  if (is.null(gmm_result)) {
    warning("GMM fitting failed for ", sn)
    next
  }
  
  message("Model: ", gmm_result$model_name, " with ", gmm_result$n_clusters, " components")
  
  # Select populations based on sample type
  if (sample_info$type == "beads") {
    selection <- select_bead_population(gmm_result)
    selected_pops <- selection$selected
    gate_level <- gate_contour_beads
    message("Selected: Pop ", selected_pops, " (", round(selection$proportion * 100, 1), "%)")
  } else {
    # Cells and unstained: merge populations, exclude outliers
    selection <- select_cell_populations(gmm_result, debris_threshold, min(fsc_max, ssc_max))
    selected_pops <- selection$selected
    gate_level <- gate_contour_cells
    if (length(selected_pops) == 0) {
      # Fallback: use largest population
      selection <- select_bead_population(gmm_result)
      selected_pops <- selection$selected
    }
    message("Selected: Pops ", paste(selected_pops, collapse = "+"), 
            " (", round(sum(gmm_result$proportions[selected_pops]) * 100, 1), "%)")
  }
  
  if (length(selected_pops) == 0) {
    warning("No valid populations for ", sn)
    next
  }
  
  # Generate gate polygon
  if (sample_info$type %in% c("cells", "unstained")) {
    # Cells/unstained: clip to detection limits
    final_gate <- create_merged_gate(gmm_result, selected_pops, gate_level, 
                                      clip_max_x = fsc_max * 0.95, clip_max_y = ssc_max * 0.95)
  } else {
    # Beads: use scale factor for larger gate
    final_gate <- create_merged_gate(gmm_result, selected_pops, gate_level, scale = bead_gate_scale)
  }
  
  # ==========================================================================
  # APPLY GATE
  # ==========================================================================
  
  all_data <- exprs(ff)
  in_gate <- sp::point.in.polygon(all_data[, fsc], all_data[, ssc], 
                                   final_gate$x, final_gate$y) > 0
  gated_data <- all_data[in_gate, ]
  n_gated <- nrow(gated_data)
  message("Events in FSC/SSC gate: ", format(n_gated, big.mark = ","), 
          " (", round(100 * n_gated / n_total, 1), "%)")
  
  # ==========================================================================
  # PLOT 1: FSC/SSC
  # ==========================================================================
  
  dt_plot <- data.table(x = data_raw[, 1], y = data_raw[, 2])
  x_lim <- c(0, plot_max)
  y_lim <- c(0, plot_max)
  
  x_breaks <- seq(0, plot_max, length.out = 201)
  y_breaks <- seq(0, plot_max, length.out = 201)
  dt_plot[, x_bin := findInterval(x, x_breaks, rightmost.closed = TRUE)]
  dt_plot[, y_bin := findInterval(y, y_breaks, rightmost.closed = TRUE)]
  dt2d <- dt_plot[x_bin >= 1 & x_bin <= 200 & y_bin >= 1 & y_bin <= 200, 
                   .(count = .N), by = .(x_bin, y_bin)]
  dt2d[, `:=`(x = x_breaks[x_bin], y = y_breaks[y_bin], fill = log10(count + 1))]
  
  p1 <- ggplot(dt2d, aes(x, y, fill = fill)) +
    geom_tile(width = diff(x_breaks)[1], height = diff(y_breaks)[1]) +
    scale_fill_gradientn(colors = c("#0000FF", "#00FFFF", "#00FF00", "#FFFF00", "#FF0000"), guide = "none") +
    geom_path(data = final_gate, aes(x, y), inherit.aes = FALSE, color = "red", linewidth = 1) +
    labs(title = paste0(sn, " - FSC/SSC"), subtitle = paste0(round(100 * n_gated / n_total, 1), "% gated"),
         x = fsc_desc, y = ssc_desc) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none", panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA)) +
    coord_cartesian(xlim = x_lim, ylim = y_lim)
  
  ggsave(paste0(out_path, "/fsc_ssc/", sn, "_fsc_ssc.png"), p1, width = 5, height = 5, dpi = 300)
  
  # ==========================================================================
  # HISTOGRAM GATING
  # ==========================================================================
  
  if (n_gated < 100 || length(fl_channels) == 0) {
    message("Skipping histogram/spectrum (too few events or no FL channels)")
    next
  }
  
  q99 <- apply(gated_data[, fl_channels, drop = FALSE], 2, quantile, probs = 0.99)
  peak_channel <- fl_channels[which.max(q99)]
  peak_desc <- pd$desc[pd$name == peak_channel]
  message("Peak channel: ", peak_channel, " (", peak_desc, ")")
  
  peak_vals <- gated_data[, peak_channel]
  vals_log <- log10(pmax(peak_vals, 1))
  
  # Find valley and peak - different logic for unstained
  if (sample_info$type == "unstained") {
    # Unstained: no minimum x, find the main peak directly
    d_full <- density(vals_log, n = 512)
    peak_x <- d_full$x[which.max(d_full$y)]
    vals_above <- vals_log  # Use all data
  } else {
    # Beads/cells: find valley above min_x
    d_full <- density(vals_log, n = 512)
    mid_idx <- which(d_full$x > histogram_min_x_log & d_full$x < quantile(vals_log, 0.98))
    if (length(mid_idx) > 0) {
      local_min_idx <- mid_idx[which.min(d_full$y[mid_idx])]
      min_x <- d_full$x[local_min_idx]
    } else {
      min_x <- histogram_min_x_log
    }
    vals_above <- vals_log[vals_log >= min_x]
    if (length(vals_above) > 10) {
      d <- density(vals_above, n = 512)
      peak_x <- d$x[which.max(d$y)]
    } else {
      peak_x <- median(vals_above)
    }
  }
  
  # Gate boundaries based on sample type and direction
  # First, check for fluorophore-specific overrides
  override_found <- FALSE
  if (length(histogram_overrides) > 0) {
    for (pattern in names(histogram_overrides)) {
      if (grepl(pattern, sn, ignore.case = TRUE)) {
        override <- histogram_overrides[[pattern]]
        pct <- override[1]
        direction <- override[2]
        override_found <- TRUE
        message("Using override for '", pattern, "': ", round(pct * 100), "% ", direction)
        break
      }
    }
  }
  
  # If no override found, use default settings based on sample type
  if (!override_found) {
    if (sample_info$type == "unstained") {
      pct <- histogram_peak_pct_cells
      direction <- histogram_direction_cells
    } else if (sample_info$type == "beads") {
      pct <- histogram_peak_pct
      direction <- histogram_direction
    } else {
      pct <- histogram_peak_pct_cells
      direction <- histogram_direction_cells
    }
  }
  
  # Calculate quantiles based on direction (always starting from center = 50th percentile)
  if (direction == "right") {
    # Capture from center going rightward (brighter events)
    # pct=0.4 → from 50th to 90th percentile
    lower_q <- 0.5
    upper_q <- 0.5 + pct
  } else if (direction == "left") {
    # Capture from center going leftward (dimmer events)
    # pct=0.4 → from 10th to 50th percentile
    lower_q <- 0.5 - pct
    upper_q <- 0.5
  } else {
    # "both" (default) → symmetric around center, half on each side
    # pct=0.4 → from 30th to 70th percentile (20% each side)
    lower_q <- 0.5 - pct / 2
    upper_q <- 0.5 + pct / 2
  }
  
  # Clamp to valid quantile range [0, 1]
  lower_q <- max(0, lower_q)
  upper_q <- min(1, upper_q)
  
  gate_min <- 10^quantile(vals_above, lower_q)
  gate_max <- 10^quantile(vals_above, upper_q)
  message("Histogram gate: ", direction, " direction, ", round(pct * 100), "% of events")
  gate_min_log <- log10(gate_min)
  gate_max_log <- log10(gate_max)
  
  in_hist_gate <- peak_vals >= gate_min & peak_vals <= gate_max
  n_hist_gated <- sum(in_hist_gate)
  final_gated_data <- gated_data[in_hist_gate, ]
  message("Events in histogram gate: ", format(n_hist_gated, big.mark = ","), 
          " (", round(100 * n_hist_gated / n_gated, 1), "% of FSC/SSC gated)")
  
  # ==========================================================================
  # PLOT 2: Histogram
  # ==========================================================================
  
  p2 <- ggplot(data.table(x = vals_log), aes(x)) +
    geom_density(fill = "grey80", color = "grey40") +
    geom_vline(xintercept = gate_min_log, color = "red", linewidth = 1) +
    geom_vline(xintercept = gate_max_log, color = "red", linewidth = 1) +
    annotate("rect", xmin = gate_min_log, xmax = gate_max_log, ymin = -Inf, ymax = Inf, 
             alpha = 0.15, fill = "red") +
    labs(title = paste0(sn, " - ", peak_desc),
         subtitle = paste0(round(100 * n_hist_gated / n_gated, 1), "% of FSC/SSC gated"),
         x = paste0("log10(", peak_desc, ")")) +
    theme_minimal() + theme(legend.position = "none")
  
  ggsave(paste0(out_path, "/histogram/", sn, "_histogram.png"), p2, width = 6.5, height = 4, dpi = 300)
  
  # ==========================================================================
  # PLOT 3: Spectrum (heatmap style - exactly like autogate_opencyto.R)
  # ==========================================================================
  
  if (nrow(final_gated_data) < 10) {
    message("Skipping spectrum (too few events after histogram gate)")
    next
  }
  
  # Order channels by laser and wavelength
  channels_desc <- pd$desc[match(fl_channels, pd$name)]
  channels_desc <- sub("-(A|H|W)$", "", sub("-Comp$", "", channels_desc))
  laser_nm <- as.integer(sub("^([0-9]+)nm.*$", "\\1", channels_desc))
  center_nm <- as.integer(sub("^([0-9]+).*$", "\\1", sub("^.*-\\s*", "", channels_desc)))
  ord <- order(laser_nm, center_nm)
  fl_channels_ord <- fl_channels[ord]
  channels_desc_ord <- channels_desc[ord]
  laser_nm_ord <- laser_nm[ord]
  
  # Extract FL data and compute log
  mat_fl <- final_gated_data[, fl_channels_ord, drop = FALSE]
  log_mat <- log10(pmax(mat_fl, 1e-3))
  min_y <- floor(min(log_mat))
  max_y <- ceiling(max(log_mat))
  breaks <- seq(min_y, max_y, length.out = 151)
  bin_mid <- (breaks[-1] + breaks[-length(breaks)]) / 2
  bin_height <- breaks[2] - breaks[1]
  
  counts_mat <- vapply(seq_len(ncol(log_mat)), function(j) 
    as.numeric(hist(log_mat[,j], breaks = breaks, plot = FALSE)$counts), numeric(length(bin_mid)))
  rownames(counts_mat) <- as.character(seq_along(bin_mid))
  colnames(counts_mat) <- as.character(seq_len(ncol(log_mat)))
  
  dt <- as.data.table(as.table(counts_mat))
  setnames(dt, c("bin_idx", "ch_idx", "count"))
  dt[, bin_idx := as.integer(as.character(bin_idx))]
  dt[, ch_idx := as.integer(as.character(ch_idx))]
  dt[, y_orig := bin_mid[bin_idx]]
  dt[, fill := log10(count + 1)]
  
  # Minimum events per bin to display (filters out artifact/noise)
  min_bin_count <- 3
  dt <- dt[count >= min_bin_count]
  
  # Power transformation to compress lower values (higher = more compression)
  y_power <- 1.5
  dt[, y := y_orig^y_power]
  
  fill_lo <- min(dt$fill)
  fill_hi <- quantile(dt$fill, 0.96)
  vlines <- which(diff(laser_nm_ord) != 0) + 0.5
  
  # Create breaks and labels for transformed y-axis
  y_breaks_orig <- 0:ceiling(max_y)
  y_breaks_trans <- y_breaks_orig^y_power
  y_labels <- sapply(y_breaks_orig, function(x) bquote(10^.(x)))
  
  p3 <- ggplot(dt, aes(ch_idx, y, fill = fill)) +
    geom_tile(width = 0.7, height = bin_height*3) +
    scale_fill_gradientn(colors = c("#0000FF", "#00FFFF", "#00FF00", "#FFFF00", "#FF0000"),
                         limits = c(fill_lo, fill_hi), oob = scales::squish) +
    { if (length(vlines) > 0) geom_vline(xintercept = vlines, color = "grey30", linewidth = 0.5) } +
    scale_x_continuous(breaks = seq_along(fl_channels_ord), labels = channels_desc_ord) +
    scale_y_continuous(limits = c(0, (max_y+0.5)^y_power), breaks = y_breaks_trans, labels = y_labels) +
    coord_cartesian(expand = FALSE) +
    labs(title = paste0(sn, " - Spectrum (", format(nrow(final_gated_data), big.mark = ","), " events)"), x = NULL, y = "Intensity") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.line = element_line(linewidth = 0.5),
      plot.title = element_text(size = 15),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.background = element_rect(fill = "white")
    )
  
  ggsave(paste0(out_path, "/spectrum/", sn, "_spectrum.png"), p3, 
         width = 300, height = 120, units = "mm", dpi = 600)
  
  # ==========================================================================
  # SAVE RESULTS
  # ==========================================================================
  
  results_list[[sn]] <- data.table(
    sample = sn,
    type = sample_info$type,
    matched_pattern = sample_info$pattern,
    n_total = n_total,
    n_fsc_ssc_gated = n_gated,
    fsc_ssc_pct = round(100 * n_gated / n_total, 2),
    peak_channel = peak_channel,
    peak_desc = peak_desc,
    n_histogram_gated = n_hist_gated,
    histogram_pct = round(100 * n_hist_gated / n_gated, 2),
    n_final = nrow(final_gated_data),
    final_pct = round(100 * nrow(final_gated_data) / n_total, 2)
  )
  
  message("Saved plots for: ", sn)
}

# Save summary
if (length(results_list) > 0) {
  results_dt <- rbindlist(results_list)
  fwrite(results_dt, paste0(out_path, "/gating_summary.csv"))
  message("\n", paste(rep("=", 70), collapse = ""))
  message("SUMMARY")
  message(paste(rep("=", 70), collapse = ""))
  print(results_dt[, .(sample, type, fsc_ssc_pct = paste0(fsc_ssc_pct, "%"), 
                       histogram_pct = paste0(histogram_pct, "%"),
                       final_pct = paste0(final_pct, "%"))])
}

message("\nDone! Plots saved to: ", out_path)
