library(flowCore)
library(ggplot2)
library(devtools)
library(data.table)

# Load the package
message("Loading package...")
devtools::load_all()

# --- STAGE 1: Create Control File ---
# This allows the user to inspect and edit fluorophore mappings
if (!file.exists("fcs_control_file.csv")) {
    message("Stage 1: Creating AutoSpectral control file...")
    control_df <- create_autospectral_control_file(
        input_folder = "scc",
        output_file = "fcs_control_file.csv"
    )
    # The user can now stop here and edit the CSV if needed.
    # For this autonomous run, we proceed.
} else {
    message("Stage 1: Using existing fcs_control_file.csv")
    control_df <- fread("fcs_control_file.csv")
}

# --- STAGE 2: Build Reference Matrix ---
message("Stage 2: Building reference matrix...")

# Map custom names from the control file for build_reference_matrix
# We convert the CSV into the named vector expected by build_reference_matrix
custom_map <- setNames(control_df$fluorophore, tools::file_path_sans_ext(control_df$filename))

opts <- gating_options(
  histogram_pct_beads = 0.98,
  histogram_direction_beads = "both",
  histogram_pct_cells = 0.35,
  histogram_direction_cells = "both"
)

M <- build_reference_matrix(
  input_folder = "scc",
  output_folder = "autogate_plots",
  custom_fluorophores = custom_map,
  gating_opts = opts
)

# Plot spectra overlay
plot_spectra(ref_matrix = M)

# --- STAGE 3: Unmixing & Quality Check ---
message("Stage 3: Unmixing raw data...")
files_raw <- list.files("raw", pattern = ".*.fcs$", full.names = TRUE)
if (length(files_raw) == 0) stop("No FCS files found in 'raw' folder")

results <- lapply(files_raw, function(f) {
  name <- tools::file_path_sans_ext(basename(f))
  message("  Processing: ", name)
  ff <- read.FCS(f, transformation = FALSE, truncate_max_range = FALSE)
  
  # Subsample for faster iteration and plot clarity
  n_cells <- nrow(ff)
  if (n_cells > 10000) {
      set.seed(42)
      ff <- ff[sample(n_cells, 10000), ]
  }
  
  # Use per-cell WLS with stable noise floor
  calc_residuals(ff, M, file_name = name, method = "WLS", background_noise = 100)
})

results_df <- rbindlist(results)
fwrite(results_df, "results_unmixed.csv")

message("Generating quality plots...")
# Plot Relative RMSE with the new percentage scaling (0-5%)
plot_scatter_rmse(results_df, 
                  metric = "Relative_RMSE", 
                  output_file = "scatter_rel_rmse.png",
                  color_limits = c(0, 5)) # Show 0% to 5% RRMSE

# Plot marker correlations with RRMSE percentage (0-10% range)
plot_marker_correlations(results_df, 
                         metric = "Relative_RMSE", 
                         output_file = "marker_correlations_rel.png",
                         y_limits = c(0, 10))

# Summary stats
summary_stats <- results_df[, .(
    Mean_RRMSE_Pct = mean(Relative_RMSE) * 100,
    Median_RRMSE_Pct = median(Relative_RMSE) * 100,
    Max_RRMSE_Pct = max(Relative_RMSE) * 100,
    Cells_Above_5Pct = sum(Relative_RMSE > 0.05)
), by = File]

fwrite(summary_stats, "summary_stats.csv")
print(summary_stats)

message("Analysis complete. Check 'scatter_rel_rmse.png' for visual QC.")