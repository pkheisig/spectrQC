library(flowCore)
library(ggplot2)
library(devtools)
library(data.table)

rm(list = c("build_reference_matrix", "calc_residuals", "plot_scatter_rmse", "plot_marker_correlations", "plot_spectra"))
devtools::load_all()

# Set gating options ONCE - used by build_reference_matrix
opts <- gating_options(
  histogram_pct_beads = 0.98,
  histogram_direction_beads = "both",
  histogram_pct_cells = 0.35,
  histogram_direction_cells = "both"
)

# 1. Build reference matrix M
M <- build_reference_matrix(
  input_folder = "scc",
  output_folder = "autogate_plots",
  custom_fluorophores = c("US_UT_0 Run 1 20251211152908_aligned" = "Unstained"),
  gating_opts = opts
)

# 2. Plot spectra overlay
plot_spectra(ref_matrix = M)

# 4. Calc residuals on raw data
files_raw <- list.files("raw", pattern = ".*.fcs$", full.names = TRUE)
raw_data <- lapply(files_raw, read.FCS)
file_names <- tools::file_path_sans_ext(basename(files_raw))

# Calculate residuals using WLS (Global)
results <- mapply(function(ff, name) calc_residuals(ff, M, file_name = name, method = "WLS"),
  raw_data, file_names,
  SIMPLIFY = FALSE
)
results_df <- do.call(rbind, results)

# Plot using Relative RMSE (normalized by total intensity)
plot_scatter_rmse(results_df, metric = "Relative_RMSE", output_file = "scatter_rel_rmse.png")
plot_marker_correlations(results_df, metric = "Relative_RMSE", output_file = "marker_correlations_rel.png")

# Optional: Try NNLS or OLS
# results_ols <- mapply(function(ff, name) calc_residuals(ff, M, file_name = name, method = "OLS"),
#   raw_data, file_names, SIMPLIFY = FALSE)

# Optional: Use AutoSpectral integration
# create_autospectral_control_file(input_folder = "scc")
# unmixed_auto <- unmix_with_autospectral(raw_data[[1]], method = "WLS")

