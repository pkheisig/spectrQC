library(devtools)
load_all()
library(data.table)

# Test the cleaning logic
control_df <- fread("fcs_control_file.csv")
custom_map <- setNames(control_df$fluorophore, tools::file_path_sans_ext(control_df$filename))

message("Running reference matrix build WITH Spectral Cleaning (5% threshold)...")
M_clean <- build_reference_matrix(
  input_folder = "scc",
  output_folder = "cleaning_test_plots",
  custom_fluorophores = custom_map,
  clean_rrmse_threshold = 0.05
)

message("\nRunning reference matrix build WITHOUT Spectral Cleaning...")
M_raw <- build_reference_matrix(
  input_folder = "scc",
  output_folder = "raw_test_plots",
  custom_fluorophores = custom_map,
  clean_rrmse_threshold = 0 # Disable
)

# Compare difference in signatures
diff_pct <- max(abs(M_clean - M_raw)) * 100
message("\nMaximum spectral signature change: ", round(diff_pct, 4), "%")

if (diff_pct > 0) {
    message("Spectral cleaning successfully shifted the signatures by removing outliers.")
} else {
    message("No significant outliers found in these controls.")
}

