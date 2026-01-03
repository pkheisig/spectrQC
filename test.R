# spectrQC Test Script
# Complete workflow demonstration with all functions and their parameters

devtools::load_all()

# =============================================================================
# PHASE 1: CHECK SINGLE COLOR CONTROLS (SCCs)
# =============================================================================
# Performs initial gating and extracts spectral signatures from SCC FCS files.
# Saves gating plots and spectrum overlays to output_dir.

# Load control file (optional but recommended for accurate fluorophore/channel mapping)
control_df <- data.table::fread("fcs_control_file.csv")

# Build reference matrix from SCC FCS files
M_initial <- build_reference_matrix(
    input_folder = "scc", # Directory containing SCC FCS files
    output_folder = "gating_and_spectrum_plots", # Output folder for gating/spectrum plots
    control_df = control_df, # Optional: data.table with filename, fluorophore, channel columns
    include_multi_af = FALSE, # Include extra AF signatures from af_dir
    af_dir = "af", # Directory for extra AF files
    default_sample_type = "beads", # Default sample type: "beads" or "cells"
    histogram_pct_beads = 0.98, # Histogram gate width for beads
    histogram_direction_beads = "both", # Histogram direction: "left", "right", or "both"
    histogram_pct_cells = 0.35, # Histogram gate width for cells
    histogram_direction_cells = "both", # Histogram direction for cells
    outlier_percentile = 0.02, # Percentile for FSC/SSC outlier removal
    debris_percentile = 0.02, # Percentile for debris gate
    bead_gate_scale = 1.3, # Scale factor for bead gate ellipse
    histogram_min_x_log = 2, # Minimum log10 value for histogram gating
    max_clusters = 6, # Maximum GMM clusters
    min_cluster_proportion = 0.03, # Minimum cluster proportion to keep
    gate_contour_beads = 0.999999999, # Contour level for bead gates
    gate_contour_cells = 0.95, # Contour level for cell gates
    subsample_n = 5000 # Number of events to subsample for GMM
)

# Generate the SCC-only QC report (PDF + individual PNGs)
generate_scc_report(
    M = M_initial, # Reference matrix from inspect_scc_spectra
    scc_dir = "scc", # SCC directory
    output_file = "SCC_QC_Report.pdf", # PDF output path
    custom_fluorophores = NULL, # Named vector mapping filenames to fluorophores
    png_dir = "spectrQC_outputs/plots/scc_audit", # Individual PNG plots folder
    control_file = "fcs_control_file.csv" # Control mapping file
)

# =============================================================================
# PHASE 2: REFINE SCCs (CLEANING)
# =============================================================================
# Cleans SCCs using RRMSE threshold to remove debris and outliers.
# Produces refined M and W matrices. Saves comparison plots.

refined <- refine_scc_matrix(
    M = M_initial, # Initial reference matrix
    scc_dir = "scc", # SCC folder
    rrmse_threshold = 0.05, # RRMSE cutoff (5%)
    output_dir = "scc_unmixed", # Folder for refined matrices and unmixed files
    control_file = "fcs_control_file.csv", # Control mapping file
    report_file = "SCC_Refinement_Report.pdf", # PDF comparison report
    png_dir = "spectrQC_outputs/plots/scc_refinement" # Individual comparison PNGs
)

M_final <- refined$M
W_final <- refined$W

# =============================================================================
# PHASE 3: UNMIX EXPERIMENTAL SAMPLES
# =============================================================================
# Applies refined signatures to experimental FCS data.
# Saves unmixed FCS files to output_dir.

unmixed_list <- unmix_samples(
    sample_dir = "samples", # Folder with experimental FCS files
    M = M_final, # Refined reference matrix
    method = "WLS", # Unmixing method: "OLS", "WLS", "NNLS", or "AutoSpectral"
    output_dir = "samples_unmixed" # Folder to save unmixed FCS files
)

# =============================================================================
# PHASE 4: QC EXPERIMENTAL SAMPLES
# =============================================================================
# Generates RRMSE scatter plots and diagnostics for unmixed samples.

generate_sample_qc(
    unmixed_list = unmixed_list, # List from unmix_samples
    M = M_final, # Reference matrix
    report_file = "Experimental_Sample_Audit.pdf", # PDF output
    png_dir = "spectrQC_outputs/plots/sample_audit", # Individual PNG plots
    pd = NULL # Optional pData for detector labels
)

# =============================================================================
# OPTIONAL: FULL COMBINED REPORT
# =============================================================================
# Generates a comprehensive QC report with all diagnostics.

results_df <- data.table::rbindlist(lapply(unmixed_list, `[[`, "data"))
generate_qc_report(
    results_df = results_df, # Combined unmixed data
    M = M_final, # Reference matrix
    output_file = "spectrQC_Report.pdf", # PDF output path
    res_list = unmixed_list[[1]], # Residuals list for detector diagnostics
    png_dir = "spectrQC_outputs/plots/report_pages", # Individual report pages
    pd = NULL # Optional pData for labels
)

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

# Calculate residuals for a single flowFrame
ff <- flowCore::read.FCS("samples/example.fcs")
res <- calc_residuals(
    flow_frame = ff, # flowFrame object with raw data
    M = M_final, # Reference matrix (Markers x Detectors)
    file_name = "example", # Optional file name for output
    method = "OLS", # Unmixing method: "OLS", "WLS", or "NNLS"
    background_noise = 25, # Baseline electronic noise for WLS
    return_residuals = FALSE # If TRUE, returns list with residuals matrix
)

# Derive unmixing matrix from reference matrix
W <- derive_unmixing_matrix(
    M = M_final, # Reference matrix (Markers x Detectors)
    method = "OLS", # "OLS" or "WLS"
    global_weights = NULL # Optional weights for WLS (one per detector)
)

# Save unmixing matrix to CSV
save_unmixing_matrix(
    W = W, # Unmixing matrix
    file = "unmixing_matrix.csv" # Output path
)

# Calculate spectral spread matrix (SSM)
ssm <- calculate_ssm(
    M = M_final, # Reference matrix
    method = "OLS", # Unmixing method
    background_noise = 100 # Background noise for WLS
)

# Get fluorophore patterns for auto-detection
patterns <- get_fluorophore_patterns() # Returns list with $unstained, $beads, $cells

# Get sorted detectors from FCS metadata
ff <- flowCore::read.FCS("scc/example.fcs")
pd <- flowCore::pData(flowCore::parameters(ff))
det_info <- get_sorted_detectors(pd) # Returns $names, $labels, $laser_nm

# =============================================================================
# PLOTTING FUNCTIONS
# =============================================================================

# Plot spectral overlays
plot_spectra(
    ref_matrix = M_final, # Reference matrix (Markers x Detectors)
    pd = NULL, # Optional pData for detector labels
    output_file = "spectra_overlay.png", # Output path
    width = 250, # Width in units
    height = 100, # Height in units
    unit = "mm", # Unit: "mm", "in", "cm"
    dpi = 300 # Resolution
)

# Plot RRMSE scatter (FSC vs SSC colored by error)
plot_scatter_rmse(
    data = results_df, # Data frame with FSC-A, SSC-A, and metric
    metric = "Relative_RMSE", # "RMSE_Score" or "Relative_RMSE"
    output_file = "scatter_rmse.png", # Output path
    width = 180, # Width
    height = 150, # Height
    unit = "mm", # Unit
    dpi = 300, # Resolution
    color_limits = c(0, 5) # Color scale limits
)

# Plot marker-RRMSE correlations
plot_marker_correlations(
    data = results_df, # Unmixed data frame
    metric = "Relative_RMSE", # "RMSE_Score" or "Relative_RMSE"
    markers = NULL, # Vector of marker names (NULL = all)
    output_file = "marker_correlations.png", # Output path
    show_smooth = TRUE, # Add GAM trend line
    y_limits = c(0, 10), # Optional y-axis limits
    width = 250, # Width
    height = 180, # Height
    unit = "mm", # Unit
    dpi = 300 # Resolution
)

# Plot spectral spread matrix
plot_ssm(
    SSM = ssm, # Matrix from calculate_ssm
    output_file = "spectral_spread_matrix.png", # Output path
    width = 200, # Width in mm
    height = 180 # Height in mm
)

# Plot unmixing matrix heatmap
plot_unmixing_matrix(
    W = W, # Unmixing matrix
    pd = NULL # Optional pData for labels
)

# =============================================================================
# INTERACTIVE GUI
# =============================================================================
# Launches the interactive web-based matrix adjustment interface.
# Requires npm to be installed for the frontend.

launch_gui(
    matrix_dir = getwd(), # Directory containing matrix CSVs
    samples_dir = NULL, # FCS samples directory (default: matrix_dir/samples)
    port = 8000, # API port
    open_browser = TRUE # Open browser automatically
)
