# specRQC: Spectral Flow Cytometry Quality Control

`specRQC` is an R package designed to validate spectral unmixing accuracy by calculating and visualizing unmixing residuals. It provides an automated pipeline to isolate single-stained control spectra and evaluate experimental data using photon-calibrated unmixing models.

## Key Features
- **Automated Gating**: Isolate positive populations from beads or cells using Gaussian Mixture Models (GMM).
- **Per-Cell WLS Unmixing**: Perform unmixing using the standard photon-counting variance model ($\sigma^2 = \text{Signal} + \text{Background}$).
- **Relative RMSE (RRMSE)**: A normalized error metric expressed as a percentage of total cell intensity.
- **AutoSpectral Integration**: Seamlessly generate control files and run unmixing through the `AutoSpectral` package.

---

## Installation

```r
# Install from GitHub
devtools::install_github("pkheisig/specRQC")
```

---

## The Workflow

The `specRQC` workflow is divided into three stages to allow for manual inspection of fluorophore mappings.

### 1. Setup Control File
Scans your single-stained control folder and guesses the fluorophore names and control types.

```r
library(specRQC)

# Generate a draft control file
create_autospectral_control_file(input_folder = "scc", output_file = "fcs_control_file.csv")

# NOTE: You can now open 'fcs_control_file.csv' in Excel/Text editor, 
# correct any fluorophore names, and save it.
```

### 2. Build Reference Matrix
Uses the control file to gate the positive populations and extract the normalized spectral signatures.

```r
# Load the mapping
control_df <- data.table::fread("fcs_control_file.csv")
custom_map <- setNames(control_df$fluorophore, tools::file_path_sans_ext(control_df$filename))

# Build matrix M (Markers x Detectors)
M <- build_reference_matrix(
  input_folder = "scc",
  output_folder = "autogate_plots",
  custom_fluorophores = custom_map
)

# Visualize spectra
plot_spectra(M)
```

### 3. Calculate Residuals & QC
Applies per-cell Weighted Least Squares (WLS) unmixing and calculates the Relative RMSE.

```r
library(flowCore)
ff <- read.FCS("raw_data.fcs")

# Calculate unmixed abundances and error metrics
results <- calc_residuals(ff, M, method = "WLS", background_noise = 100)

# Generate QC plots
plot_scatter_rmse(results, metric = "Relative_RMSE", output_file = "qc_scatter.png")
```

---

## Understanding RRMSE

**Relative RMSE (RRMSE)** is the primary quality metric in `specRQC`. It represents the "Spectral Fit Error" as a percentage of the total light collected from the cell.

$$ 	ext{RRMSE} = \frac{\sqrt{\text{mean}(\text{Residuals}^2)}}{\text{Total Intensity}} \times 100 $$

### Interpretation Guide:
- **< 1% (Excellent)**: The unmixing model perfectly explains the data. Standard for lymphocytes in well-calibrated panels.
- **1% – 3% (Good)**: High unmixing accuracy. Slight deviations may occur in very bright channels or due to autofluorescence minor mismatches.
- **3% – 10% (Warning)**: Significant residual error. This often indicates:
    - **Spectral Bottlenecks**: Two markers are too similar to resolve perfectly.
    - **Spillover Issues**: A marker is present that isn't in your reference matrix.
    - **Artifacts**: Debris, doublets, or extreme detector saturation.
- **> 10% (Poor)**: The unmixing for these cells is unreliable. These events should likely be gated out as debris or artifacts.

---

## Deriving a Static Unmixing Matrix

If you need a static unmixing matrix for use in other software (e.g., FlowJo), you can derive it from your reference matrix:

```r
# Get the 10x51 unmixing matrix (W)
W <- derive_unmixing_matrix(M, method = "OLS")

# Save to CSV
save_unmixing_matrix(W, "unmixing_matrix.csv")

# Manual unmixing math:
# Unmixed_Data = Raw_Data %*% t(W)
```

---

## Dependencies
`flowCore`, `ggplot2`, `dplyr`, `tidyr`, `data.table`, `mclust`, `sp`, `scales`, `nnls`.
`AutoSpectral` (Suggested for advanced unmixing).

---
**Author**: Paul Heisig  
**Email**: p.k.s.heisig@amsterdamumc.nl

```