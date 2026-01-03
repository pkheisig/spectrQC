# spectrQC: Full Spectrum Flow Cytometry Quality Control

`spectrQC` is an R package for validating spectral unmixing accuracy. It provides a complete pipeline to extract fluorophore signatures from single-color controls, refine them using quality metrics, and unmix experimental samples with detailed diagnostic reports.

## Key Features

- **Automated Gating**: Isolate positive populations from beads or cells using Gaussian Mixture Models
- **Background Subtraction**: Automatically subtract internal negative populations to isolate pure fluorophore signatures
- **Signature Refinement**: Remove debris and outliers using RRMSE-based filtering
- **Per-Cell WLS Unmixing**: High-accuracy unmixing using photon-counting variance weighting
- **Comprehensive Reports**: PDF reports with spectral overlays, spread matrices, and quality diagnostics
- **Interactive GUI**: Web-based interface for manual matrix adjustment

---

## Installation

```r
# Install from GitHub
devtools::install_github("pkheisig/spectrQC")
```

---

## Project Setup

Organize your data with this folder structure:

```
my_project/
├── scc/                          # Single-color control FCS files
│   ├── FITC_beads.fcs
│   ├── PE_beads.fcs
│   └── Unstained.fcs             # Autofluorescence control
├── samples/                      # Experimental FCS files
│   └── Sample1.fcs
└── fcs_control_file.csv          # Control file mapping (optional)
```

### Control File

The control file maps FCS filenames to fluorophores. Generate it using the [AutoSpectral](https://github.com/DrCytometer/AutoSpectral) package:

```r
library(AutoSpectral)
asp <- get.autospectral.param()
create.control.file("scc", asp)
```

Or create manually with columns: `filename`, `fluorophore`, `channel`, `universal.negative`

---

## Workflow

### Step 1: Build Reference Matrix

Extract spectral signatures from single-color controls:

```r
library(spectrQC)

# Load control file (optional but recommended)
control_df <- data.table::fread("fcs_control_file.csv")

# Build reference matrix from SCC files
M <- build_reference_matrix(
  input_folder = "scc",
  output_folder = "gating_plots",
  control_df = control_df,
  default_sample_type = "beads"
)
```

This saves gating and spectrum plots to `gating_plots/` and exports `reference_matrix.csv`.

#### Gating Parameters

Adjust these if auto-gating fails:

```r
M <- build_reference_matrix(
  input_folder = "scc",
  output_folder = "gating_plots",
  histogram_pct_beads = 0.98,         # Width of positive gate for beads
  histogram_pct_cells = 0.35,         # Width of positive gate for cells
  max_clusters = 6,                   # Maximum GMM clusters
  gate_contour_beads = 0.999999999,   # Contour level for bead gates
  gate_contour_cells = 0.95           # Contour level for cell gates
)
```

---

### Step 2: Generate SCC Report

Review the extracted signatures:

```r
generate_scc_report(
  M = M,
  scc_dir = "scc",
  output_file = "SCC_QC_Report.pdf"
)
```

---

### Step 3: Refine the Matrix

Remove outliers using RRMSE-based filtering:

```r
refined <- refine_scc_matrix(
  M = M,
  scc_dir = "scc",
  rrmse_threshold = 0.05,             # Remove events with >5% error
  output_dir = "scc_unmixed"
)

M_final <- refined$M                  # Refined reference matrix
W_final <- refined$W                  # Unmixing matrix
```

This exports `refined_reference_matrix.csv` and `refined_unmixing_matrix.csv`.

---

### Step 4: Unmix Experimental Samples

Apply the refined matrix to your samples:

```r
unmixed <- unmix_samples(
  sample_dir = "samples",
  M = M_final,
  method = "WLS",                     # "OLS", "WLS", or "NNLS"
  output_dir = "samples_unmixed"
)
```

**Methods:**
- **OLS**: Ordinary least squares — fast, suitable for most panels
- **WLS**: Weighted least squares — accounts for photon-counting noise, best accuracy
- **NNLS**: Non-negative least squares — forces positive abundances

---

### Step 5: QC Experimental Samples

Generate the final quality report:

```r
generate_sample_qc(
  unmixed_list = unmixed,
  M = M_final,
  report_file = "Sample_QC_Report.pdf"
)
```

---

## Interactive Matrix Adjustment

For manual fine-tuning, launch the web interface:

```r
launch_gui(
  matrix_dir = getwd(),
  samples_dir = "samples"
)
```

Open http://localhost:5174 to access the GUI.

---

## Understanding the Reports

### Reference Spectra Overlay

Shows the spectral signature of each fluorophore across all detectors.

- **Expected**: Smooth curves with distinct peaks
- **Issues**: Jagged lines or unexpected peaks indicate gating failures or contamination

### Spectral Spread Matrix

Quantifies how much noise each fluorophore introduces into other channels after unmixing.

- **Expected**: Low values (dark colors) indicate minimal interference
- **Issues**: High values between two fluorophores indicate spectral similarity, reducing sensitivity for dim co-expressed populations

### RRMSE Scatter Plots

FSC vs SSC colored by Relative Root Mean Square Error, showing unmixing quality per cell.

- **Expected**: Most cells gray or light (RRMSE < 5%)
- **Issues**: Red clusters indicate populations where unmixing failed — often debris, autofluorescence mismatches, or missing signatures

### Detector Residuals

Boxplots showing the residual signal in each detector for high-error cells.

- **Expected**: Small distributions centered at zero
- **Issues**: Large positive or negative spikes indicate unaccounted signals (contaminants or missing fluorophores)

### Marker-RRMSE Correlations

Unmixed intensity vs RRMSE for each marker.

- **Expected**: Flat trend line
- **Issues**: Upward trend indicates incorrect signature or detector non-linearity

---

## Output Directories

- `gating_plots/`: Gating and spectrum visualizations for each SCC
- `scc_unmixed/`: Refined matrices and unmixed control files
- `samples_unmixed/`: Unmixed experimental data (FCS format)
- `spectrQC_outputs/plots/`: Individual PNG exports from reports

---

**Author**: Paul Heisig  
**Email**: p.k.s.heisig@amsterdamumc.nl
