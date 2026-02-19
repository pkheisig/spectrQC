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
asp <- get.autospectral.param(cytometer = "aurora")
create.control.file("scc", asp)
```

Or create manually with columns: `filename`, `fluorophore`, `channel` (`universal.negative` is optional).

For auto-generation, `spectrQC` uses shipped dictionaries:
- `inst/extdata/fluorophore_dictionary.csv`
- `inst/extdata/marker_dictionary.csv`

Logic:
- detect `fluorophore` and `marker` from filename aliases first
- if no fluorophore is detected, fallback to peak-channel mapping (for example `YG1 -> PE`, `UV2 -> BUV395`)
- if a file is unlabeled (no marker/fluor match) but matches cytometer AF-channel behavior, it is auto-tagged as `AF`
- if no marker is detected, marker is left empty

---

## Workflow

### Typical Sequence (Field-Friendly)

You usually follow one of these two paths:

#### Path A: Fast SCC Check (`quick_unmix`)

Use this when you want a quick SCC-only run (build matrix, SCC plots, SCC unmixing, unmixing matrix plot):

```r
library(spectrQC)

quick_unmix(
  scc_dir = "scc",
  control_file = "fcs_control_file.csv",
  auto_create_control = TRUE,
  cytometer = "Aurora",
  auto_unknown_fluor_policy = "by_channel",
  output_dir = "spectrQC_outputs/quick_unmix",
  unmix_method = "WLS",
  build_qc_plots = TRUE,
  unmix_scatter_panel_size_mm = 30
)
```

If `fcs_control_file.csv` is missing and `auto_create_control = TRUE`, `quick_unmix()` auto-generates a control file (filename, marker, fluorophore, and detected peak channel), then asks for confirmation before continuing.
`cytometer` is used for channel-aware fluorophore inference via the AutoSpectral fluorophore database.

For newly auto-created files:
- `control.type` is set to `cells` only for AF rows; all non-AF rows are left empty on purpose
- `universal.negative` is left empty by default for all rows
- `quick_unmix()` pauses and asks for `y/n` confirmation so you can review/edit the file first

`quick_unmix()` writes:
- `scc_spectra.png` (reference spectra overlay)
- `scc_unmixing_matrix.png` and `scc_unmixing_matrix.csv`
- `scc_unmixing_scatter_matrix.png` (lower-triangle scatter matrix, one single-stain file per row, with x=0/y=0 guides)

Set `unmix_scatter_panel_size_mm` higher (for example `40`) if you want larger per-panel scatter plots.

`quick_unmix()` also runs a strict preflight check before processing:
- every SCC file must be mapped in `fcs_control_file.csv`
- non-AF rows must define a valid `channel`
- if `universal.negative` is present, values for active SCC rows must be empty/keyword or reference a file present in your selected SCC/AF directories

If preflight fails, fix the listed rows and rerun. A common fix is to point `universal.negative` to your in-folder unstained control file.

#### Path B: Full Production Workflow (with optional GUI adjustment)

1. Build reference matrix from SCCs.
2. If needed, open GUI and adjust matrix values.
3. Save adjusted matrix CSV from GUI.
4. Load adjusted matrix, then unmix experimental samples.
5. Generate sample QC report.

**After GUI adjustment, continue with:**

```r
# Read adjusted matrix CSV exported from GUI
M_adj_df <- read.csv("reference_matrix_adjusted.csv", stringsAsFactors = FALSE, check.names = FALSE)
M_adj <- as.matrix(M_adj_df[, -1, drop = FALSE])
rownames(M_adj) <- M_adj_df[[1]]

# Unmix experimental samples with adjusted matrix
unmixed <- unmix_samples(
  sample_dir = "samples",
  M = M_adj,
  method = "WLS",
  cytometer = "Aurora",
  output_dir = "samples_unmixed"
)

# QC report
generate_sample_qc(
  unmixed_list = unmixed,
  M = M_adj,
  report_file = "Experimental_Sample_Audit.pdf"
)
```

### Step 1: Build Reference Matrix

Extract spectral signatures from single-color controls:

```r
library(spectrQC)

# Load control file (optional but recommended)
control_df <- read.csv("fcs_control_file.csv", stringsAsFactors = FALSE, check.names = FALSE)

# Build reference matrix from SCC files
M <- build_reference_matrix(
  input_folder = "scc",
  output_folder = "gating_plots",
  control_df = control_df,
  default_sample_type = "beads",
  cytometer = "Aurora"
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
  cytometer = "Aurora",
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

For manual fine-tuning, use the web interface.

### First-time setup (once):

```bash
# From the spectrQC package directory
cd gui
npm install
```

### Launch the GUI:

```r
launch_gui(
  matrix_dir = getwd(),
  samples_dir = "samples"
)
```

This starts both the backend API and frontend automatically, opening http://localhost:5174 in your browser.

### What To Do After GUI

1. Save your adjusted matrix as a CSV (for example `reference_matrix_adjusted.csv`).
2. Load that CSV in R and convert to a matrix (`rownames = first column`).
3. Run `unmix_samples(...)` on your experimental samples using the adjusted matrix.
4. Run `generate_sample_qc(...)` for the final audit report.

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
