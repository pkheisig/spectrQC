# spectrQC: Full Spectrum Flow Cytometry Quality Control

`spectrQC` is an R package for validating spectral unmixing accuracy. It provides a complete pipeline to extract fluorophore signatures from single-color controls, refine them using quality metrics, and unmix experimental samples.

## Key Features

- **Automated Gating**: Isolate positive populations from beads or cells using Gaussian Mixture Models
- **Background Subtraction**: Automatically subtract internal negative populations to isolate pure fluorophore signatures
- **Signature Refinement**: Remove debris and outliers using RRMSE-based filtering
- **Per-Cell WLS Unmixing**: High-accuracy unmixing using photon-counting variance weighting
- **SCC Diagnostics & Visualization**: Spectra and SCC unmixing scatter outputs for control-stage QC
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

<img width="1102" height="201" alt="Screenshot 2026-02-20 at 11 00 23" src="https://github.com/user-attachments/assets/a8c5e253-e1c8-4592-8ea2-e0f404932a54" />


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

### Path A: Recommended Quick Workflow

#### Step 1: Run controls with `autounmix_controls()`:

```r
library(spectrQC)

ctrl <- autounmix_controls(
  scc_dir = "scc",
  control_file = "fcs_control_file.csv",
  auto_create_control = TRUE,
  cytometer = "Aurora",
  auto_unknown_fluor_policy = "by_channel",
  output_dir = "spectrQC_outputs/autounmix_controls",
  unmix_method = "WLS",
  build_qc_plots = TRUE,
  unmix_scatter_panel_size_mm = 30
)
```

If `fcs_control_file.csv` is missing and `auto_create_control = TRUE`, `autounmix_controls()` auto-generates a control file (filename, marker, fluorophore, and detected peak channel), then asks for confirmation before continuing.
`cytometer` is used for channel-aware fluorophore inference via the AutoSpectral fluorophore database.

For newly auto-created files:
- `control.type` is set to `cells` only for AF rows; all non-AF rows are left empty on purpose
- `universal.negative` is left empty by default for all rows
- `autounmix_controls()` pauses and asks for `y/n` confirmation so you can review/edit the file first

`autounmix_controls()` writes:
- `scc_reference_matrix.csv`
- `scc_spectra.png` (reference spectra overlay)
- `scc_unmixing_matrix.png` and `scc_unmixing_matrix.csv`
- `scc_unmixing_scatter_matrix.png` (lower-triangle scatter matrix, one single-stain file per row, with x=0/y=0 guides)

Set `unmix_scatter_panel_size_mm` higher (for example `40`) if you want larger per-panel scatter plots.

`autounmix_controls()` also runs a strict preflight check before processing:
- every SCC file must be mapped in `fcs_control_file.csv`
- non-AF rows must define a valid `channel`
- if `universal.negative` is present, values for active SCC rows must be empty/keyword or reference a file present in your selected SCC/AF directories

---

#### Step 2: Unmix samples using the unmixing matrix generated in the autounmix_controls step.

```r
# Uses saved unmixing matrix by filepath (default points to autounmix_controls output)
unmixed <- unmix_samples(
  sample_dir = "samples",
  unmixing_matrix_file = "spectrQC_outputs/autounmix_controls/scc_unmixing_matrix.csv",
  output_dir = "spectrQC_outputs/unmix_samples"
)
```

### Path B: Step-wise manual workflow

#### Step 1: Build Reference Matrix

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

This saves gating/spectrum plots to `gating_plots/` and exports the matrix to `spectrQC_outputs/reference_matrix.csv`.

---

Play around with gating parameters if auto-gating fails:

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

#### Step 2: Refine the Matrix

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

#### Step 3: Unmix Experimental Samples

Apply the refined matrix to your samples:

```r
# Option 1: dynamic unmixing directly from reference matrix (M)
unmixed <- unmix_samples(
  sample_dir = "samples",
  M = M_final,
  method = "WLS",                     # "OLS", "WLS", or "NNLS"
  cytometer = "Aurora",
  output_dir = "spectrQC_outputs/unmix_samples"
)

# Option 2: static unmixing from saved unmixing matrix (W)
unmixed_w <- unmix_samples(
  sample_dir = "samples",
  unmixing_matrix_file = "spectrQC_outputs/autounmix_controls/scc_unmixing_matrix.csv",
  output_dir = "spectrQC_outputs/unmix_samples_w"
)
```

**Methods:**
- **OLS**: Ordinary least squares — fast, suitable for most panels
- **WLS**: Weighted least squares — accounts for photon-counting noise, best accuracy
- **NNLS**: Non-negative least squares — forces positive abundances

---

### Optional: Interactive Matrix Adjustment (before Step 3)

For manual fine-tuning, use the web interface.

No terminal setup is required for end users. The bundled GUI is served directly by `launch_gui()`.

#### Launch the GUI:

```r
launch_gui(
  matrix_dir = getwd(),
  samples_dir = "samples"
)
```

This starts both the backend API and bundled frontend on one port (default `http://localhost:8000`) and opens it in your browser.

Developer mode (optional, for GUI hacking):

```r
launch_gui(
  matrix_dir = getwd(),
  samples_dir = "samples",
  dev_mode = TRUE
)
```

#### What To Do After GUI

1. Save your adjusted matrix CSV (typically `scc_unmixing_matrix.csv` or `scc_reference_matrix.csv`).
2. Run `unmix_samples(...)` on your experimental samples:
   - if you edited `scc_unmixing_matrix.csv`, pass `unmixing_matrix_file = "..."`
   - if you edited a reference matrix, load it as `M` and pass `M = ...`

### Output Directories

- `spectrQC_outputs/reference_matrix.csv`: Reference matrix written by `build_reference_matrix(...)`
- `spectrQC_outputs/autounmix_controls/`: SCC control-stage outputs (`scc_reference_matrix.csv`, `scc_unmixing_matrix.csv/.png`, `scc_spectra.png`, `scc_unmixing_scatter_matrix.png`)
- `spectrQC_outputs/autounmix_controls/scc_unmixed/`: Unmixed SCC control files (FCS format)
- `spectrQC_outputs/unmix_samples/`: Unmixed experimental data (FCS format)

If you run the manual path with `build_reference_matrix(...)`, `output_folder` (for example `gating_plots/`) is used for build-stage QC plots.

---

## Legacy Report APIs (Optional)

The package still exports report helpers for legacy workflows, but report generation is not part of the recommended quick workflow.
Reports now render directly to PDF only (no intermediate PNG files).

```r
# SCC report (singular function name)
generate_scc_report(
  M = ctrl$M,  # reference matrix from build_reference_matrix() or autounmix_controls()$M
  scc_dir = "scc",
  output_file = file.path("spectrQC_outputs", "SCC_QC_Report.pdf")
)

# Full sample-level report
r$> generate_qc_report(
      results_df = do.call(rbind, lapply(unmixed, `[[`, "data")),
      M = ctrl$M,  # matrix used for unmixing context
      output_file = file.path("spectrQC_outputs", "Sample_QC_Report.pdf")
    )
```

---

**Author**: Paul Heisig  
**Email**: p.k.s.heisig@amsterdamumc.nl
