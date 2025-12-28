library(devtools)
load_all()
library(data.table)
library(flowCore)

# Setup folders
scc_dir <- "scc"       
sample_dir <- "samples" 
af_dir <- "af"
control_file <- "fcs_control_file.csv"

# --- Step 1: Check SCCs + Multi-AF ---
message("Step 1: Inspecting SCC and Multi-AF Spectra...")
M_initial <- inspect_scc_spectra(scc_dir = scc_dir, 
                                output_dir = "gating_and_spectrum_plots", 
                                control_file = control_file,
                                include_multi_af = TRUE,
                                af_dir = af_dir)

message("  - Generating SCC initial report...")
generate_scc_report(M_initial, scc_dir = scc_dir, output_file = "SCC_Initial_Audit.pdf", control_file = control_file)


# --- Step 2: Refine signatures ---
message("\nStep 2: Refining signatures...")
refined <- refine_scc_matrix(M_initial, 
                            scc_dir = scc_dir, 
                            rrmse_threshold = 0.05, 
                            output_dir = "scc_unmixed")

M_final <- refined$M


# --- Step 3: Unmix Samples (AutoSpectral) ---
message("\nStep 3: Unmixing Experimental Samples via AutoSpectral math engine...")
unmixed_results <- unmix_samples(sample_dir = sample_dir, 
                               M = M_final, 
                               method = "AutoSpectral", 
                               output_dir = "samples_unmixed_final")


# --- Step 4: QC Samples ---
message("\nStep 4: Generating Final Sample QC...")
# Get metadata for labels from first sample
ff_qc <- read.FCS(list.files(sample_dir, pattern="fcs$", full.names=TRUE)[1])
pd_qc <- pData(parameters(ff_qc))

generate_sample_qc(unmixed_results, M = M_final, 
                  report_file = "Experimental_Sample_Audit.pdf",
                  pd = pd_qc)

message("\nWorkflow Complete.")
message("Reports updated: SCC_Initial_Audit.pdf, Experimental_Sample_Audit.pdf")
