library(devtools)
load_all()
library(data.table)
library(flowCore)

# --- spectrQC: New Modular Workflow ---

# Setup folders
scc_dir <- "scc"       # User's SCC folder
sample_dir <- "samples" # User's experimental folder
control_file <- "fcs_control_file.csv" # User's manual mapping

# --- Step 1: Check SCCs ---
message("Step 1: Inspecting SCC Spectra...")
M_initial <- inspect_scc_spectra(scc_dir = scc_dir, 
                                output_dir = "gating_and_spectrum_plots", 
                                control_file = control_file)

message("  - Generating SCC initial report...")
generate_scc_report(M_initial, scc_dir = scc_dir, output_file = "SCC_Initial_Audit.pdf")


# --- Step 2: Refine SCCs (Cleaning) ---
message("\nStep 2: Refining SCC Matrix (Removing outliers)...")
# User picks 5% cutoff based on Step 1 report
refined <- refine_scc_matrix(M_initial, 
                            scc_dir = scc_dir, 
                            rrmse_threshold = 0.05, 
                            output_dir = "scc_unmixed",
                            report_file = "SCC_Refinement_Report.pdf")

M_final <- refined$M
W_final <- refined$W


# --- Step 3: Unmix Samples ---
message("\nStep 3: Unmixing Experimental Samples...")
unmixed_results <- unmix_samples(sample_dir = sample_dir, 
                               M = M_final, 
                               method = "WLS", 
                               output_dir = "samples_unmixed")


# --- Step 4: QC Samples ---
message("\nStep 4: Generating Experimental Sample QC...")
# Get metadata for labels
ff_qc <- read.FCS(list.files(sample_dir, pattern="\\.fcs$", full.names=TRUE)[1])
pd_qc <- pData(parameters(ff_qc))

generate_sample_qc(unmixed_results, M = M_final, 
                  output_dir = "samples_unmixed_qc", 
                  report_file = "Experimental_Sample_Audit.pdf",
                  pd = pd_qc)

message("\nWorkflow Complete.")
message("Reports generated: SCC_Initial_Audit.pdf, SCC_Refinement_Report.pdf, Experimental_Sample_Audit.pdf")
