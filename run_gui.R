# run_gui.R
library(plumber)
library(devtools)

message("Loading spectrQC package...")
load_all(".")

# Check for refined matrix
if (!file.exists("refined_reference_matrix.csv")) {
    warning("refined_reference_matrix.csv not found. GUI may be empty until a matrix is loaded.")
}

message("Starting spectrQC API on http://localhost:8000")
pr <- plumb("R/gui_api.R")
pr$run(port = 8000)
