library(flowCore)
library(data.table)
library(devtools)
load_all()

# 1. Load Reference Matrix
M <- as.matrix(fread("reference_matrix.csv", select = 2:52))
rownames(M) <- fread("reference_matrix.csv")$file

# 2. Process a sample with detailed residuals
files_raw <- list.files("raw", pattern = ".*.fcs$", full.names = TRUE)[1:2]
res_list <- lapply(files_raw, function(f) {
    ff <- read.FCS(f)
    if (nrow(ff) > 5000) ff <- ff[sample(nrow(ff), 5000), ]
    calc_residuals(ff, M, file_name = basename(f), method = "WLS", return_residuals = TRUE)
})

# Combine data for NPS
all_data <- rbindlist(lapply(res_list, `[[`, "data"))

# 3. Test NPS
message("Calculating NPS...")
nps_scores <- calculate_nps(all_data)
plot_nps(nps_scores, "qc_nps.png")

# 4. Test Detector Residuals (on the first file)
message("Plotting detector residuals...")
plot_detector_residuals(res_list[[1]], top_n = 100, output_file = "qc_detector_res.png")

print("New QC tests complete. Files: qc_nps.png, qc_detector_res.png")
