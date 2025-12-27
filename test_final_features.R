library(devtools)
load_all()
library(data.table)
library(flowCore)

# 1. Load Reference Matrix
M <- as.matrix(fread("reference_matrix.csv", select = 2:52))
rownames(M) <- fread("reference_matrix.csv")$file

# 2. Test SSM
message("Testing SSM...")
ssm <- calculate_ssm(M)
plot_ssm(ssm, "qc_ssm.png")

# 3. Test Full Analysis + Report
files_raw <- list.files("raw", pattern = ".*.fcs$", full.names = TRUE)[1:2]
res_list <- lapply(files_raw, function(f) {
    ff <- read.FCS(f)
    if (nrow(ff) > 2000) ff <- ff[sample(nrow(ff), 2000), ] # Small sample for report speed
    calc_residuals(ff, M, file_name = basename(f), method = "WLS", return_residuals = TRUE)
})

all_data <- rbindlist(lapply(res_list, `[[`, "data"))

message("Generating full report...")
generate_qc_report(all_data, M, output_file = "specRQC_Full_Test_Report.pdf", res_list = res_list)

print("Final tests complete. Created: qc_ssm.png, specRQC_Full_Test_Report.pdf")
