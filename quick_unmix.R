devtools::load_all(".")

control_df <- read.csv("fcs_control_file.csv", stringsAsFactors = FALSE, check.names = FALSE)

quick_unmix(
    scc_dir = "scc",
    control_df = control_df,
    output_dir = "spectrQC_outputs/quick_unmix",
    unmix_method = "WLS",
    build_qc_plots = TRUE
)
