devtools::load_all(".")

quick_unmix(
    scc_dir = "scc",
    control_file = "fcs_control_file.csv",
    auto_create_control = TRUE,
    cytometer = "Aurora",
    auto_default_control_type = "beads",
    auto_unknown_fluor_policy = "by_channel",
    output_dir = "spectrQC_outputs/quick_unmix",
    unmix_method = "WLS",
    build_qc_plots = TRUE
)
