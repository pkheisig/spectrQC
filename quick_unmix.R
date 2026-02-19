devtools::load_all(".")

autounmix_controls(
    scc_dir = "scc",
    control_file = "fcs_control_file.csv",
    auto_create_control = TRUE,
    cytometer = "Aurora",
    auto_unknown_fluor_policy = "by_channel",
    output_dir = "spectrQC_outputs/autounmix_controls",
    unmix_method = "WLS",
    build_qc_plots = TRUE
)
