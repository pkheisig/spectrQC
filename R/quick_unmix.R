#' Quick Unmix (Compatibility Wrapper)
#'
#' Backward-compatible wrapper around [autounmix_controls()].
#' For new workflows, call [autounmix_controls()] directly.
#'
#' @inheritParams autounmix_controls
#' @return Same as [autounmix_controls()].
#' @export
#' @examples
#' \dontrun{
#' out <- quick_unmix(
#'   scc_dir = "scc",
#'   control_file = "fcs_control_file.csv",
#'   auto_create_control = TRUE
#' )
#' names(out)
#' }
quick_unmix <- function(
    scc_dir = "scc",
    control_df = NULL,
    control_file = "fcs_control_file.csv",
    auto_create_control = TRUE,
    cytometer = "Aurora",
    auto_default_control_type = "beads",
    auto_unknown_fluor_policy = c("by_channel", "empty", "filename"),
    output_dir = "spectrQC_outputs/quick_unmix",
    unmix_method = "WLS",
    build_qc_plots = FALSE,
    unmix_scatter_panel_size_mm = 30,
    ...
) {
    .Deprecated("autounmix_controls")
    autounmix_controls(
        scc_dir = scc_dir,
        control_df = control_df,
        control_file = control_file,
        auto_create_control = auto_create_control,
        cytometer = cytometer,
        auto_default_control_type = auto_default_control_type,
        auto_unknown_fluor_policy = auto_unknown_fluor_policy,
        output_dir = output_dir,
        unmix_method = unmix_method,
        build_qc_plots = build_qc_plots,
        unmix_scatter_panel_size_mm = unmix_scatter_panel_size_mm,
        ...
    )
}
