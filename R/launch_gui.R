#' Launch SpectrQC Interactive Adjustment GUI
#'
#' Starts the backend Plumber API for the interactive matrix adjustment interface.
#' The frontend must be started separately by running 'npm run dev' in the gui folder.
#'
#' @param matrix_dir Directory containing matrix CSV files (default: current working directory)
#' @param samples_dir Directory containing FCS sample files (default: "samples" subfolder of matrix_dir)
#' @param port API port (default: 8000)
#' @param open_browser Logical. Open browser automatically? (default: TRUE)
#' @return Invisibly returns NULL. This function blocks while the API is running.
#' @export
#' @examples
#' \dontrun{
#' # Start the API with default directories
#' launch_gui()
#'
#' # Start with custom matrix directory
#' launch_gui(matrix_dir = "/path/to/my/matrices")
#'
#' # Then in a separate terminal, run:
#' # cd /path/to/spectrQC/gui && npm run dev
#' # Open localhost:5174 in your browser
#' }
launch_gui <- function(matrix_dir = getwd(), samples_dir = NULL, port = 8000, open_browser = TRUE) {
    api_path <- system.file("api/gui_api.R", package = "spectrQC")
    gui_path <- system.file("gui", package = "spectrQC")

    if (api_path == "") {
        api_path <- file.path(getwd(), "inst", "api", "gui_api.R")
    }
    if (gui_path == "" || !dir.exists(gui_path)) {
        gui_path <- file.path(getwd(), "inst", "gui")
    }
    if (!dir.exists(gui_path)) {
        gui_path <- file.path(getwd(), "gui")
    }

    if (!file.exists(api_path)) stop("Could not find gui_api.R")
    if (!dir.exists(gui_path)) stop("Could not find gui folder")

    node_modules <- file.path(gui_path, "node_modules")
    if (!dir.exists(node_modules)) {
        stop("node_modules not found. Run this once in terminal:\n  cd ", gui_path, " && npm install")
    }

    matrix_dir <- normalizePath(matrix_dir, mustWork = TRUE)
    if (is.null(samples_dir)) samples_dir <- file.path(matrix_dir, "samples")
    samples_dir <- normalizePath(samples_dir, mustWork = FALSE)

    options(
        spectrqc.matrix_dir = matrix_dir,
        spectrqc.samples_dir = samples_dir
    )

    message("Starting frontend (npm run dev)...")
    old_wd <- getwd()
    setwd(gui_path)
    system2("npm", args = c("run", "dev"), wait = FALSE, stdout = FALSE, stderr = FALSE)
    setwd(old_wd)

    Sys.sleep(2)

    message("Starting spectrQC API on port ", port)
    message("Matrix directory: ", matrix_dir)
    message("Samples directory: ", samples_dir)
    message("Frontend: localhost:5174")

    if (open_browser) utils::browseURL("http://localhost:5174")

    pr <- plumber::plumb(api_path)
    pr$run(port = port, host = "0.0.0.0")

    invisible(NULL)
}
