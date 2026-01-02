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
#' # Open http://localhost:5174 in your browser
#' }
launch_gui <- function(matrix_dir = getwd(), samples_dir = NULL, port = 8000, open_browser = TRUE) {
    api_path <- system.file("api/gui_api.R", package = "spectrQC")

    if (api_path == "") {
        api_path <- file.path(getwd(), "inst", "api", "gui_api.R")
        if (!file.exists(api_path)) {
            api_path <- file.path(getwd(), "R", "gui_api.R")
        }
    }

    if (!file.exists(api_path)) {
        stop("Could not find gui_api.R. Make sure spectrQC is installed or you are in the package directory.")
    }

    matrix_dir <- normalizePath(matrix_dir, mustWork = TRUE)
    if (is.null(samples_dir)) {
        samples_dir <- file.path(matrix_dir, "samples")
    }
    samples_dir <- normalizePath(samples_dir, mustWork = FALSE)

    Sys.setenv(SPECTRQC_MATRIX_DIR = matrix_dir)
    Sys.setenv(SPECTRQC_SAMPLES_DIR = samples_dir)

    pr <- plumber::plumb(api_path)

    url <- paste0("http://localhost:", port)
    message("Starting spectrQC API at ", url)
    message("Matrix directory: ", matrix_dir)
    message("Samples directory: ", samples_dir)
    message("")
    message("IMPORTANT: Start the frontend in a separate terminal:")
    message("  cd gui && npm install && npm run dev")
    message("")
    message("Then open http://localhost:5174 in your browser.")
    message("")

    if (open_browser) {
        Sys.sleep(1)
        utils::browseURL("http://localhost:5174")
    }

    pr$run(port = port, host = "0.0.0.0")

    invisible(NULL)
}
