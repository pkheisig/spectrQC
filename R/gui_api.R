# gui_api.R
# API for the spectrQC Interactive Tuner

library(plumber)
library(data.table)
library(flowCore)

#* @filter logger
function(req) {
    cat(as.character(Sys.time()), "-", req$REQUEST_METHOD, req$PATH_INFO, "\n")
    forward()
}

#* @filter cors
function(res) {
    res$setHeader("Access-Control-Allow-Origin", "*")
    res$setHeader("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
    res$setHeader("Access-Control-Allow-Headers", "Content-Type")
    plumber::forward()
}

#* Health Check
#* @get /status
function() {
    return(list(status = "ok", time = Sys.time(), wd = getwd()))
}

#* List available matrices
#* @get /matrices
function() {
    # Look in the parent directory (project root)
    files <- list.files("..", pattern = ".*\\.csv$")
    # Exclude control file
    files <- files[!grepl("control_file", files)]
    return(as.character(files))
}

#* Load a specific matrix
#* @get /load_matrix
#* @param filename
function(filename) {
    path <- file.path("..", filename)
    if (!file.exists(path)) {
        return(list(error = paste("File not found:", path)))
    }

    df <- fread(path)
    # Ensure the first column is named 'Marker' for the frontend
    if (colnames(df)[1] != "Marker") {
        setnames(df, colnames(df)[1], "Marker")
    }
    return(df)
}

#* Get unmixed data for a subset of cells from a sample
#* @get /data
#* @param sample_name The name of the sample to load
function(sample_name = NULL) {
    # If no sample provided, take the first one in samples/
    files <- list.files("../samples", pattern = "fcs", full.names = TRUE)
    if (length(files) == 0) {
        return(list(error = "No FCS files found in samples/"))
    }

    if (is.null(sample_name)) {
        sample_path <- files[1]
    } else {
        sample_path <- file.path("../samples", paste0(sample_name, ".fcs"))
    }

    if (!file.exists(sample_path)) {
        return(list(error = paste("Sample not found:", sample_path)))
    }

    ff <- read.FCS(sample_path, transformation = FALSE, truncate_max_range = FALSE)
    raw_data <- exprs(ff)
    if (nrow(raw_data) > 5000) {
        raw_data <- raw_data[sample(nrow(raw_data), 5000), ]
    }

    pd <- pData(parameters(ff))
    det_info <- get_sorted_detectors(pd)

    return(list(
        raw_data = as.data.frame(raw_data),
        detector_names = det_info$names,
        detector_labels = det_info$labels
    ))
}

#* CORS preflight for unmix
#* @options /unmix
function(res) {
    return("")
}

#* Run unmixing
#* @post /unmix
#* @param matrix_json The matrix (M or W)
#* @param raw_data_json The raw data
#* @param type "reference" or "unmixing"
function(matrix_json, raw_data_json, type = "reference") {
    markers <- names(matrix_json)
    detectors <- names(matrix_json[[1]])

    mat <- matrix(0, nrow = length(markers), ncol = length(detectors))
    for (i in seq_along(markers)) {
        mat[i, ] <- as.numeric(unlist(matrix_json[[markers[i]]][detectors]))
    }
    rownames(mat) <- markers
    colnames(mat) <- detectors

    Y <- as.matrix(as.data.frame(raw_data_json))

    # Matching
    common_dets <- intersect(colnames(Y), colnames(mat))
    Y_sub <- Y[, common_dets, drop = FALSE]
    mat_sub <- mat[, common_dets, drop = FALSE]

    if (grepl("unmixing", tolower(type))) {
        # Provided matrix IS the unmixing matrix (W)
        unmixed <- Y_sub %*% t(mat_sub)
    } else {
        # Provided matrix IS the reference matrix (M), derive W
        W <- spectrQC::derive_unmixing_matrix(mat_sub, method = "OLS")
        unmixed <- Y_sub %*% t(W)
    }

    return(as.data.frame(unmixed))
}
