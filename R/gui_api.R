# gui_api_adjust.R
# API for the spectrQC Interactive Tuner (Adjustment/Crosstalk Correction)

library(plumber)
library(data.table)
library(flowCore)
# spectrQC must be loaded via devtools::load_all() before running this API

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
    matrix_dir <- Sys.getenv("SPECTRQC_MATRIX_DIR", getwd())
    files <- list.files(matrix_dir, pattern = "\\.csv$", full.names = FALSE)
    return(as.character(files))
}

#* Load a specific matrix
#* @get /load_matrix
#* @param filename
function(filename) {
    matrix_dir <- Sys.getenv("SPECTRQC_MATRIX_DIR", getwd())
    path <- file.path(matrix_dir, filename)
    if (!file.exists(path)) {
        return(list(error = paste("File not found:", path)))
    }

    df <- fread(path)

    # Check if it's a spillover/reference matrix (M) or unmixing matrix (W)
    # Reference/Spillover usually has Markers as rows, Detectors as cols
    # Unmixing usually has Detectors as rows, Markers as cols (or vice versa depending on convention)
    # spectrQC convention:
    # M (Reference): Rows=Markers, Cols=Detectors
    # W (Unmixing): Rows=Markers, Cols=Detectors (so Unmixed = Raw %*% t(W))

    # Ensure the first column is named 'Marker' for the frontend
    # If the CSV has a header but first col is unnamed or "V1", fix it
    if (colnames(df)[1] %in% c("V1", "")) {
        setnames(df, colnames(df)[1], "Marker")
    } else if (colnames(df)[1] != "Marker") {
        # Assume first column is Marker if it contains strings
        if (is.character(df[[1]])) {
            setnames(df, colnames(df)[1], "Marker")
        }
    }
    return(df)
}

#* Save the adjusted matrix
#* @post /save_matrix
#* @param filename The filename to save as
#* @param matrix_json The matrix data as JSON
function(filename, matrix_json) {
    dt <- rbindlist(matrix_json, fill = TRUE)
    matrix_dir <- Sys.getenv("SPECTRQC_MATRIX_DIR", getwd())
    path <- file.path(matrix_dir, filename)
    fwrite(dt, path)
    return(list(success = TRUE, path = path))
}

#* Get unmixed data for a subset of cells from a sample
#* @get /data
#* @param sample_name The name of the sample to load
function(sample_name = NULL) {
    samples_dir <- Sys.getenv("SPECTRQC_SAMPLES_DIR", file.path(getwd(), "samples"))
    files <- list.files(samples_dir, pattern = "fcs", full.names = TRUE)
    if (length(files) == 0) {
        return(list(error = paste("No FCS files found in", samples_dir)))
    }

    if (is.null(sample_name)) {
        sample_path <- files[1]
    } else {
        sample_path <- file.path(samples_dir, paste0(sample_name, ".fcs"))
    }

    if (!file.exists(sample_path)) {
        return(list(error = paste("Sample not found:", sample_path)))
    }

    ff <- read.FCS(sample_path, transformation = FALSE, truncate_max_range = FALSE)
    raw_data <- exprs(ff)

    # Subsample for speed - smaller for fast interactive updates
    n_sub <- 2000
    if (nrow(raw_data) > n_sub) {
        set.seed(123)
        raw_data <- raw_data[sample(nrow(raw_data), n_sub), ]
    }

    pd <- pData(parameters(ff))
    # Helper to get sorted detectors (copying logic from spectrQC if not exported)
    # Assuming spectrQC is loaded or we implement basic logic
    det_info <- tryCatch(
        {
            spectrQC:::get_sorted_detectors(pd)
        },
        error = function(e) {
            # Fallback if function not accessible
            fl_cols <- grep("FL", pd$name, value = TRUE)
            list(names = fl_cols, labels = pd$desc[match(fl_cols, pd$name)])
        }
    )

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

#* Run unmixing (On-demand unmixing endpoint)
#* @post /unmix
#* @param matrix_json The matrix (M or W)
#* @param raw_data_json The raw data
#* @param type "reference" (M) or "unmixing" (W)
function(matrix_json, raw_data_json, type = "reference") {
    # matrix_json format: {MarkerName: {det1: val, det2: val, ...}, ...}
    markers <- names(matrix_json)
    detectors <- names(matrix_json[[1]])

    mat <- matrix(0, nrow = length(markers), ncol = length(detectors))
    for (i in seq_along(markers)) {
        mat[i, ] <- as.numeric(unlist(matrix_json[[markers[i]]][detectors]))
    }
    rownames(mat) <- markers
    colnames(mat) <- detectors

    Y <- as.matrix(as.data.frame(raw_data_json))

    # Matching columns
    common_dets <- intersect(colnames(Y), colnames(mat))
    if (length(common_dets) == 0) {
        return(list(error = "No matching detectors found between data and matrix"))
    }

    Y_sub <- Y[, common_dets, drop = FALSE]
    mat_sub <- mat[, common_dets, drop = FALSE]

    # Perform Unmixing
    if (grepl("unmixing", tolower(type))) {
        # Provided matrix IS the unmixing matrix (W)
        # Unmixed = Raw * t(W)
        unmixed <- Y_sub %*% t(mat_sub)
    } else {
        # Provided matrix IS the reference matrix (M), derive W via OLS
        # W = (M M^t)^-1 M  ? No.
        # Simple OLS: Y ~ U * M  => U = Y * M^T * (M * M^T)^-1 ??
        # Standard OLS unmixing: U = Y * pseudoinverse(M)
        # U = Y * pinv(M) = Y * t(M) * (M * t(M))^-1  (if M is tall? No M is Markers x Detectors, usually Fat)
        # Actually in spectral: Y (Cells x Dets) = U (Cells x Markers) * M (Markers x Dets)
        # U = Y * M_pinv
        # M_pinv = M^T (M M^T)^-1
        # So W^T = M^T (M M^T)^-1
        # W = (M M^T)^-1 M

        W <- tryCatch(
            {
                spectrQC::derive_unmixing_matrix(mat_sub, method = "OLS")
            },
            error = function(e) {
                # Fallback OLS
                t(MASS::ginv(t(mat_sub)))
            }
        )

        unmixed <- Y_sub %*% t(W)
    }

    return(as.data.frame(unmixed))
}
