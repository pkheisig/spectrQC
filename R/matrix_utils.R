.as_reference_matrix <- function(M, arg_name = "M") {
    if (is.null(M)) {
        stop(arg_name, " must not be NULL.")
    }
    
    if (is.matrix(M)) {
        if (!is.numeric(M)) {
            suppressWarnings(storage.mode(M) <- "numeric")
        }
        if (!is.numeric(M)) {
            stop(arg_name, " must be numeric.")
        }
        return(M)
    }
    
    if (inherits(M, "data.frame")) {
        df <- as.data.frame(M, stringsAsFactors = FALSE)
        if (ncol(df) == 0) {
            stop(arg_name, " has no columns.")
        }
        
        lower_names <- tolower(colnames(df))
        marker_idx <- match("marker", lower_names)
        if (is.na(marker_idx)) marker_idx <- match("fluorophore", lower_names)
        if (is.na(marker_idx)) marker_idx <- match("file", lower_names)
        if (is.na(marker_idx) && ncol(df) > 0 && !is.numeric(df[[1]])) {
            marker_idx <- 1L
        }
        
        if (!is.na(marker_idx)) {
            marker_names <- as.character(df[[marker_idx]])
            df <- df[, -marker_idx, drop = FALSE]
        } else {
            marker_names <- rownames(df)
            if (is.null(marker_names) || !any(nzchar(marker_names))) {
                marker_names <- as.character(seq_len(nrow(df)))
            }
        }
        
        if (ncol(df) == 0) {
            stop(arg_name, " has no detector columns after removing marker labels.")
        }
        
        numeric_cols <- lapply(df, function(x) suppressWarnings(as.numeric(x)))
        bad_cols <- vapply(seq_along(numeric_cols), function(i) {
            original <- df[[i]]
            converted <- numeric_cols[[i]]
            any(!is.na(original)) && all(is.na(converted))
        }, logical(1))
        
        if (any(bad_cols)) {
            stop(
                arg_name,
                " contains non-numeric detector columns: ",
                paste(names(df)[bad_cols], collapse = ", ")
            )
        }
        
        mat <- as.matrix(as.data.frame(numeric_cols, stringsAsFactors = FALSE))
        rownames(mat) <- marker_names
        colnames(mat) <- colnames(df)
        return(mat)
    }
    
    stop(arg_name, " must be a numeric matrix or a data.frame with detector columns.")
}

.is_passthrough_parameter <- function(param_names) {
    if (length(param_names) == 0) {
        return(logical(0))
    }

    param_names <- trimws(as.character(param_names))
    is_scatter <- grepl("^(FSC|SSC)", param_names, ignore.case = TRUE)
    is_time <- grepl("^TIME($|[^A-Z0-9])", param_names, ignore.case = TRUE)

    is_scatter | is_time
}

.get_passthrough_parameter_names <- function(param_names, detector_names = character()) {
    param_names <- as.character(param_names)
    detector_names <- as.character(detector_names)

    keep <- param_names[.is_passthrough_parameter(param_names)]
    keep[!(keep %in% detector_names)]
}

.append_passthrough_parameters <- function(out, full_data, detector_names = character()) {
    keep <- .get_passthrough_parameter_names(colnames(full_data), detector_names = detector_names)
    if (length(keep) == 0) {
        return(out)
    }

    out[keep] <- as.data.frame(full_data[, keep, drop = FALSE], check.names = FALSE)
    out
}

.get_result_metadata_columns <- function(col_names) {
    unique(c(
        "File",
        .get_passthrough_parameter_names(col_names)
    ))
}

.get_primary_scatter_channels <- function(col_names) {
    pick_primary <- function(prefix) {
        exact <- col_names[toupper(col_names) == paste0(prefix, "-A")]
        if (length(exact) > 0) {
            return(exact[[1]])
        }

        area <- col_names[grepl(paste0("^", prefix, ".*-A$"), col_names, ignore.case = TRUE)]
        if (length(area) > 0) {
            return(area[[1]])
        }

        any_match <- col_names[grepl(paste0("^", prefix), col_names, ignore.case = TRUE)]
        if (length(any_match) > 0) {
            return(any_match[[1]])
        }

        NA_character_
    }

    list(
        fsc = pick_primary("FSC"),
        ssc = pick_primary("SSC")
    )
}
