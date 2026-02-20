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
