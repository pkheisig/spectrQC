#' Calculate unmixing residuals
#'
#' @param flow_frame A flowFrame object with raw fluorescence data
#' @param M Reference matrix (fluorophores x detectors)
#' @param file_name Optional file name to add to output
#' @param method Unmixing method: "OLS" (default), "NNLS", or "WLS" (per-cell weighted)
#' @param background_noise Baseline electronic noise variance (default: 25)
#' @param return_residuals Logical. If TRUE, returns a list containing the unmixed data and the raw residual matrix.
#' @return Data frame with unmixed abundances, RMSE, Relative_RMSE, and scatter parameters. 
#'         If return_residuals=TRUE, returns a list with [[data]] and [[residuals]].
#' @export
calc_residuals <- function(flow_frame, M, file_name = NULL, method = "OLS", 
                          background_noise = 25, return_residuals = FALSE) {
    full_data <- flowCore::exprs(flow_frame)
    detectors <- colnames(M)
    
    # Check for missing detectors
    missing <- setdiff(detectors, colnames(full_data))
    if (length(missing) > 0) {
        stop("Detectors in reference matrix not found in flow_frame: ", paste(missing, collapse = ", "))
    }
    
    Y <- full_data[, detectors, drop = FALSE]

    Mt <- t(M)
    n_cells <- nrow(Y)
    n_fluor <- nrow(M)
    
    method <- toupper(method)

    if (method == "OLS") {
        # Ordinary Least Squares: A = Y * M^T * (M * M^T)^-1
        matrix_to_invert <- M %*% Mt
        if (rcond(matrix_to_invert) < 1e-10) {
            stop("Reference Matrix is singular. You likely have collinear spectra.")
        }
        A <- Y %*% Mt %*% solve(matrix_to_invert)
    } else if (method == "NNLS") {
        # Non-Negative Least Squares (Vectorized Coordinate Descent)
        A <- solve_nnls_batch(M, Y)
    } else if (method == "WLS") {
        # Weighted Least Squares (Per-Cell)
        A <- matrix(0, nrow = n_cells, ncol = n_fluor)
        MMt_inv <- solve(M %*% Mt)
        
        for (i in seq_len(n_cells)) {
            weights_i <- 1 / (pmax(Y[i, ], 0) + background_noise)
            Wi <- diag(weights_i)
            MWMt_i <- M %*% Wi %*% Mt
            if (rcond(MWMt_i) > 1e-10) {
                A[i, ] <- Y[i, , drop = FALSE] %*% Wi %*% Mt %*% solve(MWMt_i)
            } else {
                A[i, ] <- Y[i, , drop = FALSE] %*% Mt %*% MMt_inv
            }
        }
    } else {
        stop("method must be 'OLS', 'NNLS', or 'WLS'")
    }

    Fitted <- A %*% M
    R <- Y - Fitted
    rmse <- sqrt(rowMeans(R^2))
    
    # Relative RMSE: RMSE / Total Intensity across detectors
    total_intensity <- rowSums(Y)
    relative_rmse <- rmse / pmax(total_intensity, 1)

    out <- as.data.frame(A)
    colnames(out) <- rownames(M)
    out$RMSE_Score <- rmse
    out$Relative_RMSE <- relative_rmse

    # Add FSC/SSC if available
    all_cols <- colnames(full_data)
    fsc_col <- grep("^FSC[0-9]*-A$", all_cols, value = TRUE)[1]
    ssc_col <- grep("^SSC[0-9]*-A$", all_cols, value = TRUE)[1]
    if (!is.na(fsc_col)) out[["FSC-A"]] <- full_data[, fsc_col]
    if (!is.na(ssc_col)) out[["SSC-A"]] <- full_data[, ssc_col]

    if (!is.null(file_name)) out$File <- file_name

    if (return_residuals) {
        return(list(data = out, residuals = R))
    } else {
        return(out)
    }
}

#' Batch NNLS Solver using Coordinate Descent
#'
#' @param M Reference matrix (fluorophores x detectors)
#' @param Y Observed data (cells x detectors)
#' @param max_iter Maximum iterations
#' @return Matrix of coefficients (cells x fluorophores)
#' @keywords internal
solve_nnls_batch <- function(M, Y, max_iter = 50) {
    XtX <- M %*% t(M)  # n_fluor x n_fluor
    Xty <- Y %*% t(M)  # n_cells x n_fluor

    n_cells <- nrow(Y)
    n_fluor <- nrow(M)

    # Initialization: OLS
    # Handle singularity
    if (rcond(XtX) < 1e-10) {
        A <- matrix(0, nrow = n_cells, ncol = n_fluor)
    } else {
        A <- Xty %*% solve(XtX)
        A[A < 0] <- 0
    }

    # Coordinate Descent
    for (iter in seq_len(max_iter)) {
        for (j in seq_len(n_fluor)) {
            # Update column j
            current_val <- A[, j]
            numerator <- Xty[, j] - (A %*% XtX[, j]) + current_val * XtX[j, j]
            A[, j] <- pmax(0, numerator / XtX[j, j])
        }
    }
    return(A)
}
