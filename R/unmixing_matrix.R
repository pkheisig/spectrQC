#' Derive Unmixing Matrix from Reference Matrix
#' 
#' Calculates the analytical unmixing matrix from a reference spectral matrix (M).
#' The output matrix (W) has the same dimensions as M (Markers x Detectors).
#' To unmix data manually: Unmixed_Data = Raw_Data %*% t(W)
#' 
#' @param M Reference matrix (Markers x Detectors)
#' @param method Unmixing method ("OLS" or "WLS")
#' @param global_weights Optional vector of weights for WLS (one per detector)
#' @return A matrix of unmixing coefficients (Markers x Detectors)
#' @export
derive_unmixing_matrix <- function(M, method = "OLS", global_weights = NULL) {
    # M is Markers (m) x Detectors (d)
    Mt <- t(M) # Detectors x Markers
    
    if (toupper(method) == "OLS") {
        # Standard analytical solution: W = (M %*% M^T)^-1 %*% M
        # This W is m x d.
        # Check for singularity
        MMt <- M %*% Mt
        if (rcond(MMt) < 1e-10) stop("Reference Matrix is singular (collinear spectra).")
        
        W <- solve(MMt) %*% M
        
    } else if (toupper(method) == "WLS") {
        # Weighted solution: W = (M %*% V^-1 %*% M^T)^-1 %*% M %*% V^-1
        # If no weights provided, we assume Identity (which is OLS)
        if (is.null(global_weights)) {
            warning("WLS requested but no weights provided. Using OLS.")
            return(derive_unmixing_matrix(M, method = "OLS"))
        }
        
        V_inv <- diag(global_weights)
        
        # W = (M V_inv Mt)^-1 M V_inv
        MVMt <- M %*% V_inv %*% Mt
        if (rcond(MVMt) < 1e-10) stop("Weighted Matrix is singular.")
        
        W <- solve(MVMt) %*% M %*% V_inv
    } else {
        stop("Only 'OLS' and 'WLS' (with global weights) can generate a single static matrix.")
    }
    
    rownames(W) <- rownames(M)
    colnames(W) <- colnames(M)
    
    return(W)
}

#' Save Unmixing Matrix to CSV
#' @param W Unmixing matrix
#' @param file Path to save
#' @export
save_unmixing_matrix <- function(W, file = "unmixing_matrix.csv") {
    W_df <- as.data.frame(W)
    W_df$Marker <- rownames(W)
    W_df <- W_df[, c("Marker", setdiff(colnames(W_df), "Marker"))]
    data.table::fwrite(W_df, file)
    message("Unmixing matrix (", nrow(W), "x", ncol(W), ") saved to: ", file)
}
