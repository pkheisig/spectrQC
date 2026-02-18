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
    utils::write.csv(W_df, file, row.names = FALSE, quote = TRUE)
    message("Unmixing matrix (", nrow(W), "x", ncol(W), ") saved to: ", file)
}

#' Plot Unmixing Matrix
#' @param W Unmixing matrix
#' @param pd Optional pData for descriptive labels
#' @return ggplot object
#' @export
plot_unmixing_matrix <- function(W, pd = NULL) {
    long <- as.data.frame(W)
    long$Marker <- rownames(W)
    long <- tidyr::pivot_longer(long, cols = -Marker, names_to = "Detector", values_to = "Coefficient")
    
    # Sort and label detectors
    det_names <- colnames(W)
    if (!is.null(pd)) {
        det_info <- get_sorted_detectors(pd)
        common <- intersect(det_info$names, det_names)
        levels_sorted <- common
        labels_sorted <- det_info$labels[match(common, det_info$names)]
    } else {
        nums <- as.numeric(gsub("[^0-9]", "", det_names))
        levels_sorted <- det_names[order(nums)]
        labels_sorted <- levels_sorted
    }
    
    long$Detector <- factor(long$Detector, levels = levels_sorted, labels = labels_sorted)
    
    # Cap values for color scale to avoid outliers dominating
    cap_val <- quantile(abs(long$Coefficient), 0.95, na.rm=TRUE)
    if (is.na(cap_val) || cap_val == 0) cap_val <- 1
    long$Color_Val <- pmin(pmax(long$Coefficient, -cap_val), cap_val)
    
    p <- ggplot2::ggplot(long, ggplot2::aes(Detector, Marker, fill = Color_Val)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", name = "Coefficient") +
        # Text color based on value
        ggplot2::geom_text(ggplot2::aes(label = round(Coefficient, 3), 
                                       color = abs(Color_Val) > cap_val/2), 
                           size = 1.5, show.legend = FALSE, angle = 90) +
        ggplot2::scale_color_manual(values = c("TRUE" = "white", "FALSE" = "black")) +
        ggplot2::labs(title = "Unmixing Matrix Coefficients",
                      subtitle = "Positive coefficients (red) indicate detector contribution to marker recovery.") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 5))
    
    return(p)
}
