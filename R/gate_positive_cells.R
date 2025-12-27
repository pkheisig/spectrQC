gate_positive_cells <- function(mat,
                                histogram_pct = 0.98,
                                histogram_direction = "both",
                                histogram_min_x_log = 2) {
    # Identify peak channel by variance
    peak_channel <- which.max(apply(mat, 2, var))
    peak_vals <- mat[, peak_channel]
    vals_log <- log10(pmax(peak_vals, 1))

    # Calculate gate thresholds based on direction
    if (histogram_direction == "right") {
        lower_q <- 0.5
        upper_q <- 0.5 + histogram_pct
    } else if (histogram_direction == "left") {
        lower_q <- 0.5 - histogram_pct
        upper_q <- 0.5
    } else { # "both"
        lower_q <- 0.5 - histogram_pct / 2
        upper_q <- 0.5 + histogram_pct / 2
    }
    lower_q <- max(0, lower_q)
    upper_q <- min(1, upper_q)

    # Apply thresholds on log scale, return linear
    gate_min <- 10^quantile(vals_log, lower_q)
    gate_max <- 10^quantile(vals_log, upper_q)

    peak_vals >= gate_min & peak_vals <= gate_max
}
