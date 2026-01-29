# Benchmark script for WLS optimization
# Run this to measure performance improvement

# Note: This script assumes spectrQC is installed or loaded with C++ extensions compiled.

if (!requireNamespace("flowCore", quietly = TRUE)) {
  stop("flowCore needed for benchmark")
}

# 1. Create Mock Data
set.seed(42)
n_cells <- 5000 # Adjusted for quick check
n_fluor <- 16
n_detectors <- 64

# Random reference matrix M (fluorophores x detectors)
# Ensure full rank
M <- matrix(rnorm(n_fluor * n_detectors, mean = 100, sd = 20), nrow = n_fluor, ncol = n_detectors)
M <- abs(M)
rownames(M) <- paste0("Fluor", 1:n_fluor)
colnames(M) <- paste0("Det", 1:n_detectors)

# Random data Y (cells x detectors)
Y <- matrix(rnorm(n_cells * n_detectors, mean = 50, sd = 10), nrow = n_cells, ncol = n_detectors)
Y <- abs(Y)
colnames(Y) <- colnames(M)

# Create flowFrame
ff <- flowCore::flowFrame(exprs = Y)

# 2. Benchmark WLS
cat("Running WLS benchmark with", n_cells, "cells, ", n_fluor, " fluorophores, ", n_detectors, " detectors...\n")

# To benchmark properly, one would need to compare against the old R version.
# Since we replaced it, this just measures current performance.

start_time <- Sys.time()
# Assuming spectrQC is loaded
tryCatch({
    res <- spectrQC::calc_residuals(ff, M, method = "WLS", background_noise = 25)
    end_time <- Sys.time()

    duration <- as.numeric(end_time - start_time, units = "secs")
    cat("WLS Duration:", duration, "seconds\n")
    cat("Average time per cell:", (duration / n_cells) * 1000, "ms\n")
}, error = function(e) {
    cat("Error running benchmark (likely package not installed/compiled): ", e$message, "\n")
})
