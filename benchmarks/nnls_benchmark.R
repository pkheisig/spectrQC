# Benchmark for NNLS Optimization
# To run this: Rscript benchmarks/nnls_benchmark.R

set.seed(123)

# 1. Setup Data
n_cells <- 10000  # Number of cells
n_fluor <- 15     # Number of fluorophores
n_detectors <- 40 # Number of detectors

message("Generating data for ", n_cells, " cells, ", n_fluor, " fluorophores, ", n_detectors, " detectors...")

# Reference Matrix M (fluor x detectors)
M <- matrix(runif(n_fluor * n_detectors, 0, 1000), nrow = n_fluor, ncol = n_detectors)

# True Abundances (sparse, non-negative)
A_true <- matrix(rexp(n_cells * n_fluor), nrow = n_cells, ncol = n_fluor)
A_true[sample(length(A_true), length(A_true) * 0.5)] <- 0

# Observed Data Y = A * M + Noise
Y <- A_true %*% M + matrix(rnorm(n_cells * n_detectors, 0, 10), nrow = n_cells, ncol = n_detectors)

# Transpose M for nnls package
Mt <- t(M)

# 2. Define Old Method (Row-wise NNLS)
run_old_method <- function(M, Y) {
  if (!requireNamespace("nnls", quietly = TRUE)) {
    stop("Package 'nnls' required for benchmark.")
  }
  n_cells <- nrow(Y)
  n_fluor <- nrow(M)
  Mt <- t(M)
  A_out <- matrix(0, nrow = n_cells, ncol = n_fluor)
  for (i in seq_len(n_cells)) {
    A_out[i, ] <- nnls::nnls(Mt, Y[i, ])$x
  }
  return(A_out)
}

# 3. Define New Method (Vectorized Coordinate Descent)
solve_nnls_batch <- function(M, Y, max_iter = 50) {
  # Y: n_cells x n_detectors
  # M: n_fluor x n_detectors

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
          # Gradient direction: (A %*% XtX)[,j] - Xty[,j]
          # Newton step: A_new = A_old - Gradient / XtX[j,j]
          # Simplify: A_new = (Xty[,j] - (A %*% XtX)[,j] + A[,j]*XtX[j,j]) / XtX[j,j]

          # Note: (A %*% XtX)[, j] = A %*% XtX[, j]
          current_val <- A[, j]
          numerator <- Xty[, j] - (A %*% XtX[, j]) + current_val * XtX[j, j]
          A[, j] <- pmax(0, numerator / XtX[j, j])
      }
  }
  return(A)
}

# 4. Benchmark
message("\nStarting Benchmark...")

start_time <- Sys.time()
A_old <- run_old_method(M, Y)
time_old <- Sys.time() - start_time
message("Old Method Time: ", round(as.numeric(time_old, units = "secs"), 4), " seconds")

start_time <- Sys.time()
A_new <- solve_nnls_batch(M, Y, max_iter = 50)
time_new <- Sys.time() - start_time
message("New Method Time: ", round(as.numeric(time_new, units = "secs"), 4), " seconds")

message("\nSpeedup: ", round(as.numeric(time_old)/as.numeric(time_new), 2), "x")

# 5. Verify Correctness
# Compare Frobenius norm of difference
diff_norm <- sqrt(sum((A_old - A_new)^2)) / sqrt(sum(A_old^2))
message("\nRelative Difference (Frobenius Norm): ", formatC(diff_norm, format = "e", digits = 4))

if (diff_norm < 0.05) {
    message("SUCCESS: Results match closely.")
} else {
    message("WARNING: Results differ significantly.")
}
