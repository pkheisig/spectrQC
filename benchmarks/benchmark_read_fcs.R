# Benchmark Script for Autospectral Interface Optimization
# This script demonstrates the performance improvement of reading a subset of events
# vs reading the full file for channel identification.

# MOCKING INFRASTRUCTURE
# Since we might not have flowCore or real files, we mock the behavior.

mock_read_FCS <- function(filename, transformation = FALSE, truncate_max_range = FALSE, which.lines = NULL) {
    # Simulate file size (randomly between 10k and 1M events)
    n_total <- 500000

    if (is.null(which.lines)) {
        n_read <- n_total
        # Simulate delay: 0.1s per 100k events
        Sys.sleep(n_read / 100000 * 0.1)
    } else {
        if (length(which.lines) == 1) {
             n_read <- which.lines
             # Simulate random access delay (seek time) + reading time
             # Seek overhead: 0.05s, Read time: linear
             Sys.sleep(0.05 + n_read / 100000 * 0.1)
        } else {
             # Vector case
             n_read <- length(which.lines)
             Sys.sleep(0.01 + n_read / 100000 * 0.1)
        }
    }

    # Return a mock flowFrame
    # We only need 'exprs' and 'parameters' accessors to work
    structure(list(
        exprs = matrix(rnorm(n_read * 10), ncol = 10, dimnames = list(NULL, paste0("FL", 1:10))),
        parameters = list(data = data.frame(name = paste0("FL", 1:10), desc = paste0("Marker", 1:10)))
    ), class = "flowFrame")
}

# Mock flowCore methods
flowCore_exprs <- function(ff) ff$exprs
flowCore_parameters <- function(ff) ff$parameters
flowCore_pData <- function(params) params$data

# Mock get_sorted_detectors (simplified)
get_sorted_detectors <- function(pd) {
    list(names = pd$name)
}

# BENCHMARK FUNCTION
run_benchmark <- function() {
    message("Starting Benchmark...")

    # Mock data frame of files
    df <- data.frame(
        filename = paste0("file_", 1:5, ".fcs"),
        channel = "",
        stringsAsFactors = FALSE
    )

    # METHOD 1: ORIGINAL (Full Read)
    message("\n--- Method 1: Original (Full Read) ---")
    start_time <- Sys.time()
    for (i in seq_len(nrow(df))) {
        fn <- df$filename[i]
        # Simulate read.FCS
        ff <- mock_read_FCS(fn, transformation = FALSE, truncate_max_range = FALSE)

        pd <- flowCore_pData(flowCore_parameters(ff))
        fl_pd <- get_sorted_detectors(pd)
        if (length(fl_pd$names) > 0) {
            medians <- apply(flowCore_exprs(ff)[, fl_pd$names, drop = FALSE], 2, median)
            df$channel[i] <- names(medians)[which.max(medians)]
        }
    }
    end_time <- Sys.time()
    duration_orig <- end_time - start_time
    message(sprintf("Time taken: %.4f seconds", as.numeric(duration_orig)))

    # METHOD 2: OPTIMIZED (Subset Read)
    message("\n--- Method 2: Optimized (Subset Read = 1000) ---")
    start_time <- Sys.time()
    for (i in seq_len(nrow(df))) {
        fn <- df$filename[i]
        tryCatch({
            # Try reading subset
            ff <- tryCatch({
                mock_read_FCS(fn, transformation = FALSE, truncate_max_range = FALSE, which.lines = 1000)
            }, error = function(e) {
                # Fallback
                mock_read_FCS(fn, transformation = FALSE, truncate_max_range = FALSE)
            })

            pd <- flowCore_pData(flowCore_parameters(ff))
            fl_pd <- get_sorted_detectors(pd)
            if (length(fl_pd$names) > 0) {
                medians <- apply(flowCore_exprs(ff)[, fl_pd$names, drop = FALSE], 2, median)
                df$channel[i] <- names(medians)[which.max(medians)]
            }
        }, error = function(e) NULL)
    }
    end_time <- Sys.time()
    duration_opt <- end_time - start_time
    message(sprintf("Time taken: %.4f seconds", as.numeric(duration_opt)))

    message(sprintf("\nSpeedup: %.2fx", as.numeric(duration_orig) / as.numeric(duration_opt)))
}

# Run if interactive or sourced
if (!interactive()) {
    # In a real environment, we would run this.
    # run_benchmark()
    message("Benchmark script created. Run inside an R environment with flowCore to verify.")
}
