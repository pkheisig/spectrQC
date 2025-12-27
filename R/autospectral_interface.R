' Create AutoSpectral Control File
#' 
#' Generates the CSV file required by the AutoSpectral package.
#' 
#' @param input_folder Directory containing single-stained control FCS files.
#' @param output_file Path where the CSV will be saved (default: "fcs_control_file.csv").
#' @param custom_fluorophores Optional named vector to map filenames to fluorophore names.
#' @return A data frame containing the control file information.
#' @export
create_autospectral_control_file <- function(input_folder = "scc", 
                                          output_file = "fcs_control_file.csv",
                                          custom_fluorophores = NULL) {
    fcs_files <- list.files(input_folder, pattern = "\\.fcs$", full.names = FALSE)
    if (length(fcs_files) == 0) stop("No FCS files found in ", input_folder)
    
    # Expanded fluorophore list from get_fluorophore_patterns
    patterns_list <- get_fluorophore_patterns()
    specific_dyes <- c(patterns_list$beads, patterns_list$cells)
    
    # Exclude generic terms from specific list to avoid high-priority matches
    generic_terms <- c("LIVE/DEAD", "Fixable Viability", "Viability", "FVD")
    specific_dyes <- setdiff(specific_dyes, generic_terms)
    
    # Sort specific dyes by length descending
    specific_dyes <- specific_dyes[order(-nchar(specific_dyes))]
    # Patterns to check in order
    fluor_patterns <- c(specific_dyes, generic_terms)
    
    df <- data.frame(
        filename = fcs_files,
        fluorophore = "Unknown",
        marker = "",
        channel = "",
        control.type = "beads",
        universal.negative = FALSE,
        large.gate = FALSE,
        is.viability = FALSE,
        stringsAsFactors = FALSE
    )
    
    for (i in seq_len(nrow(df))) {
        fn <- df$filename[i] 
        
        # Custom mapping first
        if (!is.null(custom_fluorophores) && fn %in% names(custom_fluorophores)) {
            df$fluorophore[i] <- custom_fluorophores[fn]
        } else {
            # Try patterns
            for (p in fluor_patterns) {
                if (grepl(p, fn, fixed = FALSE, ignore.case = TRUE)) {
                    df$fluorophore[i] <- p
                    break
                }
            }
        }
        
        # Check for Unstained
        if (grepl("Unstained|US_UT", fn, ignore.case = TRUE)) {
            df$fluorophore[i] <- "AF"
            df$control.type[i] <- "cells"
        }
        
        # Check for cells/beads in name
        if (grepl("cells", fn, ignore.case = TRUE)) {
            df$control.type[i] <- "cells"
        } else if (grepl("beads", fn, ignore.case = TRUE)) {
            df$control.type[i] <- "beads"
        }
        
        # Check for viability
        if (grepl("Viability|Live.Dead|7-AAD|DAPI|FVD|eFluor 506|eFluor 780", fn, ignore.case = TRUE)) {
            df$is.viability[i] <- TRUE
            df$control.type[i] <- "cells"
        }
    }
    
    # Final cleanup of names (remove duplicates if necessary, though make.unique handles it)
    df$fluorophore <- make.unique(df$fluorophore)
    
    write.csv(df, output_file, row.names = FALSE)
    message("AutoSpectral control file saved to: ", output_file)
    return(df)
}

#' Unmix using AutoSpectral
#' 
#' A wrapper function to use AutoSpectral's unmixing methods.
#' 
#' @param flow_frame A flowFrame object.
#' @param control_file Path to the AutoSpectral control file.
#' @param control_dir Directory containing control FCS files.
#' @param method Unmixing method ("WLS", "OLS", or "Poisson").
#' @param cytometer Cytometer type (e.g., "Aurora", "ID7000").
#' @export
unmix_with_autospectral <- function(flow_frame, 
                                   control_file = "fcs_control_file.csv", 
                                   control_dir = "scc",
                                   method = "WLS",
                                   cytometer = "Aurora") {
    if (!requireNamespace("AutoSpectral", quietly = TRUE)) {
        stop("Package 'AutoSpectral' required. Install from GitHub.")
    }
    
    # 1. Setup AutoSpectral parameters
    if (cytometer == "Aurora") {
        asp <- AutoSpectral::get.autospectral.param.aurora()
    } else if (cytometer == "ID7000") {
        asp <- AutoSpectral::get.autospectral.param.id7000()
    } else {
        asp <- AutoSpectral::get.autospectral.param.minimal()
    }
    
    # 2. Define flow control
    flow.control <- AutoSpectral::define.flow.control(control.dir = control_dir, 
                                                    control.def.file = control_file, 
                                                    asp = asp)
    
    # 3. Get spectra
    spectra <- AutoSpectral::get.fluorophore.spectra(flow.control, asp)
    
    # 4. Prepare data
    raw_data <- flowCore::exprs(flow_frame)
    detectors <- colnames(spectra)
    
    # Check if detectors match
    missing <- setdiff(detectors, colnames(raw_data))
    if (length(missing) > 0) {
        stop("Detectors in AutoSpectral spectra not found in flow_frame: ", paste(missing, collapse = ", "))
    }
    
    Y <- raw_data[, detectors, drop = FALSE]
    
    # 5. Unmix
    method <- toupper(method)
    if (method == "WLS") {
        unmixed <- AutoSpectral::unmix.wls(Y, spectra)
    } else if (method == "OLS") {
        unmixed <- AutoSpectral::unmix.ols(Y, spectra)
    } else if (method == "POISSON") {
        unmixed <- AutoSpectral::unmix.poisson(Y, spectra, asp = asp)
    } else {
        stop("Unsupported AutoSpectral method: ", method)
    }
    
    return(unmixed)
}
