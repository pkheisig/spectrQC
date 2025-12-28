#' Gating options
#' 
#' @param histogram_pct_beads Percentage of beads to gate in histogram
#' @param histogram_direction_beads Direction to gate beads ("both", "left", "right")
#' @param histogram_pct_cells Percentage of cells to gate in histogram
#' @param histogram_direction_cells Direction to gate cells ("both", "left", "right")
#' @export
gating_options <- function(histogram_pct_beads = 0.98,
                          histogram_direction_beads = "both",
                          histogram_pct_cells = 0.35,
                          histogram_direction_cells = "both") {
    list(
        histogram_pct_beads = histogram_pct_beads,
        histogram_direction_beads = histogram_direction_beads,
        histogram_pct_cells = histogram_pct_cells,
        histogram_direction_cells = histogram_direction_cells
    )
}

#' Get standard fluorophore patterns
#' @export
get_fluorophore_patterns <- function() {
    list(
        unstained = c("US_UT", "Unstained", "unstained", "blank", "AF only", "AF"),
        beads = c(
            # Alexa dyes
            "Alexa 350", "Alexa 405", "Alexa 430", "Alexa 488", "Alexa 514",
            "Alexa 532", "Alexa 546", "Alexa 555", "Alexa 568", "Alexa 594",
            "Alexa 610", "Alexa 633", "Alexa 647", "Alexa 660", "Alexa 680",
            "Alexa 700", "Alexa 750", "Alexa 790",
            # Proteins
            "FITC", "GFP", "YFP", "CFP", "mCherry", "DsRed", "tdTomato",
            "PE-Cy5.5", "PE-Cy5", "PE-Cy7", "PE-CF594", "PE-Texas Red",
            "PE-Dazzle 594", "PE-Fire 640", "PE-Fire 700", "PE-Fire 810", "PE",
            "APC-Cy7", "APC-H7", "APC-R700", "APC-Fire 750", "APC-Fire 810", "APC",
            "PerCP-Cy5.5", "PerCP-eFluor 710", "PerCP",
            # Brilliant Violet
            "BV421", "BV480", "BV510", "BV570", "BV605", "BV650", "BV711", "BV750", "BV785",
            # Brilliant UltraViolet
            "BUV395", "BUV496", "BUV563", "BUV615", "BUV661", "BUV737", "BUV805",
            # Others
            "Pacific Blue", "Pacific Orange", "Cy5.5", "Cy5", "Cy7", "Cy3",
            "eFluor 450", "eFluor 520", "eFluor 660", "V450", "V500",
            "SB436", "SB600", "SB645", "SB702", "SB780"
        ),
        cells = c(
            "DAPI", "Hoechst", "PI", "7-AAD", "FVD",
            "eFluor 506", "eFluor 780", "eFluor 455UV",
            "Zombie Aqua", "Zombie NIR", "Zombie UV", "Zombie Violet",
            "Zombie Green", "Zombie Red", "Zombie Yellow",
            "LIVE/DEAD", "Fixable Viability", "Viability", "Calcein", "SYTOX"
        )
    )
}

#' Get Sorted Detectors by Laser and Wavelength
#' 
#' @param pd pData from flowFrame parameters
#' @return A list with [[names]] (FL...) and [[labels]] (405nm - 420/10) sorted by laser.
#' @export
get_sorted_detectors <- function(pd) {
    # 1. Identify spectral detectors
    # Common patterns: FL[0-9]+-A, [A-Z][0-9]+-A, etc.
    # Exclude FSC, SSC, Time
    exclude_patterns <- "^FSC|^SSC|^Time|^Event|^ID"
    matches <- grep(exclude_patterns, pd$name, ignore.case = TRUE, invert = TRUE)
    
    # Further filter for Area channels if they exist, otherwise use all
    area_matches <- grep("-A$", pd$name[matches])
    if (length(area_matches) > 0) {
        matches <- matches[area_matches]
    }
    
    fl_pd <- pd[matches, ]
    
    # 2. Extract descriptions/labels
    # Use 'desc' if it looks like a filter name, otherwise fallback to 'name'
    desc <- as.character(fl_pd$desc)
    names <- as.character(fl_pd$name)
    
    # Final labels for plotting
    labels <- ifelse(!is.na(desc) & desc != "" & desc != names, desc, names)
    
    # 3. Parse laser and wavelength for sorting
    # Try to find laser nm (e.g., "405nm" or "405-")
    laser_nm <- as.integer(gsub(".*?([0-9]{3})\\s*nm.*", "\\1", labels))
    # If failed, try name (e.g., V1, B2)
    if (all(is.na(laser_nm))) {
        laser_code <- substr(names, 1, 1)
        laser_nm <- ifelse(laser_code == "U", 355,
                    ifelse(laser_code == "V", 405,
                    ifelse(laser_code == "B", 488,
                    ifelse(laser_code == "Y" | laser_code == "G", 561,
                    ifelse(laser_code == "R", 640, 999)))))
    }
    
    # Laser order: UV (~355), V (~405), B (~488), YG (~561), R (~640)
    laser_priority <- ifelse(laser_nm < 360, 1, # UV
                      ifelse(laser_nm < 420, 2, # V
                      ifelse(laser_nm < 500, 3, # B
                      ifelse(laser_nm < 600, 4, # YG
                      5)))) # R
    
    # Wavelength: look for numbers after the laser name
    wavelength <- as.integer(gsub(".*?([0-9]{3}).*", "\\1", sub("^[0-9]{3}nm", "", labels)))
    if (all(is.na(wavelength))) wavelength <- seq_along(names) # Fallback to index
    
    # 4. Sort
    ord <- order(laser_priority, wavelength, na.last = TRUE)
    
    return(list(
        names = names[ord],
        labels = labels[ord],
        laser_nm = laser_nm[ord]
    ))
}
