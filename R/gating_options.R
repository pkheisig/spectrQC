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
