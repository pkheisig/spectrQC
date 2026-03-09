testthat::test_that("create_control_file recognizes fluor and control type from filename variants", {
    testthat::skip_if_not_installed("spectreasy")

    scc_dir <- tempfile("spectreasy_scc_")
    dir.create(scc_dir, recursive = TRUE, showWarnings = FALSE)

    files <- c(
        "LIVE DEAD NIR (Cells).fcs",
        "PE-CF594 (Beads).fcs",
        "PECF594 (beads).fcs",
        "pe cy7 (bEaDs).fcs",
        "PE-Fire 700 (BEADS).fcs",
        "PE (Beads).fcs"
    )
    created <- file.create(file.path(scc_dir, files))
    testthat::expect_true(all(created))

    out_csv <- tempfile(fileext = ".csv")
    df <- spectreasy::create_control_file(
        input_folder = scc_dir,
        include_af_folder = FALSE,
        output_file = out_csv
    )

    by_file <- split(df, df$filename)

    testthat::expect_equal(by_file[["LIVE DEAD NIR (Cells).fcs"]]$fluorophore[[1]], "LIVE/DEAD NIR")
    testthat::expect_equal(by_file[["LIVE DEAD NIR (Cells).fcs"]]$control.type[[1]], "cells")
    testthat::expect_equal(by_file[["LIVE DEAD NIR (Cells).fcs"]]$is.viability[[1]], "TRUE")

    testthat::expect_equal(by_file[["PE-CF594 (Beads).fcs"]]$fluorophore[[1]], "PE-CF594")
    testthat::expect_equal(by_file[["PE-CF594 (Beads).fcs"]]$control.type[[1]], "beads")

    testthat::expect_equal(by_file[["PECF594 (beads).fcs"]]$fluorophore[[1]], "PE-CF594")
    testthat::expect_equal(by_file[["PECF594 (beads).fcs"]]$control.type[[1]], "beads")

    testthat::expect_equal(by_file[["pe cy7 (bEaDs).fcs"]]$fluorophore[[1]], "PE-Cy7")
    testthat::expect_equal(by_file[["pe cy7 (bEaDs).fcs"]]$control.type[[1]], "beads")

    testthat::expect_equal(by_file[["PE-Fire 700 (BEADS).fcs"]]$fluorophore[[1]], "PE-Fire 700")
    testthat::expect_equal(by_file[["PE-Fire 700 (BEADS).fcs"]]$control.type[[1]], "beads")

    testthat::expect_equal(by_file[["PE (Beads).fcs"]]$fluorophore[[1]], "PE")
    testthat::expect_equal(by_file[["PE (Beads).fcs"]]$control.type[[1]], "beads")
    testthat::expect_equal(by_file[["PE (Beads).fcs"]]$is.viability[[1]], "")
})
