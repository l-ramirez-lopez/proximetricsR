# Tests for proxiscout_write_data()
#
# The function writes up to three comma-separated files into `path`:
#   * <prefix>_spectra.csv    -- always; sample name, device id and spectra
#   * <prefix>_properties.csv -- only when `properties` are given
#   * <prefix>_metadata.csv   -- only when columns remain that are written
#                                nowhere else
# The tests below exercise all three files, the automatic column detection /
# renaming, and the error paths.

# -----------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------

# A small, hand-made proxiscout-style data.frame. The spectra live in a single
# matrix column "spc" (as produced by proxiscout_read_data()) with numeric
# wavelength headers; sample and device columns are detected by name; and there
# are extra columns ("operator", "batch") that belong in the metadata file.
make_pcs_data <- function() {
  df <- data.frame(
    sampleName = c("A_1", "A_2", "B_1"),
    deviceId = c("dev1", "dev1", "dev2"),
    protein = c(10.5, 10.5, 20.1),
    operator = c("alice", "alice", "bob"),
    batch = c("b1", "b1", "b2"),
    stringsAsFactors = FALSE
  )
  spc <- matrix(seq_len(3 * 4) / 100, nrow = 3)
  colnames(spc) <- c("1000", "1100", "1200", "1300")
  df$spc <- spc
  df
}

# read one of the written csv files back, preserving numeric headers
read_out <- function(path) {
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

# -----------------------------------------------------------------------
# Return value and file creation
# -----------------------------------------------------------------------

test_that("all three file paths are returned and exist when applicable", {
  out <- withr::local_tempdir()
  paths <- proxiscout_write_data(make_pcs_data(), path = out, properties = "protein")

  expect_type(paths, "character")
  expect_length(paths, 3L)
  expect_true(all(file.exists(paths)))
  expect_setequal(
    basename(paths),
    c(
      "proxiscout_export_spectra.csv",
      "proxiscout_export_properties.csv",
      "proxiscout_export_metadata.csv"
    )
  )
})

test_that("the return value is invisible", {
  out <- withr::local_tempdir()
  expect_output(proxiscout_write_data(make_pcs_data(), path = out), NA)
})

test_that("file_prefix is honoured for every file", {
  out <- withr::local_tempdir()
  proxiscout_write_data(make_pcs_data(), path = out, file_prefix = "myrun", properties = "protein")

  expect_true(file.exists(file.path(out, "myrun_spectra.csv")))
  expect_true(file.exists(file.path(out, "myrun_properties.csv")))
  expect_true(file.exists(file.path(out, "myrun_metadata.csv")))
})

# -----------------------------------------------------------------------
# Spectra file
# -----------------------------------------------------------------------

test_that("spectra file holds sample name, device id and spectra scaled by 100", {
  out <- withr::local_tempdir()
  x <- make_pcs_data()
  proxiscout_write_data(x, path = out, properties = "protein")

  spectra <- read_out(file.path(out, "proxiscout_export_spectra.csv"))

  expect_equal(names(spectra), c("sampleName", "deviceId", "1000", "1100", "1200", "1300"))
  expect_equal(spectra[["sampleName"]], x$sampleName)
  expect_equal(spectra[["deviceId"]], x$deviceId)
  # spectra are multiplied by 100 on write
  expect_equal(unname(as.matrix(spectra[, c("1000", "1100", "1200", "1300")])), unname(x$spc * 100))
})

test_that("spectra file never contains metadata or property columns", {
  out <- withr::local_tempdir()
  proxiscout_write_data(make_pcs_data(), path = out, properties = "protein")

  spectra <- read_out(file.path(out, "proxiscout_export_spectra.csv"))
  expect_false(any(c("operator", "batch", "protein") %in% names(spectra)))
})

# -----------------------------------------------------------------------
# Properties file
# -----------------------------------------------------------------------

test_that("no properties file is written when properties is NULL or empty", {
  out <- withr::local_tempdir()

  paths_null <- proxiscout_write_data(make_pcs_data(), path = out)
  expect_false(file.exists(file.path(out, "proxiscout_export_properties.csv")))
  expect_false(any(grepl("_properties\\.csv$", paths_null)))

  out2 <- withr::local_tempdir()
  proxiscout_write_data(make_pcs_data(), path = out2, properties = character(0))
  expect_false(file.exists(file.path(out2, "proxiscout_export_properties.csv")))
})

test_that("properties file strips repetition suffixes and de-duplicates rows", {
  out <- withr::local_tempdir()
  proxiscout_write_data(make_pcs_data(), path = out, properties = "protein")

  props <- read_out(file.path(out, "proxiscout_export_properties.csv"))

  expect_equal(names(props), c("sampleName", "protein"))
  # "A_1"/"A_2" collapse to a single "A" row; "B_1" -> "B"
  expect_equal(props$sampleName, c("A", "B"))
  expect_equal(props$protein, c(10.5, 20.1))
})

test_that("properties file drops rows where all properties are NA", {
  out <- withr::local_tempdir()
  x <- data.frame(
    sampleName = c("A_1", "C_2", "B_1"),
    deviceId = "dev1",
    p1 = c(1, NA, 3),
    p2 = c(NA, NA, 4),
    stringsAsFactors = FALSE
  )
  x$spc <- matrix(seq_len(3 * 4) / 100, nrow = 3)

  proxiscout_write_data(x, path = out, properties = c("p1", "p2"))
  props <- read_out(file.path(out, "proxiscout_export_properties.csv"))

  # the all-NA row (C_2) is dropped; A_1 and B_1 remain
  expect_equal(nrow(props), 2L)
  expect_equal(props$sampleName, c("A", "B"))
  # NA values are written as empty strings, so a partially-NA cell reads as NA
  expect_true(is.na(props$p2[props$sampleName == "A"]))
})

test_that("an error is raised when a requested property is missing from x", {
  out <- withr::local_tempdir()
  expect_error(
    proxiscout_write_data(make_pcs_data(), path = out, properties = c("protein", "nope")),
    "Properties not found in 'x': nope"
  )
})

# -----------------------------------------------------------------------
# Metadata file
# -----------------------------------------------------------------------

test_that("metadata holds sampleName plus only the leftover columns", {
  out <- withr::local_tempdir()
  proxiscout_write_data(make_pcs_data(), path = out, properties = "protein")

  meta <- read_out(file.path(out, "proxiscout_export_metadata.csv"))

  expect_equal(names(meta), c("sampleName", "operator", "batch"))
  # spectra, device id and the property column live in the other files
  expect_false("deviceId" %in% names(meta))
  expect_false("protein" %in% names(meta))
  expect_false(any(names(meta) %in% c("1000", "1100", "1200", "1300")))
})

test_that("metadata keeps one row per input row, sample names verbatim", {
  out <- withr::local_tempdir()
  x <- make_pcs_data()
  proxiscout_write_data(x, path = out, properties = "protein")

  meta <- read_out(file.path(out, "proxiscout_export_metadata.csv"))
  expect_equal(nrow(meta), nrow(x))
  # unlike the properties file, the repetition suffix is preserved and rows are
  # not de-duplicated
  expect_equal(meta$sampleName, x$sampleName)
  expect_equal(meta$operator, x$operator)
})

test_that("without properties the metadata absorbs the leftover property column", {
  out <- withr::local_tempdir()
  proxiscout_write_data(make_pcs_data(), path = out) # properties = NULL

  meta <- read_out(file.path(out, "proxiscout_export_metadata.csv"))
  # protein is no longer written to a properties file, so it becomes metadata
  expect_equal(names(meta), c("sampleName", "protein", "operator", "batch"))
})

test_that("no metadata file is written when no leftover columns remain", {
  out <- withr::local_tempdir()
  x <- data.frame(
    sampleName = c("A_1", "A_2"),
    deviceId = c("dev1", "dev1"),
    protein = c(1.5, 2.5),
    stringsAsFactors = FALSE
  )
  x$spc <- matrix(seq_len(2 * 4) / 100, nrow = 2)

  paths <- proxiscout_write_data(x, path = out, properties = "protein")

  # only spectra + properties -> no metadata (consistent with the properties file)
  expect_length(paths, 2L)
  expect_false(file.exists(file.path(out, "proxiscout_export_metadata.csv")))
  expect_false(any(grepl("_metadata\\.csv$", paths)))
})

test_that("metadata respects a custom spc column selection by index", {
  out <- withr::local_tempdir()
  # spectra given as a plain block of numeric columns selected via `spc`
  x <- data.frame(
    sampleName = c("A_1", "B_1"),
    deviceId = c("dev1", "dev2"),
    note = c("x", "y"),
    "1000" = c(0.1, 0.2),
    "1100" = c(0.3, 0.4),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  proxiscout_write_data(x, path = out, spc = 4:5)
  meta <- read_out(file.path(out, "proxiscout_export_metadata.csv"))

  # the spectra columns (indices 4:5) are excluded, "note" remains
  expect_equal(names(meta), c("sampleName", "note"))
  expect_false(any(names(meta) %in% c("1000", "1100")))
})

# -----------------------------------------------------------------------
# Automatic column detection and renaming
# -----------------------------------------------------------------------

test_that("an 'ID' column is used as the sample name when no sample column exists", {
  out <- withr::local_tempdir()
  x <- data.frame(
    ID = c("A_1", "B_1"),
    deviceId = c("dev1", "dev2"),
    stringsAsFactors = FALSE
  )
  x$spc <- matrix(seq_len(2 * 4) / 100, nrow = 2)

  proxiscout_write_data(x, path = out)
  spectra <- read_out(file.path(out, "proxiscout_export_spectra.csv"))
  expect_equal(spectra[["sampleName"]], x$ID)
})

test_that("an 'SNR'/'SRN' column is used as the device id when no device column exists", {
  out <- withr::local_tempdir()
  x <- data.frame(
    sampleName = c("A_1", "B_1"),
    SNR = c("dev1", "dev2"),
    stringsAsFactors = FALSE
  )
  x$spc <- matrix(seq_len(2 * 4) / 100, nrow = 2)

  proxiscout_write_data(x, path = out)
  spectra <- read_out(file.path(out, "proxiscout_export_spectra.csv"))
  expect_equal(spectra[["deviceId"]], x$SNR)
})

test_that("a 'scanner' column is accepted as the device column", {
  out <- withr::local_tempdir()
  x <- data.frame(
    sampleName = c("A_1", "B_1"),
    scannerId = c("s1", "s2"),
    stringsAsFactors = FALSE
  )
  x$spc <- matrix(seq_len(2 * 4) / 100, nrow = 2)

  expect_no_error(proxiscout_write_data(x, path = out))
  spectra <- read_out(file.path(out, "proxiscout_export_spectra.csv"))
  expect_equal(spectra[["deviceId"]], x$scannerId)
})

# -----------------------------------------------------------------------
# Error paths
# -----------------------------------------------------------------------

test_that("an error is raised when no sample column can be detected", {
  out <- withr::local_tempdir()
  x <- data.frame(deviceId = c("dev1", "dev2"), stringsAsFactors = FALSE)
  x$spc <- matrix(seq_len(2 * 4) / 100, nrow = 2)

  expect_error(proxiscout_write_data(x, path = out), "No sample column detected")
})

test_that("an error is raised when no device or scanner column can be detected", {
  out <- withr::local_tempdir()
  x <- data.frame(sampleName = c("A_1", "B_1"), stringsAsFactors = FALSE)
  x$spc <- matrix(seq_len(2 * 4) / 100, nrow = 2)

  expect_error(proxiscout_write_data(x, path = out), "No device or scanner column detected")
})
