#' @docType data
#' @name proxiscoutCannabis
#' @aliases proxiscoutCannabis
#' @title proxiscoutCannabis
#' @format A `data.frame` of class \code{"proxiscout_data"} containing 192
#' observations of eleven response variables, with their corresponding spectral
#' data.
#' @usage
#' data("proxiscoutCannabis")
#' @description
#' Selected samples of cannabis NIR measurements for demo purposes. Note these 
#' samples are different from the ones in \code{\link{proximateCannabis}}.
#' The dataset contains absorbance spectra of 192 cannabis samples measured with
#' three ProxiScout devices at the 257 standard ProxiScout wavenumbers, ranging
#' from 3921.569 \eqn{cm^{-1}} to 7407.407 \eqn{cm^{-1}} in steps of around
#' 13.61655 \eqn{cm^{-1}} (equivalent to a spectral range of 1350 nm to 2550 nm),
#' see \code{\link{get_proxiscout_wavenumbers}}. A total number of eleven
#' reference vectors is included: \code{"CBDA"} (Cannabidiolic acid),
#' \code{"CBG"} (Cannabigerol), \code{"CBD"} (Cannabidiol), \code{"CBN"}
#' (Cannabinol), \code{"THC"} (Tetrahydrocannabinol), \code{"THCA"}
#' (Tetrahydrocannabinolic acid), \code{"CBC"} (Cannabichromene), \code{"CBGA"}
#' (Cannabigerolic acid), \code{"CBD_total"} (total Cannabidiol),
#' \code{"THC_total"} (total Tetrahydrocannabinol) and \code{"CBGt"} (total
#' Cannabigerol).
#' @details
#' This dataset is an example for a typical data file for ProxiScout
#' applications, with a total of 192 cannabis samples, selected as a subset of a
#' larger database. It was obtained by reading a spectra file together with its
#' corresponding property file using \code{\link{proxiscout_read_data}}. As the
#' data contains 257 spectral columns, the object is of class
#' \code{"proxiscout_data"}.
#'
#' Since the samples were measured with three different devices, the dataset is
#' also suitable for illustrating the grouped validation plots of
#' \code{\link{calibrate_models}} and \code{\link{plot.spectral_model}}, using
#' \code{deviceId} as the grouping variable.
#'
#' The dataset contains one clear outlier, a single sample whose absorbance is
#' markedly higher than that of the rest of the set across the whole spectral
#' range. This is a legitimate measurement and has been retained intentionally,
#' as it is useful for demonstrating the detection of outliers, for example with
#' the leverage and spectral residual diagnostics returned by
#' \code{\link{validate_prediction}}.
#'
#' It contains the following rows for each observation:
#' \itemize{
#'     \item \strong{\code{sampleName}:} Characters for the name of the sample.
#'     Repeated scans of the same sample are indicated by a numerical suffix, see
#'     \code{\link{proxiscout_repetition_pattern}}.
#'     \item \strong{\code{deviceId}:} Characters of the identifier of the
#'     involved ProxiScout device.
#'     \item \strong{\code{train}:} Logicals, indicating whether the particular
#'     observation belongs to the training set used for the construction of the
#'     model.
#'     \item \strong{\code{type}:} Characters for the type of the sample.
#'     \item \strong{\code{CBDA}:} Numerics for the reference values of
#'     Cannabidiolic acid.
#'     \item \strong{\code{CBG}:} Numerics for the reference values of
#'     Cannabigerol.
#'     \item \strong{\code{CBD}:} Numerics for the reference values of
#'     Cannabidiol.
#'     \item \strong{\code{CBN}:} Numerics for the reference values of
#'     Cannabinol.
#'     \item \strong{\code{THC}:} Numerics for the reference values of
#'     Tetrahydrocannabinol.
#'     \item \strong{\code{THCA}:} Numerics for the reference values of
#'     Tetrahydrocannabinolic acid.
#'     \item \strong{\code{CBC}:} Numerics for the reference values of
#'     Cannabichromene.
#'     \item \strong{\code{CBGA}:} Numerics for the reference values of
#'     Cannabigerolic acid.
#'     \item \strong{\code{CBD_total}:} Numerics for the reference values of total
#'     Cannabidiol.
#'     \item \strong{\code{THC_total}:} Numerics for the reference values of total
#'     Tetrahydrocannabinol.
#'     \item \strong{\code{CBGt}:} Numerics for the reference values of total
#'     Cannabigerol.
#'     \item \strong{\code{.repetition_group}:} Integers identifying the rows that
#'     correspond to repeated scans of the same sample, added by
#'     \code{\link{proxiscout_read_data}}.
#'     \item \strong{\code{spc}:} A numerical matrix of the absorbance spectra,
#'     corresponding to each individual observation.
#' }
#' @author Marçal Plans
#' @seealso \code{\link{proximateCannabis}}
#' @source BUCHI Labortechnik AG.
#' @keywords datasets
NULL
