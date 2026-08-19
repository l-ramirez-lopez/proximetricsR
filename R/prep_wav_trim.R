#' @title Wavelength trimming constructor for spectral preprocessing
#'
#' @description
#'
#' \loadmathjax
#'
#' Creates a preprocessing constructor for trimming spectral data to a
#' specified wavelength band. The constructor is intended to be passed to
#' \code{\link{preprocess_recipe}} and executed via \code{\link{process}}.
#'
#' @usage
#' prep_wav_trim(
#'   band,
#'   trim_constant_edges = FALSE
#' )
#'
#' @param band A numeric vector of length 2 giving the minimum and maximum
#' wavenumber/wavelength to retain. Columns of \code{X} outside this range are dropped.
#' Pass \code{c()} (empty vector) to skip band trimming and only apply
#' \code{trim_constant_edges}.
#' @param trim_constant_edges A logical. If \code{TRUE}, constant or zero-valued
#' columns at the left and right edges are removed after band trimming. Default
#' is \code{FALSE}.
#'
#' @return An object of class \code{preprocessing} to be used in
#' \code{\link{preprocess_recipe}} and executed by \code{\link{process}}.
#'
#' @details
#' Band trimming retains only those columns whose names (coerced to numeric)
#' fall within \code{[min(band), max(band)]}. If no columns fall within the
#' band the original matrix is returned with a warning.
#'
#' Constant edge trimming scans inward from each edge and drops columns that
#' are identical to their immediate neighbour or are all zero. If trimming
#' would leave fewer than two columns the step is skipped with a warning.
#'
#' Because constant edge trimming depends on the data values, it is resolved to
#' the exact set of wavelengths retained on the training data when the step is
#' used in \code{\link{calibrate}}, so that predictions trim \code{newdata} to
#' exactly the same columns.
#'
#' @author Claudio Orellano and Leonardo Ramirez-Lopez
#'
#' @examples
#' data("NIRcannabis")
#' X <- NIRcannabis$spc
#'
#' tr <- prep_wav_trim(band = c(1000, 1800))
#' recipe <- preprocess_recipe(tr, device = "proxiscout")
#' X_trim <- process(X, recipe)
#'
#' @seealso \code{\link{preprocess_recipe}}, \code{\link{process}}
#' @export
prep_wav_trim <- function(band, trim_constant_edges = FALSE) {
  if (missing(band)) {
    stop("'band' is required. Use c() to skip band trimming.")
  }

  if (length(band) != 2 && length(band) != 0) {
    stop("'band' must be of length 0 (no band trimming) or 2 (min and max wavelength).")
  }
  if (any(is.na(band))) {
    stop("values in 'band' cannot be NA.")
  }
  if (length(band) == 2 && min(band) >= max(band)) {
    stop("'band[1]' must be strictly less than 'band[2]'.")
  }
  if (!is.logical(trim_constant_edges) || length(trim_constant_edges) != 1 || is.na(trim_constant_edges)) {
    stop("'trim_constant_edges' must be a single logical value.")
  }

  structure(
    list(
      method = "prep_wav_trim",
      band = band,
      trim_constant_edges = trim_constant_edges,
      compatible_devices = "proxiscout"
    ),
    class = c("preprocessing", "list")
  )
}

#' @keywords internal
.exec_wav_trim <- function(X, step) {
  X_trim <- X

  # Fitted step: select the exact wavelengths retained during calibration (see
  # .freeze_trim_steps), in training order, and skip the constant-edge scan.
  if (!is.null(step$resolved_wavs)) {
    wav <- as.numeric(colnames(X_trim))
    if (any(is.na(wav))) {
      warning("Column names are not numeric wavelengths; band trimming skipped.")
      return(X_trim)
    }
    idx <- match(step$resolved_wavs, wav)
    if (anyNA(idx)) {
      miss <- step$resolved_wavs[is.na(idx)]
      stop("'newdata' is missing ", length(miss), " wavelength(s) required by the model.")
    }
    return(X_trim[, idx, drop = FALSE])
  }

  if (length(step$band) == 2) {
    wav <- as.numeric(colnames(X))
    if (any(is.na(wav))) {
      warning("Column names are not numeric wavelengths; band trimming skipped.")
    } else {
      in_range <- which(wav >= min(step$band) & wav <= max(step$band))
      if (length(in_range) < 1) {
        warning("Band trimming would drop all columns; step ignored.")
      } else {
        X_trim <- X_trim[, in_range, drop = FALSE]
      }
    }
  }

  if (step$trim_constant_edges && ncol(X_trim) > 3) {
    left <- 1
    while (
      left < ncol(X_trim) &&
        (isTRUE(all(X_trim[, left] == X_trim[, left + 1])) || isTRUE(all(X_trim[, left] == 0)))
    ) {
      left <- left + 1
    }

    right <- ncol(X_trim)
    while (
      right > left &&
        (isTRUE(all(X_trim[, right] == X_trim[, right - 1])) || isTRUE(all(X_trim[, right] == 0)))
    ) {
      right <- right - 1
    }

    if (right - left < 1) {
      warning("Constant edge trimming would leave fewer than two columns; step ignored.")
    } else {
      X_trim <- X_trim[, left:right, drop = FALSE]
    }
  }
  X_trim
}

#' Freeze data-dependent constant-edge trim steps into a fixed wavelength band
#'
#' Records, on each \code{prep_wav_trim} step using
#' \code{trim_constant_edges = TRUE}, the exact wavelengths retained on the
#' training data (\code{resolved_wavs}) so \code{.exec_wav_trim} reapplies them
#' deterministically to \code{newdata} instead of re-deriving the trimming from
#' its own values.
#'
#' @param recipe A \code{preprocess_recipe} object.
#' @param processed_wavs A \code{processed_wavs} object holding the wavelengths
#' retained after each step (from the \code{"processed_wavs"} attribute of
#' \code{\link{process}}).
#'
#' @return The \code{recipe} with \code{resolved_wavs} set on each fitted step.
#'
#' @keywords internal
.freeze_trim_steps <- function(recipe, processed_wavs) {
  for (j in seq_along(recipe$steps)) {
    step <- recipe$steps[[j]]
    if (isTRUE(step$method == "prep_wav_trim") && isTRUE(step$trim_constant_edges)) {
      retained <- processed_wavs[[paste0("step_", j)]]
      if (length(retained) > 0) {
        recipe$steps[[j]]$resolved_wavs <- retained
      }
    }
  }
  recipe
}
