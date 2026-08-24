#' @title Computes the Mahalanobis (GH) distance from scaled PLS scores
#' @description
#' Computes the squared Mahalanobis distance of each row of \code{x} to the
#' origin, using the covariance matrix estimated from \code{reference},
#' normalized by the number of components used. The origin represents the
#' center of the calibration score space. For calibration statistics,
#' \code{x} and \code{reference} are the same (calibration) scaled scores. For
#' new/prediction samples, \code{reference} must be the calibration model's
#' scaled scores, so that the covariance is estimated from the calibration set
#' and merely applied to (possibly few) prediction samples, rather than
#' re-estimated from them.
#' @param x a matrix of scores scaled by their calibration standard
#' deviations, for which the distance is computed.
#' @param reference a matrix of scores scaled by their calibration standard
#' deviations, used to estimate the covariance matrix. Typically the
#' calibration model's scaled scores.
#' @param ncomp a vector of component counts to compute the distance for.
#' @return A matrix with one row per sample and one column per entry in \code{ncomp}.
#' @keywords internal
.gh_distance <- function(x, reference, ncomp) {
  get_diss <- function(factors) {
    ref <- reference[, 1:factors, drop = FALSE]
    obs <- x[, 1:factors, drop = FALSE]
    stats::mahalanobis(obs, center = rep(0, factors), cov = cov(ref))
  }
  
  diss <- matrix(
    unlist(lapply(ncomp, get_diss), use.names = FALSE),
    nrow = nrow(x), ncol = length(ncomp)
  )
  sweep(diss, MARGIN = 2, STATS = ncomp, FUN = "/")
}

#' @title Leverage (Mahalanobis) limit for a new observation
#' @description
#' Upper control limit for the \code{ncomp}-normalized Mahalanobis leverage
#' produced by \code{\link{.gh_distance}}, for an observation that was \emph{not}
#' part of the calibration set. For a model with \code{ncomp} components fitted
#' on \code{n} calibration samples, the squared Mahalanobis distance of a new
#' observation has upper limit

#'  \mjeqn{\frac{p(n+1)(n-1)}{n(n-p)} F_{\alpha; p,\, n-p}}{p(n+1)(n-1)/(n(n-p)) F}.

#' Dividing by \code{p = ncomp} (the leverage normalization) gives the limit
#' returned here. Unlike an in-sample limit, this grows as \code{ncomp}
#' approaches \code{n}, reflecting the wider spread expected of new samples.
#' @param n the number of calibration observations.
#' @param ncomp the number of components.
#' @param conf the confidence level of the limit (default \code{0.95}).
#' @references
#' De Maesschalck, Roy, Delphine Jouan-Rimbaud, and Désiré L. Massart. "The mahalanobis distance."
#' Chemometrics and intelligent laboratory systems 50.1 (2000): 1-18.
#' @return A single numeric upper limit, or \code{NA_real_} if it cannot be
#' computed (non-finite inputs, or \code{ncomp} not in \code{1:(n - 1)}).
#' @keywords internal
.leverage_limit <- function(n, ncomp, conf = 0.95) {
  if (!is.finite(n) || !is.finite(ncomp) || ncomp < 1 || ncomp >= n) {
    return(NA_real_)
  }
  (n^2 - 1) / (n * (n - ncomp)) * stats::qf(conf, ncomp, n - ncomp)
}

#' @title Computes the NIRWise QVAL statistic
#' @description
#' QVAL indicates how different the predicted response variable (y) in
#' cross-validation deviates from the fitted version of y (i.e. the fitted y
#' values obtained when all calibration observations are used to fit the model).
#' @param y a matrix of one column with the response variable.
#' @param fitted_y a matrix with the estimated response variable for each
#' component.
#' @param predicted_y_in_cv the cross-validation estimates of the response
#' variable for every component.
#' @param scaled_scores a matrix of the scaled scores of the model.
#' @param ncomp a vector for each included component.
#' @seealso
#' \code{\link{calibrate}}
#' @return A list containing calibration statistics including residuals, predicted values, Mahalanobis distance, and Q-values.
#' @keywords internal
.calibration_statistics <- function(y, fitted_y, predicted_y_in_cv = NULL,
                                    scaled_scores, ncomp) {
  residual <- sweep(-fitted_y[, ncomp, drop = FALSE], MARGIN = 1, STATS = y, FUN = "+")

  if (!is.null(predicted_y_in_cv)) {
    cv_residual <- sweep(-predicted_y_in_cv[, ncomp, drop = FALSE], MARGIN = 1, STATS = y, FUN = "+")
    sd_res <- apply(
      residual,
      MARGIN = 2,
      FUN = function(x) {
        sqrt(sd(x)^2 * (length(x) - 1) / length(x))
      }
    )
    q_value <- abs(sweep(residual - cv_residual, MARGIN = 2, STATS = sd_res, FUN = "/"))
  } else {
    cv_residual <- q_value <- NULL
  }

  gh <- .gh_distance(x = scaled_scores, reference = scaled_scores, ncomp)

  # rownames(gh) <- rownames(scaled_scores)
  calibration_results <- list(
    Sample_index = seq_along(y),
    Target = y,
    fitted_y = fitted_y[, ncomp, drop = FALSE],
    residual = residual,
    predicted_y_in_cv = predicted_y_in_cv[, ncomp, drop = FALSE],
    cv_residual = cv_residual,
    Mahalanobis = gh,
    Q_value = q_value
  )
  calibration_results
}

#' @title Test if a string can be coerced to a numeric
#' @description
#' based on the code found at # https://stackoverflow.com/a/21154566/2292993
#' @return A logical vector indicating whether each element can be coerced to numeric.
#' @keywords internal
is_numeric_like <- function(x, na_strings = c("", ".", "NA", "na", "N/A", "n/a", "NaN", "nan")) {
  x <- trimws(x, "both")
  x[x %in% na_strings] <- NA_character_
  # https://stackoverflow.com/a/21154566/2292993
  result <- grepl("^[\\-\\+]?[0-9]+[\\.,]?[0-9]*$|^[\\-\\+]?[0-9]+[L]?$|^[\\-\\+]?[0-9]+[\\.,]?[0-9]*[eE][0-9]+$", x, perl = TRUE)
  result
}
