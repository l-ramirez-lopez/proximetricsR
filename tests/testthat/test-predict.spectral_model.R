data("NIRcannabis", package = "proximetricsR")

# Setup for checking if model predictions are as expected.
dat <- NIRcannabis[41:80, ] # reduce number of samples
X <- dat$spc # reduced number of samples
rownames(X) <- NULL
Y <- matrix(dat$THC, dimnames = list(41:80, "THC"))
method <- fit_plsr(10, "standard")
control <- calibration_control(validation_type = "kfold", number = 3, folds = "random", tuning_parameter = "rsq", seed = 42)
pretreats <- preprocess_recipe(
  prep_resample(c(1100, 1600, 5)),
  prep_derivative(m = 1, w = 5, p = 11, algorithm = "nwp"),
  prep_snv(),
  prep_smooth(w = 3, algorithm = "moving-average"),
  device = "unspecified"
)
# Index 50 cannot be skipped, will give a warning.
skiped_ind <- c(6, 14, 27, 31, 33, 50)

# model with formula; to check the predictions in detail
expect_warning(
  model_form <- calibrate(
    THC ~ spc,
    data = dat, preprocess = pretreats, method = method, control = control,
    skip_indices = skiped_ind, verbose = FALSE, metadata = add_model_metadata(unit = "%")
  ),
  "Unable to skip index 50, as it lies outside the considered indices."
)

# model with X, Y; to check the predictions in detail
expect_warning(
  model_xy <- calibrate(
    X, Y,
    data = dat, preprocess = pretreats, method = method, control = control,
    skip_indices = skiped_ind, verbose = FALSE, metadata = add_model_metadata()
  ),
  "Unable to skip index 50, as it lies outside the considered indices."
)

##################################
# CHECK PREDICTIONS ON SAME DATA #
##################################

test_that("Predictions on same dataset should be the same as the fitted values in the model", {
  # For model_form, data.frame as newdat
  expect_equal(
    model_form$final_model$model$fitted_y,
    predict(model_form, dat[-skiped_ind, ], ncomp = 1:10, verbose = FALSE)$predictions,
    tolerance = 1e-07
  )
  # For model_form, matrix as newdat
  expect_equal(
    unname(model_form$final_model$model$fitted_y),
    unname(predict(model_form, X[-skiped_ind, ], ncomp = 1:10, verbose = FALSE)$predictions),
    tolerance = 1e-07
  )
  # For model_xy, data.frame as newdat
  expect_equal(
    unname(model_xy$final_model$model$fitted_y),
    unname(predict(model_form, dat[-skiped_ind, ], ncomp = 1:10, verbose = FALSE)$predictions),
    tolerance = 1e-07
  )
  # For model_form, data.frame as newdat
  expect_equal(
    model_xy$final_model$model$fitted_y,
    predict(model_form, X[-skiped_ind, ], ncomp = 1:10, verbose = FALSE)$predictions,
    tolerance = 1e-07
  )
})

####################################################
# GENERATE PREDICTIONS USING DIFFERENT DATAFORMATS #
####################################################

# Predictions for model with formula, done with newdat being a data.frame, 6 components
predictions_df_form <- predict(
  model_form,
  dat[skiped_ind[1:(length(skiped_ind) - 1)], , drop = FALSE],
  ncomp = 1:6,
  verbose = FALSE
)

# Predictions for model with formula, done with newdat being a matrix, one component only
predictions_mat_form <- predict(
  model_form,
  X[skiped_ind[1:(length(skiped_ind) - 1)], , drop = FALSE],
  ncomp = 6,
  verbose = FALSE
)

# Predictions for model with X, Y, done with newdat being a data.frame, 1 component
predictions_df_mat <- predict(
  model_xy,
  dat[skiped_ind[1:(length(skiped_ind) - 1)], , drop = FALSE],
  ncomp = 4,
  verbose = FALSE
)

# Predictions for model with X, Y, done with newdat being a matrix, 7 components
predictions_mat_mat <- predict(
  model_xy,
  X[skiped_ind[1:(length(skiped_ind) - 1)], , drop = FALSE],
  ncomp = 4:10,
  verbose = FALSE
)

############################################
# CHECKING PREDICTIONS predictions_df_form #
############################################

test_that("Predictions with formula from a data.frame are correct", {
  withr::local_options(digits = 8)
  expect_snapshot(predictions_df_form$predictions)
})

test_that("Predictions with formula from a data.frame is of class 'spectral_prediction", {
  expect_true(inherits(predictions_df_form, "spectral_prediction"))
})

test_that("Predictions with formula from a data.frame is named correctly", {
  expect_named(predictions_df_form, c(
    "predictions", "scores", "mahalanobis", "q_residual",
    "q_limit", "leverage_limit", "control_limit_conf", "model_information"
  ))
})

test_that("Model information of predictions with formula from a data.frame is correct", {
  expect_identical(predictions_df_form$model_information$target_var, "THC")
  expect_identical(predictions_df_form$model_information$preprocess_recipe, model_form$preprocess)
  expect_identical(predictions_df_form$model_information$unit, "%")
  expect_identical(predictions_df_form$model_information$opt_comp, 10L)
  withr::local_options(digits = 8)
  expect_snapshot(predictions_df_form$model_information$model_grid)
})

test_that("Predictions for dataframes from a formula are correctly printed", {
  withr::local_options(digits = 8)
  expect_snapshot(print(predictions_df_form))
})

test_that("Mahalanobis distance of predictions with formula from a data.frame has the correct shape and values", {
  expect_true(is.matrix(predictions_df_form$mahalanobis))
  expect_identical(dim(predictions_df_form$mahalanobis), dim(predictions_df_form$predictions))
  expect_identical(dimnames(predictions_df_form$mahalanobis), dimnames(predictions_df_form$predictions))
  expect_true(all(predictions_df_form$mahalanobis >= 0))
})

test_that("Mahalanobis distance collapses to the squared scaled score for a single component", {
  predictions_1comp <- predict(
    model_form,
    dat[skiped_ind[1:(length(skiped_ind) - 1)], , drop = FALSE],
    ncomp = 1,
    verbose = FALSE
  )
  ref_scaled_scores <- model_form$final_model$model$scaled_scores[, 1]
  new_scaled_scores <- predictions_1comp$scores[, 1] / model_form$final_model$model$sd_scores[1]
  expected <- new_scaled_scores^2 / var(ref_scaled_scores)
  expect_equal(unname(predictions_1comp$mahalanobis[, 1]), unname(expected), tolerance = 1e-10)
})

test_that("Spectral (Q) residual of predictions with formula from a data.frame has the correct shape and values", {
  expect_true(is.matrix(predictions_df_form$q_residual))
  expect_identical(dim(predictions_df_form$q_residual), dim(predictions_df_form$predictions))
  expect_identical(dimnames(predictions_df_form$q_residual), dimnames(predictions_df_form$predictions))
  expect_true(all(predictions_df_form$q_residual >= 0))
})

test_that("Spectral (Q) residual matches a hand-reconstructed residual sum of squares", {
  newdata_raw <- X[skiped_ind[1:(length(skiped_ind) - 1)], , drop = FALSE]
  newdata_processed <- process(newdata_raw, model_xy$preprocess)
  new_data_centered <- scale(newdata_processed, center = model_xy$final_model$model$x_means, scale = FALSE)
  nc <- 4
  scores_new <- new_data_centered %*% t(model_xy$final_model$model$projection_m)
  x_hat <- scores_new[, 1:nc, drop = FALSE] %*% model_xy$final_model$model$x_loadings[1:nc, , drop = FALSE]
  expected_q <- rowSums((new_data_centered - x_hat)^2)
  expect_equal(unname(predictions_mat_mat$q_residual[, "ncomp_4"]), unname(expected_q), tolerance = 1e-8)
})

test_that("Predicting a single new sample does not error and returns proper matrices", {
  single_pred <- predict(
    model_form,
    dat[skiped_ind[1], , drop = FALSE],
    ncomp = 1:6,
    verbose = FALSE
  )
  expect_true(is.matrix(single_pred$mahalanobis))
  expect_identical(dim(single_pred$mahalanobis), c(1L, 6L))
  expect_true(is.matrix(single_pred$q_residual))
  expect_identical(dim(single_pred$q_residual), c(1L, 6L))
})

####################################################
# CONTROL LIMITS (Q / leverage) CARRIED BY PREDICT #
####################################################

test_that("Control limits are one per requested component and named accordingly", {
  expect_length(predictions_df_form$q_limit, ncol(predictions_df_form$predictions))
  expect_length(predictions_df_form$leverage_limit, ncol(predictions_df_form$predictions))
  expect_named(predictions_df_form$q_limit, colnames(predictions_df_form$predictions))
  expect_named(predictions_df_form$leverage_limit, colnames(predictions_df_form$predictions))
})

test_that("Control limits are finite, positive and use the 99% confidence", {
  expect_true(all(is.finite(predictions_df_form$q_limit)))
  expect_true(all(is.finite(predictions_df_form$leverage_limit)))
  expect_true(all(predictions_df_form$q_limit > 0))
  expect_true(all(predictions_df_form$leverage_limit > 0))
  expect_identical(predictions_df_form$control_limit_conf, 0.99)
})

test_that("Q limit is the calibration 99th percentile; leverage limit is the new-observation Mahalanobis limit", {
  cal_q <- model_form$final_model$model$x_residuals
  n_cal <- model_form$final_model$model$n_observations
  for (nc in 1:6) {
    expect_equal(
      unname(predictions_df_form$q_limit[nc]),
      unname(quantile(cal_q[, nc], 0.99))
    )
    expect_equal(
      unname(predictions_df_form$leverage_limit[nc]),
      .leverage_limit(n_cal, nc, 0.99)
    )
  }
})

test_that("Control limits only cover the requested components (aligned to ncomp)", {
  # ncomp = 4:10 -> the fourth column's limit equals the 4-component calibration limit.
  cal_q <- model_xy$final_model$model$x_residuals
  expect_equal(
    unname(predictions_mat_mat$q_limit["ncomp_4"]),
    unname(quantile(cal_q[, 4], 0.99))
  )
  expect_named(predictions_mat_mat$q_limit, colnames(predictions_mat_mat$predictions))
})

test_that("control_limit_conf is configurable and flows into the limits", {
  cal_q <- model_form$final_model$model$x_residuals
  n_cal <- model_form$final_model$model$n_observations
  pred90 <- predict(
    model_form, dat[skiped_ind[1:5], , drop = FALSE],
    ncomp = 1:6, verbose = FALSE, control_limit_conf = 0.90
  )
  expect_identical(pred90$control_limit_conf, 0.90)
  for (nc in 1:6) {
    expect_equal(unname(pred90$q_limit[nc]), unname(quantile(cal_q[, nc], 0.90)))
    expect_equal(unname(pred90$leverage_limit[nc]), .leverage_limit(n_cal, nc, 0.90))
  }
  # A lower confidence must give tighter (smaller) limits than the default 0.99.
  expect_true(all(pred90$leverage_limit < predictions_df_form$leverage_limit))
})

test_that("control_limit_conf must be a single number strictly between 0 and 1", {
  nd <- dat[skiped_ind[1:5], , drop = FALSE]
  expect_error(predict(model_form, nd, verbose = FALSE, control_limit_conf = 1), "between 0 and 1")
  expect_error(predict(model_form, nd, verbose = FALSE, control_limit_conf = 0), "between 0 and 1")
  expect_error(predict(model_form, nd, verbose = FALSE, control_limit_conf = c(0.9, 0.95)), "between 0 and 1")
  expect_error(predict(model_form, nd, verbose = FALSE, control_limit_conf = "0.95"), "between 0 and 1")
})

#########################################################################
# PREDICTIONS OF predictions_mat_form CORRESPOND TO predictions_df_form #
#########################################################################

test_that("Predictions are the same for matrix/data.frame as newdata", {
  expect_identical(
    predictions_mat_form$model_information$model_grid,
    predictions_df_form$model_information$model_grid[6, , drop = FALSE]
  )
  expect_equal(
    unname(predictions_mat_form$predictions),
    unname(predictions_df_form$predictions[, 6, drop = FALSE]),
    tolerance = 1e-14
  )
})

#######################################################################
# PREDICTIONS OF predictions_df_mat CORRESPOND TO predictions_df_form #
#######################################################################

test_that("Predictions are the same for model with matrices/model with formula", {
  expect_identical(
    predictions_df_mat$model_information$model_grid,
    predictions_df_form$model_information$model_grid[4, , drop = FALSE]
  )
  expect_equal(
    predictions_df_mat$predictions,
    predictions_df_form$predictions[, 4, drop = FALSE],
    tolerance = 1e-14
  )
})

############################################
# CHECKING PREDICTIONS predictions_mat_mat #
############################################

test_that("Predictions with matrices from a matrix are correct", {
  withr::local_options(digits = 8)
  expect_snapshot(predictions_mat_mat$predictions)
})

test_that("Predictions with matrices from a matrix is of class 'spectral_prediction", {
  expect_true(inherits(predictions_mat_mat, "spectral_prediction"))
})

test_that("Predictions with matrices from a matrix is named correctly", {
  expect_named(predictions_mat_mat, c(
    "predictions", "scores", "mahalanobis", "q_residual",
    "q_limit", "leverage_limit", "control_limit_conf", "model_information"
  ))
})

test_that("Model information of predictions with matrices from a matrix is correct", {
  expect_identical(predictions_mat_mat$model_information$target_var, "THC")
  expect_identical(predictions_mat_mat$model_information$preprocess_recipe, model_form$preprocess)
  expect_identical(predictions_mat_mat$model_information$unit, "")
  expect_identical(predictions_mat_mat$model_information$opt_comp, 10L)
  withr::local_options(digits = 8)
  expect_snapshot(predictions_mat_mat$model_information$model_grid)
})

test_that("Predictions of a dataframe from a matrix are correctly printed", {
  withr::local_options(digits = 8)
  expect_snapshot(print(predictions_df_mat))
})

test_that("Mahalanobis distance of predictions with matrices from a matrix has the correct shape and values", {
  expect_true(is.matrix(predictions_mat_mat$mahalanobis))
  expect_identical(dim(predictions_mat_mat$mahalanobis), dim(predictions_mat_mat$predictions))
  expect_identical(dimnames(predictions_mat_mat$mahalanobis), dimnames(predictions_mat_mat$predictions))
  expect_true(all(predictions_mat_mat$mahalanobis >= 0))
})

#################
# SANITY CHECKS #
#################

test_that("For predictions, newdata must be given", {
  expect_error(predict(model_form, verbose = FALSE))
})

test_that("For predictions, the model given must be of class 'spectral_model'", {
  expect_error(
    predict.spectral_model(list(), newdata = NIRcannabis),
    "'object' must be of class 'spectral_model."
  )
})

test_that("For predictions, the number of components must be in the computed components of the model", {
  expect_error(
    predict(model_form, newdata = NIRcannabis, ncomp = c(20, 21), verbose = FALSE),
    "The maximum of 'ncomp' is larger than the number of components in the model"
  )
  expect_error(
    predict(model_form, newdata = NIRcannabis, ncomp = 20, verbose = FALSE),
    "'ncomp' is larger than the number of components in the model"
  )
})

test_that("For predictions, newdata must be a data.frame or matrix", {
  expect_error(predict(model_form, newdata = c(NIRcannabis), verbose = FALSE))
})

test_that("For predictions with formula, the predictor variables must be contained in newdata", {
  expect_error(
    predict(model_form, newdata = unname(NIRcannabis), verbose = FALSE),
    "The following predictor variables are missing: spc"
  )
})

test_that("Verbose must be a logical", {
  expect_error(
    predict(model_form, newdata = NIRcannabis, verbose = "1"),
    "'verbose' must a logical."
  )
})

test_that("ncomp must be a vector of numerics", {
  expect_error(
    predict(model_form, newdata = NIRcannabis, ncomp = "1")
  )
})

test_that("Predictor variables (newdata) must have the same wavelengths as the model", {
  dt <- dat
  s_model <- calibrate(
    THC ~ spc,
    data = dt, preprocess = preprocess_recipe(), method = method, control = control, verbose = FALSE, metadata = add_model_metadata()
  )
  colnames(dt$spc) <- 1:ncol(dt$spc)
  expect_error(
    predict(s_model, newdata = dt$spc, verbose = FALSE),
    "Missing predictor variables"
  )
})

##################################################
# FROZEN CONSTANT-EDGE TRIMMING (end-to-end)     #
##################################################

test_that("constant-edge trimming is frozen at calibration and keeps newdata aligned", {
  skip_on_cran()

  # Force a constant left edge so training edge-trimming drops the first columns.
  Xtr <- X
  Xtr[, 1:3] <- Xtr[, 4]

  trim_recipe <- preprocess_recipe(
    prep_wav_trim(band = c(), trim_constant_edges = TRUE),
    prep_snv(),
    device = "proxiscout"
  )

  m <- suppressWarnings(calibrate(
    Xtr, Y,
    preprocess = trim_recipe, method = fit_plsr(5, "standard"),
    control = control, verbose = FALSE
  ))

  # (1) The recipe stored on the model is frozen to the exact retained wavelengths.
  frozen_wavs <- m$preprocess$steps[[1]]$resolved_wavs
  expect_false(is.null(frozen_wavs))
  expect_equal(as.numeric(m$predictor_variables), frozen_wavs)

  # (2) newdata with NON-constant edges: a re-derived scan would keep the first
  #     three columns and misalign. The frozen recipe drops them and predicts.
  Xnew <- NIRcannabis$spc[1:15, ]
  pred <- predict(m, Xnew, verbose = FALSE)
  expect_s3_class(pred, "spectral_prediction")
  expect_equal(nrow(pred$predictions), nrow(Xnew))

  # (3) Missing a required wavelength errors clearly instead of misaligning.
  expect_error(
    predict(m, Xnew[, -10], verbose = FALSE),
    "missing"
  )
})

test_that("predict errors clearly when newdata lacks model spectral variables", {
  skip_on_cran()

  m0 <- suppressWarnings(calibrate(
    X, Y,
    preprocess = preprocess_recipe(), method = fit_plsr(5, "standard"),
    control = control, verbose = FALSE
  ))

  # Named matrix missing one required wavelength.
  expect_error(
    predict(m0, X[, -5], verbose = FALSE),
    "is missing 1 spectral variable"
  )

  # Unnamed matrix: colnames() is NULL, so the guard must still fire.
  expect_error(
    predict(m0, unname(X), verbose = FALSE),
    "spectral variable"
  )
})
