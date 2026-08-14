reference <- matrix(rnorm(50 * 5), 50, 5)

test_that(".gh_distance returns a proper matrix for multiple samples and multiple components", {
  x <- matrix(rnorm(5 * 5), 5, 5)
  out <- .gh_distance(x, reference, 1:3)
  expect_true(is.matrix(out))
  expect_identical(dim(out), c(5L, 3L))
})

test_that(".gh_distance returns a proper matrix for a single component count", {
  x <- matrix(rnorm(5 * 5), 5, 5)
  out <- .gh_distance(x, reference, 1)
  expect_true(is.matrix(out))
  expect_identical(dim(out), c(5L, 1L))
})

test_that(".gh_distance returns a proper matrix for a single sample", {
  x <- matrix(rnorm(1 * 5), 1, 5)
  out <- .gh_distance(x, reference, 1:3)
  expect_true(is.matrix(out))
  expect_identical(dim(out), c(1L, 3L))
})

test_that(".gh_distance returns a proper matrix for a single sample and a single component count", {
  x <- matrix(rnorm(1 * 5), 1, 5)
  out <- .gh_distance(x, reference, 1)
  expect_true(is.matrix(out))
  expect_identical(dim(out), c(1L, 1L))
})

test_that(".gh_distance values match a hand-computed Mahalanobis distance", {
  x <- matrix(rnorm(5 * 5), 5, 5)
  out <- .gh_distance(x, reference, 3)
  expected <- stats::mahalanobis(
    x[, 1:3, drop = FALSE],
    center = rep(0, 3),
    cov = cov(reference[, 1:3, drop = FALSE])
  ) / 3
  expect_equal(unname(out[, 1]), unname(expected))
})
