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

###################
# .leverage_limit #
###################

test_that(".leverage_limit reproduces the new-observation leverage formula", {
  n <- 50
  for (p in c(1, 5, 20)) {
    expected <- (n^2 - 1) / (n * (n - p)) * qf(0.95, p, n - p)
    expect_equal(.leverage_limit(n, p), expected)
  }
})

test_that(".leverage_limit honours the requested confidence level", {
  expect_gt(.leverage_limit(50, 5, conf = 0.99), .leverage_limit(50, 5, conf = 0.95))
  expect_equal(.leverage_limit(50, 5, conf = 0.99), (50^2 - 1) / (50 * 45) * qf(0.99, 5, 45))
})

test_that(".leverage_limit grows as the components approach the sample size", {
  n <- 50
  # More components -> fewer residual df -> a wider (larger) new-observation limit.
  expect_gt(.leverage_limit(n, 45), .leverage_limit(n, 5))
})

test_that(".leverage_limit is a positive scalar", {
  out <- .leverage_limit(40, 3)
  expect_length(out, 1L)
  expect_true(is.finite(out))
  expect_gt(out, 0)
})

test_that(".leverage_limit returns NA when it cannot be computed", {
  expect_true(is.na(.leverage_limit(NA, 3)))
  expect_true(is.na(.leverage_limit(50, NA)))
  expect_true(is.na(.leverage_limit(50, 0)))   # ncomp < 1
  expect_true(is.na(.leverage_limit(50, 50)))  # ncomp == n
  expect_true(is.na(.leverage_limit(50, 60)))  # ncomp > n
})
