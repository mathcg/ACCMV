test_that("single-primary estimators recover the paper target", {
  dat <- simulate_accmv_single(5000, seed = 47)
  for (method in c("ra", "mr")) {
    fit <- estimate_accmv_single(dat$x, dat$y, method = method)
    expect_equal(fit$estimate, 89 / 96, tolerance = 0.15)
  }
  expect_true(is.finite(estimate_accmv_single(dat$x, dat$y, method = "ipw")$estimate))
})

test_that("multiple-primary product estimators recover the paper target", {
  dat <- simulate_accmv_multiple(6000, seed = 91)
  for (method in c("ipw", "ra", "mr")) {
    fit <- estimate_accmv_multiple(dat$x, dat$y, method = method, target = "product")
    expect_equal(fit$estimate, 175 / 128, tolerance = 0.2)
  }
})

test_that("bootstrap is reproducible", {
  dat <- simulate_accmv_single(800, seed = 9)
  one <- estimate_accmv_single(dat$x, dat$y, method = "ra", n_boot = 9, seed = 3)
  two <- estimate_accmv_single(dat$x, dat$y, method = "ra", n_boot = 9, seed = 3)
  expect_equal(one$bootstrap, two$bootstrap)
  expect_gt(one$std.error, 0)
})

test_that("R agrees with Python on the shared fixture", {
  fixture <- read.csv(system.file("extdata", "parity_single.csv", package = "accmv"), na.strings = "NA")
  expected <- c(ipw = 0.8258232120833385, ra = 0.9233429113247089, mr = 0.7096728452585138)
  for (method in names(expected)) {
    fit <- estimate_accmv_single(fixture[, c("x1", "x2")], fixture$y, method = method)
    expect_equal(fit$estimate, unname(expected[method]), tolerance = 1e-5)
  }
})
