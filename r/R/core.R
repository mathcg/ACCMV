# Core estimators for the ACCMV package.

.as_matrix <- function(value, name) {
  value <- as.matrix(value)
  storage.mode(value) <- "double"
  if (length(dim(value)) != 2L) stop(name, " must be a numeric vector or matrix")
  if (any(is.infinite(value))) stop(name, " may contain NA but not infinity")
  value
}

.stable_patterns <- function(mask) {
  keys <- apply(mask, 1L, paste0, collapse = "")
  unique_keys <- unique(keys)
  first <- match(unique_keys, keys)
  list(patterns = mask[first, , drop = FALSE], index = match(keys, unique_keys))
}

#' Prepare ACCMV data
#'
#' @param x Numeric secondary-variable matrix; missing values are `NA`.
#' @param y Numeric primary-variable vector or matrix; missing values are `NA`.
#' @return An `accmv_data` object.
#' @export
accmv_data <- function(x, y) {
  x <- .as_matrix(x, "x")
  y <- .as_matrix(y, "y")
  if (nrow(x) != nrow(y)) stop("x and y must have the same number of rows")
  if (nrow(x) < 2L) stop("at least two observations are required")
  xp <- .stable_patterns(!is.na(x))
  yp <- .stable_patterns(!is.na(y))
  structure(list(x = x, y = y, missing_pattern_x = xp$patterns,
                 missing_pattern_y = yp$patterns, R = xp$index, A = yp$index,
                 num_pattern_x = nrow(xp$patterns),
                 num_pattern_y = nrow(yp$patterns)), class = "accmv_data")
}

#' @export
single_data_preparation <- function(x, y) accmv_data(x, y)

#' @export
multiple_data_preparation <- function(x, y) accmv_data(x, y)

.available <- function(patterns, observed) {
  if (!length(observed)) return(seq_len(nrow(patterns)))
  which(rowSums(patterns[, observed, drop = FALSE]) == length(observed))
}

.design <- function(x, rows, cols, y = NULL, ycols = integer()) {
  pieces <- list(`(Intercept)` = rep(1, sum(rows)))
  if (length(cols)) pieces$x <- x[rows, cols, drop = FALSE]
  if (!is.null(y) && length(ycols)) pieces$y <- y[rows, ycols, drop = FALSE]
  do.call(cbind, pieces)
}

.quadratic_design <- function(values) {
  values <- as.matrix(values); output <- cbind(`(Intercept)` = 1, values)
  if (ncol(values)) {
    for (left in seq_len(ncol(values))) {
      for (right in left:ncol(values)) {
        output <- cbind(output, values[, left] * values[, right])
      }
    }
  }
  output
}

.logit_fit <- function(design, response) {
  if (length(unique(response)) != 2L) {
    stop("each fitted odds model needs observations in both comparison groups")
  }
  fit <- suppressWarnings(glm.fit(design, response, family = binomial()))
  if (anyNA(fit$coefficients)) {
    fit$coefficients[is.na(fit$coefficients)] <- 0
  }
  fit$coefficients
}

.outcome_fit <- function(design, response, binary = FALSE) {
  if (binary) .logit_fit(design, response) else qr.solve(design, response, tol = 1e-10)
}

.predict_fit <- function(design, coefficients, binary = FALSE) {
  eta <- drop(design %*% coefficients)
  if (binary) plogis(eta) else eta
}

.single_y <- function(data) {
  if (ncol(data$y) != 1L) stop("single-primary estimators require one y column")
  drop(data$y)
}

#' Single-primary regression-adjustment estimator
#' @param data An object from [accmv_data()].
#' @param fun Transformation defining the estimand.
#' @param binary Whether the transformed response is binary.
#' @export
single_regression_adjustment <- function(data, fun = identity, binary = FALSE) {
  y <- .single_y(data); observed <- !is.na(y); outcome <- fun(y); n <- length(y)
  total <- sum(outcome[observed])
  for (rid in seq_len(data$num_pattern_x)) {
    recipient <- !observed & data$R == rid
    if (!any(recipient)) next
    cols <- which(data$missing_pattern_x[rid, ])
    donor <- observed & data$R %in% .available(data$missing_pattern_x, cols)
    coef <- .outcome_fit(.design(data$x, donor, cols), outcome[donor], binary)
    total <- total + sum(.predict_fit(.design(data$x, recipient, cols), coef, binary))
  }
  total / n
}

#' Single-primary inverse-probability-weighted estimator
#' @inheritParams single_regression_adjustment
#' @param delta Optional exponential-tilt sensitivity parameter.
#' @export
single_ipw <- function(data, fun = identity, delta = NULL) {
  y <- .single_y(data); observed <- !is.na(y); outcome <- fun(y); n <- length(y)
  contributions <- weights <- numeric(n)
  contributions[observed] <- outcome[observed]; weights[observed] <- 1
  for (rid in seq_len(data$num_pattern_x)) {
    recipient <- !observed & data$R == rid
    if (!any(recipient)) next
    cols <- which(data$missing_pattern_x[rid, ])
    donor <- observed & data$R %in% .available(data$missing_pattern_x, cols)
    selected <- donor | recipient
    coef <- .logit_fit(.design(data$x, selected, cols), as.numeric(recipient[selected]))
    odds <- exp(pmin(pmax(.predict_fit(.design(data$x, donor, cols), coef), -30), 30))
    if (!is.null(delta)) odds <- odds * exp(pmin(pmax(delta * y[donor], -30), 30))
    contributions[donor] <- contributions[donor] + odds * outcome[donor]
    weights[donor] <- weights[donor] + odds
  }
  sum(contributions) / if (is.null(delta)) n else sum(weights)
}

#' @export
single_ipw_sensitivity <- function(data, delta, fun = identity) {
  single_ipw(data, fun = fun, delta = delta)
}

#' Single-primary multiply-robust estimator
#' @inheritParams single_regression_adjustment
#' @export
single_multiply_robust <- function(data, fun = identity, binary = FALSE) {
  y <- .single_y(data); observed <- !is.na(y); outcome <- fun(y); n <- length(y)
  total <- sum(outcome[observed])
  for (rid in seq_len(data$num_pattern_x)) {
    recipient <- !observed & data$R == rid
    if (!any(recipient)) next
    cols <- which(data$missing_pattern_x[rid, ])
    donor <- observed & data$R %in% .available(data$missing_pattern_x, cols)
    selected <- donor | recipient
    odds_coef <- .logit_fit(.design(data$x, selected, cols), as.numeric(recipient[selected]))
    odds <- exp(pmin(pmax(.predict_fit(.design(data$x, donor, cols), odds_coef), -30), 30))
    out_coef <- .outcome_fit(.design(data$x, donor, cols), outcome[donor], binary)
    m0 <- .predict_fit(.design(data$x, donor, cols), out_coef, binary)
    m1 <- .predict_fit(.design(data$x, recipient, cols), out_coef, binary)
    total <- total + sum((outcome[donor] - m0) * odds) + sum(m1)
  }
  total / n
}

.complete_id <- function(data) {
  found <- which(rowSums(data$missing_pattern_y) == ncol(data$y))
  if (length(found) != 1L) stop("data must contain complete primary-variable cases")
  found
}

.target_values <- function(y, target, threshold = NULL) {
  if (ncol(y) != 2L) stop("multiple-primary targets require exactly two y columns")
  switch(target,
    average = 0.5 * (y[, 1] + y[, 2]),
    product = y[, 1] * y[, 2],
    indicator = {
      if (is.null(threshold)) stop("threshold is required for the indicator target")
      as.numeric(y[, 1] <= threshold & y[, 2] <= threshold)
    }, stop("target must be 'average', 'product', or 'indicator'"))
}

.components <- function(data) {
  complete_id <- .complete_id(data); complete <- data$A == complete_id; out <- list(); k <- 0L
  for (aid in seq_len(data$num_pattern_y)) {
    if (aid == complete_id) next
    ycols <- which(data$missing_pattern_y[aid, ])
    for (rid in seq_len(data$num_pattern_x)) {
      recipient <- data$A == aid & data$R == rid
      if (!any(recipient)) next
      xcols <- which(data$missing_pattern_x[rid, ])
      donor <- complete & data$R %in% .available(data$missing_pattern_x, xcols)
      if (!any(donor)) stop("a missingness pattern has no available complete-case donors")
      k <- k + 1L
      out[[k]] <- list(donor = donor, recipient = recipient, xcols = xcols, ycols = ycols)
    }
  }
  out
}

.conditional_target <- function(data, component, target, threshold) {
  donor <- component$donor; recipient <- component$recipient
  xcols <- component$xcols; ycols <- component$ycols; missing <- setdiff(1:2, ycols)
  d0 <- .design(data$x, donor, xcols, data$y, ycols)
  d1 <- .design(data$x, recipient, xcols, data$y, ycols)
  if (!length(ycols)) {
    response <- .target_values(data$y[donor, , drop = FALSE], target, threshold)
    binary <- target == "indicator"; coef <- .outcome_fit(d0, response, binary)
    if (target == "product") {
      d0 <- .quadratic_design(data$x[donor, xcols, drop = FALSE])
      d1 <- .quadratic_design(data$x[recipient, xcols, drop = FALSE])
      coef <- .outcome_fit(d0, response, FALSE)
    }
    return(list(m0 = .predict_fit(d0, coef, binary), m1 = .predict_fit(d1, coef, binary)))
  }
  if (target == "indicator") {
    response <- as.numeric(data$y[donor, missing] <= threshold)
    coef <- .outcome_fit(d0, response, TRUE)
    return(list(m0 = as.numeric(data$y[donor, ycols] <= threshold) * .predict_fit(d0, coef, TRUE),
                m1 = as.numeric(data$y[recipient, ycols] <= threshold) * .predict_fit(d1, coef, TRUE)))
  }
  coef <- .outcome_fit(d0, data$y[donor, missing])
  pred0 <- .predict_fit(d0, coef); pred1 <- .predict_fit(d1, coef)
  obs0 <- data$y[donor, ycols]; obs1 <- data$y[recipient, ycols]
  if (target == "average") list(m0 = .5 * (obs0 + pred0), m1 = .5 * (obs1 + pred1))
  else list(m0 = obs0 * pred0, m1 = obs1 * pred1)
}

.multiple_ipw_target <- function(data, target, threshold = NULL, delta = NULL, fun = NULL) {
  complete <- data$A == .complete_id(data); n <- nrow(data$y)
  outcome <- if (is.null(fun)) .target_values(data$y, target, threshold) else fun(data$y[, 1], data$y[, 2])
  contributions <- weights <- numeric(n); contributions[complete] <- outcome[complete]; weights[complete] <- 1
  for (z in .components(data)) {
    selected <- z$donor | z$recipient
    coef <- .logit_fit(.design(data$x, selected, z$xcols, data$y, z$ycols), as.numeric(z$recipient[selected]))
    odds <- exp(pmin(pmax(.predict_fit(.design(data$x, z$donor, z$xcols, data$y, z$ycols), coef), -30), 30))
    if (!is.null(delta)) {
      missing <- setdiff(1:2, z$ycols)
      tilt <- rowMeans(data$y[z$donor, missing, drop = FALSE])
      odds <- odds * exp(pmin(pmax(delta * tilt, -30), 30))
    }
    contributions[z$donor] <- contributions[z$donor] + odds * outcome[z$donor]
    weights[z$donor] <- weights[z$donor] + odds
  }
  sum(contributions) / if (is.null(delta)) n else sum(weights)
}

#' Multiple-primary IPW estimator
#' @param data An object from [accmv_data()].
#' @param fun Optional two-argument estimand function.
#' @param target Built-in target: `"average"`, `"product"`, or `"indicator"`.
#' @param threshold Threshold for the indicator target.
#' @param delta Optional exponential-tilt sensitivity parameter.
#' @export
multiple_ipw <- function(data, fun = NULL, target = "average", threshold = NULL, delta = NULL) {
  .multiple_ipw_target(data, target, threshold, delta, fun)
}

#' @export
multiple_ipw_sensitivity <- function(data, delta, fun = NULL, target = "average", threshold = NULL) {
  multiple_ipw(data, fun = fun, target = target, threshold = threshold, delta = delta)
}

.multiple_ra <- function(data, target, threshold = NULL) {
  complete <- data$A == .complete_id(data)
  total <- sum(.target_values(data$y, target, threshold)[complete])
  for (z in .components(data)) total <- total + sum(.conditional_target(data, z, target, threshold)$m1)
  total / nrow(data$y)
}

.multiple_mr <- function(data, target, threshold = NULL) {
  complete <- data$A == .complete_id(data); outcome <- .target_values(data$y, target, threshold)
  total <- sum(outcome[complete])
  for (z in .components(data)) {
    selected <- z$donor | z$recipient
    coef <- .logit_fit(.design(data$x, selected, z$xcols, data$y, z$ycols), as.numeric(z$recipient[selected]))
    odds <- exp(pmin(pmax(.predict_fit(.design(data$x, z$donor, z$xcols, data$y, z$ycols), coef), -30), 30))
    m <- .conditional_target(data, z, target, threshold)
    total <- total + sum((outcome[z$donor] - m$m0) * odds) + sum(m$m1)
  }
  total / nrow(data$y)
}

#' @export
multiple_ra_average <- function(data) .multiple_ra(data, "average")
#' @export
multiple_ra_product <- function(data) .multiple_ra(data, "product")
#' @export
multiple_ra_indicator <- function(data, a) .multiple_ra(data, "indicator", a)
#' @export
multiple_mr_average <- function(data) .multiple_mr(data, "average")
#' @export
multiple_mr_product <- function(data) .multiple_mr(data, "product")
#' @export
multiple_mr_indicator <- function(data, a) .multiple_mr(data, "indicator", a)

#' ACCMV weights for marginal regression
#' @param data An object from [accmv_data()].
#' @return Numeric weights, zero for incomplete primary-variable cases.
#' @export
accmv_ipw_weights <- function(data) {
  complete <- data$A == .complete_id(data); weights <- numeric(nrow(data$y)); weights[complete] <- 1
  for (z in .components(data)) {
    selected <- z$donor | z$recipient
    coef <- .logit_fit(.design(data$x, selected, z$xcols, data$y, z$ycols), as.numeric(z$recipient[selected]))
    weights[z$donor] <- weights[z$donor] + exp(pmin(pmax(.predict_fit(.design(data$x, z$donor, z$xcols, data$y, z$ycols), coef), -30), 30))
  }
  weights
}

#' @export
ipw_regression <- function(data) accmv_ipw_weights(data)

.regression_columns <- function(data, response, predictors) {
  response <- as.integer(response); predictors <- as.integer(predictors)
  if (length(response) != 1L || response < 1L || response > ncol(data$y)) {
    stop("response must be a valid one-based y column index")
  }
  if (!length(predictors) || any(predictors < 1L | predictors > ncol(data$y))) {
    stop("predictors must contain valid one-based y column indices")
  }
  if (response %in% predictors) stop("response may not also be a predictor")
  if (anyDuplicated(predictors)) stop("predictor indices must be unique")
  list(response = response, predictors = predictors)
}

.weighted_regression <- function(data, response, predictors) {
  columns <- .regression_columns(data, response, predictors)
  weights <- accmv_ipw_weights(data)
  complete <- rowSums(is.na(data$y)) == 0L
  design <- cbind(`(Intercept)` = 1, data$y[complete, columns$predictors, drop = FALSE])
  colnames(design) <- c("(Intercept)", paste0("y", columns$predictors))
  root_weight <- sqrt(weights[complete])
  coefficients <- qr.solve(design * root_weight,
                           data$y[complete, columns$response] * root_weight,
                           tol = 1e-10)
  names(coefficients) <- colnames(design)
  list(coefficients = drop(coefficients), weights = weights)
}

#' Fit an ACCMV-weighted marginal linear regression
#' @param x Numeric secondary-variable matrix.
#' @param y Numeric primary-variable matrix.
#' @param response One-based response column in `y`.
#' @param predictors One-based predictor columns in `y`.
#' @param n_boot Number of bootstrap replicates; zero disables bootstrap.
#' @param level Confidence level.
#' @param seed Optional bootstrap seed.
#' @return An `accmv_regression_result` object.
#' @export
fit_accmv_regression <- function(x, y, response = 2, predictors = 1,
                                 n_boot = 0, level = .95, seed = NULL) {
  data <- accmv_data(x, y)
  point <- .weighted_regression(data, response, predictors)
  output <- list(coefficients = point$coefficients, response = as.integer(response),
                 predictors = as.integer(predictors), nobs = nrow(data$y),
                 weights = point$weights)
  if (n_boot > 0) {
    if (!is.null(seed)) set.seed(seed)
    values <- list(); attempts <- 0L
    while (length(values) < n_boot && attempts < max(10L * n_boot, 100L)) {
      attempts <- attempts + 1L
      sampled <- .resample_data(data, sample.int(nrow(data$y), replace = TRUE))
      value <- try(.weighted_regression(sampled, response, predictors)$coefficients,
                   silent = TRUE)
      if (!inherits(value, "try-error") && all(is.finite(value))) {
        values[[length(values) + 1L]] <- value
      }
    }
    if (length(values) != n_boot) stop("too many bootstrap samples lacked estimable comparisons")
    bootstrap <- do.call(rbind, values); alpha <- 1 - level
    output$bootstrap <- bootstrap
    output$std.error <- apply(bootstrap, 2L, sd)
    output$conf.int <- t(apply(bootstrap, 2L, quantile,
                               probs = c(alpha / 2, 1 - alpha / 2), names = FALSE))
  } else if (n_boot < 0) stop("n_boot must be nonnegative")
  structure(output, class = "accmv_regression_result")
}

#' Bootstrap an ACCMV-weighted regression
#' @param fm Regression formula using columns of the primary-variable matrix.
#' @param method A model-fitting function such as [stats::lm()].
#' @inheritParams bootstrap_accmv_single
#' @export
bootstrap_regression <- function(data, method, fm, n_B = 1000, ...) {
  n <- nrow(data$y); output <- NULL
  for (i in seq_len(n_B)) {
    index <- sample.int(n, replace = TRUE); sampled <- .resample_data(data, index)
    weight <- accmv_ipw_weights(sampled); frame <- as.data.frame(sampled$y)
    fit <- method(formula = fm, data = frame, weights = weight, ...)
    output <- rbind(output, coef(fit))
  }
  output
}

.resample_data <- function(data, index) accmv_data(data$x[index, , drop = FALSE], data$y[index, , drop = FALSE])

#' Bootstrap a single-primary estimator
#' @param data An object from [accmv_data()].
#' @param method One of `"ipw"`, `"ra"`, or `"mr"`, or an estimator function.
#' @param n_B Number of accepted bootstrap replicates.
#' @param seed Optional random seed.
#' @param ... Arguments passed to the estimator.
#' @export
bootstrap_accmv_single <- function(data, method = "mr", n_B = 999, seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  fn <- if (is.function(method)) method else switch(method, ipw = single_ipw, ra = single_regression_adjustment, mr = single_multiply_robust, stop("unknown method"))
  values <- numeric(); attempts <- 0L
  while (length(values) < n_B && attempts < max(10L * n_B, 100L)) {
    attempts <- attempts + 1L; index <- sample.int(nrow(data$y), replace = TRUE)
    value <- try(fn(.resample_data(data, index), ...), silent = TRUE)
    if (!inherits(value, "try-error") && is.finite(value)) values <- c(values, value)
  }
  if (length(values) != n_B) stop("too many bootstrap samples lacked estimable comparisons")
  values
}

#' @export
single_bootstrap <- function(data, method, n_B = 1000, ...) bootstrap_accmv_single(data, method, n_B, ...)

#' @export
multiple_bootstrap <- function(data, method, n_B = 1000, ...) {
  values <- numeric(); attempts <- 0L
  while (length(values) < n_B && attempts < max(10L * n_B, 100L)) {
    attempts <- attempts + 1L; index <- sample.int(nrow(data$y), replace = TRUE)
    value <- try(method(.resample_data(data, index), ...), silent = TRUE)
    if (!inherits(value, "try-error") && is.finite(value)) values <- c(values, value)
  }
  if (length(values) != n_B) stop("too many bootstrap samples lacked estimable comparisons")
  values
}

#' Bootstrap a multiple-primary estimator
#' @inheritParams bootstrap_accmv_single
#' @param target Built-in estimand target.
#' @param threshold Threshold for the indicator target.
#' @export
bootstrap_accmv_multiple <- function(data, method = "mr", target = "average", threshold = NULL,
                                     n_B = 999, seed = NULL, delta = NULL) {
  if (!is.null(seed)) set.seed(seed)
  fn <- switch(method, ipw = function(d) multiple_ipw(d, target = target, threshold = threshold, delta = delta),
               ra = function(d) .multiple_ra(d, target, threshold),
               mr = function(d) .multiple_mr(d, target, threshold), stop("unknown method"))
  values <- numeric(); attempts <- 0L
  while (length(values) < n_B && attempts < max(10L * n_B, 100L)) {
    attempts <- attempts + 1L; index <- sample.int(nrow(data$y), replace = TRUE)
    value <- try(fn(.resample_data(data, index)), silent = TRUE)
    if (!inherits(value, "try-error") && is.finite(value)) values <- c(values, value)
  }
  if (length(values) != n_B) stop("too many bootstrap samples lacked estimable comparisons")
  values
}

.result <- function(estimate, method, target, n, boot = NULL, level = .95) {
  out <- list(estimate = estimate, method = method, target = target, nobs = n)
  if (!is.null(boot)) {
    alpha <- 1 - level; out$std.error <- sd(boot)
    out$conf.int <- unname(quantile(boot, c(alpha / 2, 1 - alpha / 2))); out$bootstrap <- boot
  }
  structure(out, class = "accmv_result")
}

#' Fit a single-primary ACCMV estimator
#' @param x Numeric secondary-variable matrix.
#' @param y Numeric primary-variable vector.
#' @param method `"ipw"`, `"ra"`, or `"mr"`.
#' @param fun Transformation defining the expectation.
#' @param binary Whether `fun(y)` is binary.
#' @param delta Optional IPW sensitivity parameter.
#' @param n_boot Number of bootstrap replicates; zero disables bootstrap.
#' @param level Confidence level.
#' @param seed Optional bootstrap seed.
#' @export
estimate_accmv_single <- function(x, y, method = "mr", fun = identity, binary = FALSE,
                                  delta = NULL, n_boot = 0, level = .95, seed = NULL) {
  data <- accmv_data(x, y)
  estimate <- switch(method,
    ipw = single_ipw(data, fun, delta),
    ra = {if (!is.null(delta)) stop("delta is only available for IPW"); single_regression_adjustment(data, fun, binary)},
    mr = {if (!is.null(delta)) stop("delta is only available for IPW"); single_multiply_robust(data, fun, binary)},
    stop("method must be 'ipw', 'ra', or 'mr'"))
  boot <- NULL
  if (n_boot > 0) {
    boot_method <- switch(method,
      ipw = function(d, ...) single_ipw(d, fun = fun, delta = delta),
      ra = function(d, ...) single_regression_adjustment(d, fun = fun, binary = binary),
      mr = function(d, ...) single_multiply_robust(d, fun = fun, binary = binary))
    boot <- bootstrap_accmv_single(data, boot_method, n_boot, seed)
  }
  .result(estimate, method, "expectation", nrow(data$y), boot, level)
}

#' Fit a two-primary ACCMV estimator
#' @inheritParams estimate_accmv_single
#' @param target `"average"`, `"product"`, or `"indicator"`.
#' @param threshold Threshold for the indicator target.
#' @export
estimate_accmv_multiple <- function(x, y, method = "mr", target = "average", threshold = NULL,
                                    delta = NULL, n_boot = 0, level = .95, seed = NULL) {
  data <- accmv_data(x, y)
  estimate <- switch(method,
    ipw = multiple_ipw(data, target = target, threshold = threshold, delta = delta),
    ra = {if (!is.null(delta)) stop("delta is only available for IPW"); .multiple_ra(data, target, threshold)},
    mr = {if (!is.null(delta)) stop("delta is only available for IPW"); .multiple_mr(data, target, threshold)},
    stop("method must be 'ipw', 'ra', or 'mr'"))
  boot <- if (n_boot > 0) bootstrap_accmv_multiple(data, method, target, threshold, n_boot, seed, delta) else NULL
  .result(estimate, method, target, nrow(data$y), boot, level)
}

#' @export
print.accmv_result <- function(x, ...) {
  cat("ACCMV", toupper(x$method), "estimate for", x$target, "\n")
  cat("Estimate:", format(x$estimate), "\n")
  if (!is.null(x$std.error)) cat("Bootstrap SE:", format(x$std.error), "\n")
  invisible(x)
}

#' @export
print.accmv_regression_result <- function(x, ...) {
  cat("ACCMV IPW weighted linear regression\n")
  print(x$coefficients)
  invisible(x)
}
