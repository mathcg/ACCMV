.accmv_sigma <- function(d) 0.5 * diag(d) + 0.5 * matrix(1, d, d)

#' Simulate the paper's single-primary experiment
#' @param n Sample size.
#' @param seed Optional random seed.
#' @return A list with `x` and `y`; the true mean is 89/96.
#' @export
simulate_accmv_single <- function(n = 2000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  patterns <- rbind(c(FALSE, FALSE), c(FALSE, TRUE), c(TRUE, FALSE), c(TRUE, TRUE))
  a <- sample(0:1, n, TRUE); r <- sample(1:4, n, TRUE)
  x <- matrix(NA_real_, n, 2); y <- rep(NA_real_, n)
  means <- list(`1` = 1, `2` = c(1, -1), `3` = c(0, -1, -1))
  for (i in seq_len(n)) {
    cols <- which(patterns[r[i], ])
    if (a[i] == 1) {
      d <- length(cols) + 1L
      draw <- drop(means[[as.character(d)]] + t(chol(.accmv_sigma(d))) %*% rnorm(d))
      y[i] <- draw[1]; if (length(cols)) x[i, cols] <- draw[-1]
    } else if (length(cols)) {
      d <- length(cols)
      draw <- drop(means[[as.character(d)]] + t(chol(.accmv_sigma(d))) %*% rnorm(d))
      x[i, cols] <- draw
    }
  }
  list(x = x, y = y)
}

#' Simulate the paper's multiple-primary experiment
#' @inheritParams simulate_accmv_single
#' @return A list with `x` and `y`; the true product moment is 175/128.
#' @export
simulate_accmv_multiple <- function(n = 2000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  patterns <- rbind(c(FALSE, FALSE), c(FALSE, TRUE), c(TRUE, FALSE), c(TRUE, TRUE))
  a <- sample(1:4, n, TRUE); r <- sample(1:4, n, TRUE)
  x <- y <- matrix(NA_real_, n, 2)
  for (i in seq_len(n)) {
    xcols <- which(patterns[r[i], ]); ycols <- which(patterns[a[i], ])
    if (a[i] == 4) {
      d <- 2L + length(xcols); draw <- drop(rep(1, d) + t(chol(.accmv_sigma(d))) %*% rnorm(d))
      y[i, ] <- draw[1:2]; if (length(xcols)) x[i, xcols] <- draw[-(1:2)]
    } else if (length(ycols) == 1L) {
      d <- 1L + length(xcols); mu <- rep(if (d == 1L) .5 else 1, d)
      draw <- drop(mu + t(chol(.accmv_sigma(d))) %*% rnorm(d))
      y[i, ycols] <- draw[1]; if (length(xcols)) x[i, xcols] <- draw[-1]
    } else if (length(xcols)) {
      d <- length(xcols); mu <- rep(if (d == 1L) .5 else 1, d)
      x[i, xcols] <- drop(mu + t(chol(.accmv_sigma(d))) %*% rnorm(d))
    }
  }
  list(x = x, y = y)
}

#' Simulate the paper's marginal-regression experiment
#' @inheritParams simulate_accmv_single
#' @return A list with `x` and `y`; the true coefficients are -1 and 0.5.
#' @export
simulate_accmv_regression <- function(n = 2000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (n < 1L) stop("n must be positive")
  full <- matrix(c(1, 0, -1), n, 3L, byrow = TRUE) +
    matrix(rnorm(3L * n), n, 3L) %*% chol(.accmv_sigma(3L))
  ratio <- exp(pmin(pmax(.5 * full[, 1], -30), 30))
  base <- 1 / (5 + 3 * ratio)
  category <- integer(n)
  for (i in seq_len(n)) {
    category[i] <- sample.int(8L, 1L,
      prob = c(rep(base[i], 4L), rep(base[i] * ratio[i], 3L), base[i]))
  }
  patterns <- rbind(c(FALSE, FALSE), c(FALSE, TRUE), c(TRUE, FALSE), c(TRUE, TRUE))
  r <- as.integer(category >= 5L)
  a <- ifelse(category <= 4L, category, ifelse(category <= 7L, category - 4L, 4L))
  x <- full[, 1, drop = FALSE]; y <- full[, 2:3, drop = FALSE]
  x[r == 0L, 1] <- NA_real_
  for (i in seq_len(n)) y[i, !patterns[a[i], ]] <- NA_real_
  list(x = x, y = y)
}
