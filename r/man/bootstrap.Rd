\name{bootstrap_accmv_single}
\alias{bootstrap_accmv_single}
\alias{bootstrap_accmv_multiple}
\alias{single_bootstrap}
\alias{multiple_bootstrap}
\alias{bootstrap_regression}
\title{Bootstrap ACCMV Estimators}
\usage{
bootstrap_accmv_single(data, method = "mr", n_B = 999, seed = NULL, ...)
bootstrap_accmv_multiple(data, method = "mr", target = "average",
  threshold = NULL, n_B = 999, seed = NULL, delta = NULL)
single_bootstrap(data, method, n_B = 1000, ...)
multiple_bootstrap(data, method, n_B = 1000, ...)
bootstrap_regression(data, method, fm, n_B = 1000, ...)
}
\arguments{
\item{data}{An object from \code{accmv_data()}.}
\item{method}{Estimator name or a function; for regression bootstrap, a model-fitting function.}
\item{n_B}{Number of accepted bootstrap replicates.}
\item{seed}{Optional random seed.}
\item{target}{Built-in multiple-primary target.}
\item{threshold}{Threshold for the indicator target.}
\item{delta}{Optional exponential-tilt sensitivity parameter.}
\item{...}{Additional estimator arguments.}
\item{fm}{Regression formula using primary-variable columns.}
}
\value{A numeric vector of bootstrap estimates.}
\description{Performs a nonparametric row bootstrap, retrying draws lacking estimable pattern comparisons.}
