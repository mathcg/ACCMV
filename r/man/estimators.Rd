\name{estimate_accmv_single}
\alias{estimate_accmv_single}
\alias{estimate_accmv_multiple}
\alias{single_regression_adjustment}
\alias{single_ipw}
\alias{single_ipw_sensitivity}
\alias{single_multiply_robust}
\alias{multiple_ipw}
\alias{multiple_ipw_sensitivity}
\alias{multiple_ra_average}
\alias{multiple_ra_product}
\alias{multiple_ra_indicator}
\alias{multiple_mr_average}
\alias{multiple_mr_product}
\alias{multiple_mr_indicator}
\alias{accmv_ipw_weights}
\alias{ipw_regression}
\title{ACCMV Estimators}
\usage{
estimate_accmv_single(x, y, method = "mr", fun = identity,
  binary = FALSE, delta = NULL, n_boot = 0, level = 0.95, seed = NULL)
estimate_accmv_multiple(x, y, method = "mr", target = "average",
  threshold = NULL, delta = NULL, n_boot = 0, level = 0.95, seed = NULL)
single_regression_adjustment(data, fun = identity, binary = FALSE)
single_ipw(data, fun = identity, delta = NULL)
single_ipw_sensitivity(data, delta, fun = identity)
single_multiply_robust(data, fun = identity, binary = FALSE)
multiple_ipw(data, fun = NULL, target = "average", threshold = NULL,
  delta = NULL)
multiple_ipw_sensitivity(data, delta, fun = NULL, target = "average",
  threshold = NULL)
multiple_ra_average(data)
multiple_ra_product(data)
multiple_ra_indicator(data, a)
multiple_mr_average(data)
multiple_mr_product(data)
multiple_mr_indicator(data, a)
accmv_ipw_weights(data)
ipw_regression(data)
}
\arguments{
\item{x}{Numeric secondary-variable matrix.}
\item{y}{Numeric primary-variable vector or two-column matrix.}
\item{method}{One of \code{"ipw"}, \code{"ra"}, or \code{"mr"}.}
\item{fun}{Transformation defining the estimand.}
\item{binary}{Whether the transformed outcome is binary.}
\item{delta}{Optional exponential-tilt sensitivity parameter.}
\item{n_boot}{Number of bootstrap replicates; zero disables bootstrap.}
\item{level}{Confidence level.}
\item{seed}{Optional random seed.}
\item{target}{One of \code{"average"}, \code{"product"}, or \code{"indicator"}.}
\item{threshold}{Threshold for the indicator target.}
\item{data}{An object from \code{accmv_data()}.}
\item{a}{Threshold for the indicator target.}
}
\value{High-level functions return an \code{accmv_result}; low-level functions return an estimate or weights.}
\description{Fits the estimators described in Cheng et al. (2022).}
