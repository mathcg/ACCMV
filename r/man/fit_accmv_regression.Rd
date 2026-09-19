\name{fit_accmv_regression}
\alias{fit_accmv_regression}
\title{Fit an ACCMV-Weighted Marginal Linear Model}
\usage{
fit_accmv_regression(x, y, response = 2, predictors = 1,
  n_boot = 0, level = 0.95, seed = NULL)
}
\arguments{
\item{x}{Numeric secondary-variable matrix.}
\item{y}{Numeric primary-variable matrix.}
\item{response}{One-based response column in \code{y}.}
\item{predictors}{One-based predictor columns in \code{y}.}
\item{n_boot}{Number of bootstrap replicates; zero disables bootstrap.}
\item{level}{Confidence level.}
\item{seed}{Optional bootstrap seed.}
}
\value{An \code{accmv_regression_result} with coefficients, weights, and optional bootstrap inference.}
\description{Fits a weighted least-squares marginal model using the ACCMV IPW weights from Section 5.}
