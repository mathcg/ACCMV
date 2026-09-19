\name{simulate_accmv_single}
\alias{simulate_accmv_single}
\alias{simulate_accmv_multiple}
\title{Paper Simulation Designs}
\usage{
simulate_accmv_single(n = 2000, seed = NULL)
simulate_accmv_multiple(n = 2000, seed = NULL)
}
\arguments{
\item{n}{Sample size.}
\item{seed}{Optional random seed.}
}
\value{A list containing \code{x} and \code{y}.}
\description{Generates data from Sections 7.1 and 7.2 of the ACCMV paper.}
