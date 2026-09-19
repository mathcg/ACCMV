\name{accmv_data}
\alias{accmv_data}
\alias{single_data_preparation}
\alias{multiple_data_preparation}
\title{Prepare ACCMV Data}
\usage{
accmv_data(x, y)
single_data_preparation(x, y)
multiple_data_preparation(x, y)
}
\arguments{
\item{x}{Numeric secondary-variable matrix with missing values represented by \code{NA}.}
\item{y}{Numeric primary-variable vector or matrix.}
}
\value{An object containing the data and stable missing-pattern encodings.}
\description{Validates data and identifies observed-value patterns.}
