% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mnl_model_functions.R
\name{transform_varcov_matrix}
\alias{transform_varcov_matrix}
\title{Transform variance covariance matrix.}
\usage{
transform_varcov_matrix(Sigma, q)
}
\arguments{
\item{Sigma}{Original variance covariance matrix of size \code{(m, m)}.}

\item{q}{Index of the vector that is subtracted to the first q-1 entries.}
}
\value{
Matrix of size \code{(m-1, m-1)}.
}
\description{
This function transforms the variance covariance matrix of a multivariate normal vector to the identifiable space for a Scheffé model.
}
