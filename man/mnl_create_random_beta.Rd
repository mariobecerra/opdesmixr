% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mnl_model_functions.R
\name{mnl_create_random_beta}
\alias{mnl_create_random_beta}
\title{Creation of a random parameter vector using the MNL model}
\usage{
mnl_create_random_beta(q, order = 3, seed = NULL)
}
\arguments{
\item{q}{integer specifying the number of mixture ingredients}

\item{order}{integer corresponding to a Scheffé model order (1, 2, 3).}

\item{seed}{integer used for reproducibility}
}
\value{
Returns a list in which the first element of the list is
a numerical vector with the parameters and the second element is a matrix with
the indices as in \emph{Bayesian D-optimal choice designs for mixtures} by Ruseckaite, Goos & Fok (2017).
}
\description{
Creates a random parameter vector with the appropriate dimensions as in \emph{Bayesian D-optimal choice designs for mixtures} by Ruseckaite, Goos & Fok (2017).
}
\examples{
mnl_create_random_beta(q = 3, order = 3, seed = 3)

}
