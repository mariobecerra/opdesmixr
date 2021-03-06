% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mnl_model_functions.R
\name{mnl_get_opt_crit_value}
\alias{mnl_get_opt_crit_value}
\title{Computes optimality criterion value using the MNL model}
\usage{
mnl_get_opt_crit_value(X, beta, order, opt_crit = "D", transform_beta = T)
}
\arguments{
\item{X}{3 dimensional array of size \code{(q, J, S)} where: \itemize{
\item q is the number of ingredient proportions,
\item J is the number of alternatives within a choice set,
\item S is the number of choice sets.
}}

\item{beta}{numeric vector containing the parameters or numeric matrix containing draws of the prior distribution of the parameters.}

\item{order}{integer corresponding to a Scheffé model order (1, 2, 3).}

\item{opt_crit}{optimality criterion: 0 or "D" is D-optimality; while 1 or "I" is I-optimality.}

\item{transform_beta}{boolean parameter. Should the beta vector/matrix be transformed by subtracting the \code{q}-th element?}
}
\value{
Returns the value of the optimality criterion for this particular design and this \code{beta} vector
}
\description{
Computes optimality criterion value for a design array \code{X} and \code{beta} parameter.
}
