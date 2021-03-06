% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mnl_model_functions.R
\name{mnl_mixture_coord_exch}
\alias{mnl_mixture_coord_exch}
\title{Coordinate exchange algorithm for a Multinomial Logit Scheffé model.}
\usage{
mnl_mixture_coord_exch(
  q = NULL,
  J = NULL,
  S = NULL,
  n_random_starts = 100,
  X = NULL,
  beta,
  transform_beta = T,
  order = 3,
  opt_method = "B",
  max_it = 10,
  tol = 1e-04,
  n_cox_points = NULL,
  plot_designs = F,
  verbose = 1,
  opt_crit = 0,
  seed = NULL,
  n_cores = 1,
  save_all_designs = F
)
}
\arguments{
\item{q}{number of mixture ingredient proportions.}

\item{J}{number of alternatives within a choice set.}

\item{S}{number of choice sets.}

\item{n_random_starts}{number or random starts. Defaults to 100.}

\item{X}{If an initial design is to be supplied, then it must be a 3-dimensional array of size \code{(q, J, S)}, with q, J, and S defined above.}

\item{beta}{Prior parameters. For a locally optimal design, it should be a numeric vector of length m = (q^3 + 5*q)/6. For a pseudo-Bayesian design, it must be a matrix with prior simulations of size (nxm) where m is previously defined and m is the number of prior draws, i.e., there is a prior draw per row.}

\item{transform_beta}{boolean parameter. Should the \code{beta} vector/matrix be transformed by subtracting the q-th element?}

\item{order}{integer corresponding to a Scheffé model order (1, 2, 3).}

\item{opt_method}{Optimization method in each step of the coordinate exchange algorithm.
It can be "B" (Brent's algorithm) or "D" (discretization of Cox direction)}

\item{max_it}{integer for maximum number of iterations that the coordinate exchange algorithm will do}

\item{tol}{A positive error tolerance in Brent's method.}

\item{n_cox_points}{number of points to use in the discretization of Cox direction. Ignored if opt_method is Brent.}

\item{plot_designs}{boolean. If TRUE, shows a plot of the initial and the final design. Only works if q is 3 or 4.}

\item{verbose}{level of verbosity. See below for details.}

\item{opt_crit}{optimality criterion: D-optimality ("D" or 0) or I-optimality ("I" or 1).}

\item{seed}{Seed for reproducibility.}

\item{n_cores}{Number of cores for parallel processing.}

\item{save_all_designs}{Whether the function should return a list with all the designs created at random or only the best.}
}
\value{
The function returns a list with 11 elements: \enumerate{
\item \code{X_orig}: The original design. A 3-dimensional array of size \code{(q, J, S)}.
\item \code{X}: The optimized design. A 3-dimensional array of size \code{(q, J, S)}.
\item \code{beta}: The original \code{beta} vector or matrix.
\item \code{opt_crit_value_orig}: efficiency of the original design.
\item \code{opt_crit_value}: efficiency of the optimized design.
\item \code{n_iter}: Number of iterations performed.
\item \code{efficiency_value_per_iteration}: Efficiency value in each iteration of the algorithm.
\item \code{opt_crit}: The optimality criterion used.
\item \code{q}: Number of mixture ingredients.
\item \code{seed}: seed used to generate the final design. If a design was used as input by the user, this will be NA.
}
}
\description{
Performs the coordinate exchange algorithm for a Multinomial Logit Scheffé model as described in
\emph{Bayesian D-optimal choice designs for mixtures} by Ruseckaite, Goos & Fok (2017).
}
\details{
Verbosity levels: each level prints the previous plus additional things:
\enumerate{
\item Print the efficiency value in each iteration and a final summary
\item Print the values of k, s, i, and efficiency value in each subiteration
\item Print the resulting X after each iteration, i.e., after each complete pass on the data
\item Print efficiency value for each point in the Cox direction discretization
\item Print the resulting X and information matrix after each subiteration
\item Print the resulting X or each point in the Cox direction discretization
}
}
