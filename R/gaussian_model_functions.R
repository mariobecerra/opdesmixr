
#' TODO: write doc
#' @export
create_random_initial_design_gaussian = function(n_runs, q, seed = NULL){
  X = matrix(rep(NA_real_, n_runs*q), nrow = n_runs)

  if(!is.null(seed)) set.seed(seed)

  for(i in 1:nrow(X)){
    rands = runif(q)

    # rows in X sum to 1:
    X[i,] = rands/sum(rands)
  }

  return(X)
}




#' Coordinate exchange algorithm for a mixture model assuming Gaussian iid errors.
#'
#' \code{mixture_coord_ex_gaussian} Performs the coordinate exchange algorithm for a Scheffé mixture model. It can have many different random starts, or the coordinate exchange algorithm can be performed in a user-supplied matrix.
#' @param n_runs number of runs
#' @param q number of ingredient proportions
#' @param n_random_starts number or random starts. Defaults to 100.
#' @param X User supplied design matrix. Must be of size (n_runs, q) where:
#'     n_runs is the number of runs
#'     q is the number of ingredient proportions.
#' @param order Order of the Scheffé model (1, 2, or 3).
#' @param opt_method Optimization method in each step of the coordinate exchange algorithm.
#'      It can be "B" (Brent's algorithm) or "D" (discretization of Cox direction)
#' @param max_it integer for maximum number of iterations that the coordinate exchange algorithm will do
#' @param tol A positive error tolerance in Brent's method.
#' @param n_cox_points number of points to use in the discretization of Cox direction
#' @param plot_designs boolean. If TRUE, shows a plot of the initial and the final design. Only works if q is 3.
#' @param verbose level of verbosity.
#' @param opt_crit optimality criterion: D-optimality ("D" or 0) or I-optimality ("I" or 1)
#' @param seed Seed for reproducibility
#' @param n_cores Number of cores for parallel processing
#' @return list with 6 elements. See below for details.
#'
#'
#' Return list has 6 elements:
#' \itemize{
#'     \item X_orig: The original design. Matrix of size (n_runs, q).
#'     \item X: The optimized design. Matrix of size (n_runs, q).
#'     \item opt_crit_value_orig: efficiency of the original design.
#'     \item opt_crit_value: efficiency of the optimized design.
#'     \item n_iter: Number of iterations performed.
#'     \item opt_crit: The optimality criterion used.
#'  }
#'
#' @export
mixture_coord_ex_gaussian = function(
  n_runs = NULL,
  q = NULL,
  n_random_starts = 100,
  X = NULL,
  order = 1,
  opt_method = "B",
  max_it = 10,
  tol = 0.0001,
  n_cox_points = NULL,
  plot_designs = F,
  verbose = 1,
  opt_crit = 0,
  seed = NULL,
  n_cores = 1){


  #############################################
  ## Check that optimality criterion is okay
  #############################################
  if(!(opt_crit %in% c(0, 1) | opt_crit %in% c("D", "I"))){
    stop('Unknown optimality criterion. Must be either "D" or 0 for D-optimality, or "I" or 1 for I-optimality.' )
  }

  # Recode opt_crit
  if(opt_crit == "D") opt_crit = 0
  if(opt_crit == "I") opt_crit = 1



  if(!(opt_method %in% c("B", "D"))){
    stop('Unknown optimization method. Must be either "B" for Brent or "D" for discretization of Cox direction.' )
  }



  #############################################
  ## Check that optimization method is okay
  #############################################

  if(opt_method == "B" & !is.null(n_cox_points)){
    warning("n_cox_points provided but ignoring because optimization method is Brent.")
  }

  # Make n_cox_points an integer otherwise C++ will throw an error
  if(is.null(n_cox_points)) n_cox_points = 2

  if(opt_method == "B") opt_method = 0
  if(opt_method == "D") opt_method = 1

  # if(opt_method == "D" & !is.null(tol)){
  #   warning("tol provided but ignoring because optimization method is discretization of Cox direction.")
  # }


  #############################################
  ## Create random initial designs or check that
  ## the provide design is okay.
  #############################################
  if(is.null(X)){
    if(!is.null(seed)) set.seed(seed)

    seeds_designs = sample.int(1e9, n_random_starts)

    designs = lapply(seeds_designs, function(x){
      des = create_random_initial_design_gaussian(n_runs, q, seed = x)
      return(des)
    })
  } else{
    n_runs = nrow(X)
    q = ncol(X)

    # Check that rows in X sum to 1
    row_sums = apply(X, 1, sum)

    if(sum(abs(row_sums - rep(1, n_runs)) < 1e-10) != n_runs){
      stop("Rows in X must sum 1")
    }

    designs = list(X)
  }



  #############################################
  ## Moments matrices
  #############################################

  # If criterion is D-optimality, send a matrix with only one zero element as a moment matrix
  # Maybe a NULL value would be better. Gotta check.
  if(opt_crit == 0){
    # "D-optimality"
    W = matrix(0.0, nrow = 1)
  } else{
    # "I-optimality")
    W = create_moment_matrix_gaussian(q)
  }

  #############################################
  ## Check operating system for parallel processing
  #############################################

  if(.Platform$OS.type != "unix") n_cores = 1


  #############################################
  ## Apply the coordinate exchange algorithm
  #############################################

  results = parallel::mclapply(seq_along(designs), function(i){
    X = designs[[i]]

    if(verbose > 0) cat("\nDesign", i, "\n")

    out = try(
      mixtureCoordinateExchangeGaussian(
        X_orig = X,
        order = order,
        max_it = max_it,
        verbose = verbose,
        opt_crit = opt_crit,
        W = W,
        opt_method = opt_method,
        lower = 0,
        upper = 1,
        tol = tol,
        n_cox_points = n_cox_points),
      silent = T)


    # If there was an error with this design, return whatever
    if(class(out) == "try-error"){

      if(verbose > 0) cat("Design", i, "encountered an error.\n")
      warning("Design ", i, " encountered an error.")

      out = list(
        X_orig = X,
        X = X,
        opt_crit_value_orig = Inf,
        opt_crit_value = Inf,
        n_iter = NA
      )
    }

    return(out)
  }, mc.cores = n_cores)

  # Get optimality values for all designs
  optimality_values = unlist(lapply(results, function(x) x$opt_crit_value))

  # Return the result with the best optimality criterion
  X_result = results[[which.min(optimality_values)]]

  # # Coordinate exchanges:
  # X_result = mixtureCoordinateExchangeGaussian(X, order, n_cox_points, max_it, verbose, opt_crit, W)

  out_list = list(
    X_orig = X_result$X_orig,
    X = X_result$X,
    opt_crit_value_orig = X_result$opt_crit_value_orig,
    opt_crit_value = X_result$opt_crit_value,
    n_iter = X_result$n_iter,
    opt_crit = ifelse(X_result$opt_crit == 0, "D-optimality", "I-optimality")
  )

  if(plot_designs) {
    if(q == 3) gaussian_plot_result(out_list)
    else warning("Could not plot results because q != 3")
  }
  return(out_list)

}




#' TODO: write doc
#' @export
gaussian_plot_result = function(res_alg){
  # res_alg: output of a call to mixture_coord_ex_gaussian() function

  ggtern::grid.arrange(
    res_alg$X_orig %>%
      dplyr::as_tibble() %>%
      purrr::set_names(c("c1", "c2", "c3")) %>%
      ggtern::ggtern(ggtern::aes(c1, c2, c3)) +
      geom_point(shape = "x", size = 4) +
      theme_minimal() +
      ggtern::theme_nomask() +
      ggtitle(
        label = paste0("Criterion: ", res_alg$opt_crit),
        subtitle = paste0("Value = ", round(res_alg$opt_crit_value_orig, 3)))
    ,
    res_alg$X %>%
      dplyr::as_tibble() %>%
      purrr::set_names(c("c1", "c2", "c3")) %>%
      ggtern::ggtern(ggtern::aes(c1, c2, c3)) +
      geom_point(shape = "x", size = 4) +
      theme_minimal() +
      ggtern::theme_nomask() +
      ggtitle(label = paste0("Criterion: ", res_alg$opt_crit),
              subtitle = paste0("Value = ", round(res_alg$opt_crit_value, 3)))
    ,
    ncol = 2
  )


}






#' TODO: write doc
#' @export
create_moment_matrix_gaussian = function(q){

  m = (q^3+ 5*q)/6

  f = lapply(1:m, function(x) rep(0, q))

  counter = 0
  # Fill indicators of first part of the model expansion
  for(i in 1:q){
    counter = counter + 1
    f[[counter]][i] = 1
  }

  # Fill indicators of second part of the model expansion
  for(i in 1:(q-1)){
    for(j in (i+1):q){
      counter = counter + 1
      f[[counter]][i] = 1
      f[[counter]][j] = 1
    }
  }


  # Fill indicators of third part of the model expansion
  for(i in 1:(q-2)){
    for(j in (i+1):(q-1)){
      for(k in (j+1):q){
        counter = counter + 1
        f[[counter]][i] = 1
        f[[counter]][j] = 1
        f[[counter]][k] = 1
      }
    }
  }


  W = matrix(rep(NA_real_, m*m), ncol = m)

  for(i in 1:m){
    for(j in 1:m){

      aux_ij = f[[i]] + f[[j]]
      num_ij = prod(factorial(aux_ij))
      denom_ij = factorial(2 + sum(aux_ij))
      W[i,j] = num_ij/denom_ij
    }
  }

  return(W)
}









#' TODO: write doc
#' @export
get_opt_crit_value_Gaussian = function(X, order = 1, opt_crit = 0){

  q = dim(X)[2]

  if(opt_crit == 0){
    # "D-optimality"
    W = matrix(0.0, nrow = 1)
  } else{
    # "I-optimality")
    W = create_moment_matrix_gaussian(q)
  }

  return(getOptCritValueGaussian(X = X, order = order, q = q, opt_crit = opt_crit, W = W))
}









# #' R wrapper for BrentCoxScheffeGaussian in C++.
# #' Minimizes the efficiency function for Scheffé model using Brent's method for local optima.
# #' The function also checks the edge cases (in normal mixtures it's 0 and 1) to see if the minimizers are there.
# #' TODO: write doc
# #' @export
# brent_cox_scheffe_gaussian = function(X, j, i, order, opt_crit,
#                                       lower = 0, upper = 1, tol = 0.0001){
#
#   q = dim(X)[2]
#
#   if(opt_crit == 0){
#     # "D-optimality"
#     W = matrix(0.0, nrow = 1)
#   } else{
#     # "I-optimality")
#     W = create_moment_matrix_gaussian(q)
#   }
#
#   return(
#     BrentCoxScheffeGaussian(X, j, i, order, opt_crit, W,
#                             lower,upper, tol)
#   )
# }
#
#
#
#
#
# #' R wrapper for BrentGloCoxScheffeGaussian in C++.
# #' Minimizes the efficiency function for Scheffé model using Brent's method for global optima.
# #' Needs a bound for the second derivative.
# #' TODO: write doc
# #' @export
# brent_global_cox_scheffe_gaussian = function(
#   X, j, i, order, opt_crit,
#   lower = 0,
#   upper = 1,
#   initial_guess = 0.5,
#   hessian_bound = 1e5,
#   abs_err_tol = 0.0001,
#   tol = 0.0001){
#
#   q = dim(X)[2]
#
#   if(opt_crit == 0){
#     # "D-optimality"
#     W = matrix(0.0, nrow = 1)
#   } else{
#     # "I-optimality")
#     W = create_moment_matrix_gaussian(q)
#   }
#
#   return(
#     BrentGloCoxScheffeGaussian(X, j, i, order, opt_crit, W,
#                                lower,
#                                upper,
#                                initial_guess,
#                                hessian_bound,
#                                abs_err_tol,
#                                tol)
#   )
# }




