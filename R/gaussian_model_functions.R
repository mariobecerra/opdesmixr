
#' Creation of a random initial design using the linear classical Gaussian model.
#'
#' @param n_runs Number of runs in the design.
#' @param q Number of mixture ingredients.
#' @param seed integer used for reproducibility
#'
#' @return Matrix of size \code{(n_runs, q)}.
#'
#' @examples
#' gaussian_create_random_initial_design(10, 3, 2, seed = 3)
#'
#'
#' @export
gaussian_create_random_initial_design = function(n_runs, q, seed = NULL){

  if(!is.null(seed)) set.seed(seed)

  # Mixture variables
  X = matrix(rep(NA_real_, n_runs*q), nrow = n_runs)

  for(i in 1:nrow(X)){
    rands = runif(q)

    # rows in X sum to 1:
    X[i,] = rands/sum(rands)
  }

  return(X)
}




#' Coordinate exchange algorithm for a mixture model assuming Gaussian iid errors.
#'
#' Performs the coordinate exchange algorithm for a Scheffé mixture model. It can have many different random starts, or the coordinate exchange algorithm can be performed in a user-supplied matrix.
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
#' @return List with the following elements: \itemize{
#'     \item \code{X_orig}: The original design. Matrix of size (n_runs, q).
#'     \item \code{X}: The optimized design. Matrix of size (n_runs, q).
#'     \item \code{opt_crit_value_orig}: efficiency of the original design.
#'     \item \code{opt_crit_value}: efficiency of the optimized design.
#'     \item \code{n_iter}: Number of iterations performed.
#'     \item \code{efficiency_value_per_iteration}:
#'     \item \code{opt_crit}: The optimality criterion used.
#'     \item \code{q}: Number of mixture ingredients.
#'     \item \code{seed}: seed used to generate the final design. If a design was used as input by the user, this will be NA.
#'  }
#'
#' @export
gaussian_mixture_coord_exch = function(
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
  n_cores = 1
){


  #############################################
  ## Check that the order is okay
  #############################################
  if(order != 1 & order != 2 & order != 3){
    stop("Inadmissible value for order. Must be 1, 2, or 3")
  }


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

  seeds_designs = NULL
  #############################################
  ## Create random initial designs or check that
  ## the provide design is okay.
  #############################################
  if(is.null(X)){
    if(!is.null(seed)) set.seed(seed)

    seeds_designs = sample.int(1e9, n_random_starts)

    designs = lapply(seeds_designs, function(x){
      des = gaussian_create_random_initial_design(n_runs = n_runs, q = q, seed = x)
      return(des)
    })
  } else{
    n_runs = nrow(X)
    q = ncol(X)

    # Check that rows in first q columns of X sum to 1
    row_sums = apply(X[, 1:q], 1, sum)

    if(sum(abs(row_sums - rep(1, n_runs)) < 1e-10) != n_runs){
      stop("Rows in first q columns of X must sum up to 1.")
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
    W = gaussian_create_moment_matrix(q = q, order = order)
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
        n_cox_points = n_cox_points,
        n_pv = 0 # This refers to process variables. They are not used in this version of the package, but the building blocks are already in the C++ code, so they must be declared here in R.
      ),
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

  out_list = list(
    X_orig = X_result$X_orig,
    X = X_result$X,
    opt_crit_value_orig = X_result$opt_crit_value_orig,
    opt_crit_value = X_result$opt_crit_value,
    n_iter = X_result$n_iter,
    efficiency_value_per_iteration = X_result$efficiency_value_per_iteration,
    opt_crit = ifelse(opt_crit == 0, "D-optimality", "I-optimality"),
    q = q,
    seed = seeds_designs[which.min(optimality_values)]
  )

  if(plot_designs) {
    if(q == 3) gaussian_plot_result(out_list)
    else warning("Could not plot results because q != 3")
  }
  return(out_list)

}




#' Plots the result of the list that results from the \code{gaussian_mixture_coord_exch()} function.
#'
#' @param res_alg List resulting from a call to the \code{gaussian_mixture_coord_exch()} function.
#'
#' @details
#' Only works for \code{q = 3}.
#'
#' @export
gaussian_plot_result = function(res_alg){

  stopifnot(res_alg$q == 3)

  ggtern::grid.arrange(
    res_alg$X_orig[, 1:res_alg$q] %>%
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
    res_alg$X[, 1:res_alg$q] %>%
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











#' Computes optimality criterion value using the linear classical Gaussian model.
#'
#' Computes optimality criterion value for a design matrix \code{X}.
#'
#' @param X Design matrix.
#' @param order integer corresponding to a Scheffé model order (1, 2, or 3).
#' @param opt_crit optimality criterion: 0 or "D" is D-optimality; while 1 or "I" is I-optimality.
#'
#' @return Returns the value of the optimality criterion for this particular design.
#'
#' @export
gaussian_get_opt_crit_value = function(X, order = 1, opt_crit = 0){

  #############################################
  ## Check that the order is okay
  #############################################
  if(order != 1 & order != 2 & order != 3){
    stop("Inadmissible value for order. Must be 1, 2, or 3")
  }


  q = dim(X)[2]

  #############################################
  ## Check that optimality criterion is okay
  #############################################
  if(!(opt_crit %in% c(0, 1) | opt_crit %in% c("D", "I"))){
    stop('Unknown optimality criterion. Must be either "D" or 0 for D-optimality, or "I" or 1 for I-optimality.' )
  }

  # Recode opt_crit
  if(opt_crit == "D") opt_crit = 0
  if(opt_crit == "I") opt_crit = 1




  if(opt_crit == 0){
    # "D-optimality"
    W = matrix(0.0, nrow = 1)
  } else{
    # "I-optimality")
    W = gaussian_create_moment_matrix(q = q, order = order)
  }

  # n_pv refers to process variables. They are not used in this version of the package, but the building blocks are already in the C++ code, so they must be declared here in R.
  return(getOptCritValueGaussian(X = X, order = order, q = q, opt_crit = opt_crit, W = W, n_pv = 0))
}










#' Computes moment matrix using the linear classical Gaussian model.
#'
#' @param q Number of mixture ingredients
#' @param order integer corresponding to a Scheffé model order (1, 2, or 3).
#'
#' @export
gaussian_create_moment_matrix = function(q, order = 3){

  stopifnot(order %in% 1:3)


  if(order == 1){
    m = q
  } else{
    if(order == 2){
      m = q*(q-1)/2 + q
    } else{
      if(order == 3){
        m = (q^3+ 5*q)/6 # = q + q*(q-1)/2 + q*(q-1)*(q-2)/6
      } else{
        stop()
      }
    }
  }

  f_matrix = matrix(rep(0L, m*(q)), ncol = q)

  counter = 0
  # Fill indicators of first part of the model expansion
  for(i in 1:q){
    counter = counter + 1
    f_matrix[counter, i] = 1
  }


  # Fill indicators of second part of the model expansion
  if(order >= 2){
    for(i in 1:(q-1)){
      for(j in (i+1):q){
        counter = counter + 1
        f_matrix[counter, i] = 1
        f_matrix[counter, j] = 1
      }
    }
  }


  # Fill indicators of third part of the model expansion
  if(order == 3){
    for(i in 1:(q-2)){
      for(j in (i+1):(q-1)){
        for(k in (j+1):q){
          counter = counter + 1
          f_matrix[counter, i] = 1
          f_matrix[counter, j] = 1
          f_matrix[counter, k] = 1
        }
      }
    }
  }





  W = matrix(rep(NA_real_, m*m), ncol = m)

  for(i in 1:m){
    for(j in 1:m){

      aux_ij_1 = f_matrix[i, 1:q] + f_matrix[j, 1:q]
      num_ij_1 = prod(factorial(aux_ij_1))
      denom_ij_1 = factorial(q - 1 + sum(aux_ij_1))


      aux_ij_2 = 0.0
      num_ij_2 = 1.0

      denom_ij_2 = prod(1 + aux_ij_2)


      W[i,j] = (num_ij_1 * num_ij_2)/(denom_ij_1 * denom_ij_2)

    }
  }

  return(W)
}


