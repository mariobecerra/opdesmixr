
#' TODO: write doc
#' @export
gaussian_create_random_initial_design = function(n_runs, q, n_pv = 0, pv_bounds = NULL, seed = NULL){
  # gaussian_create_random_initial_design(10, 3, 2, seed = 3)
  # gaussian_create_random_initial_design(10, 3, 2, pv_bounds = list(c(-1, 0), c(-1, 0)), 3)

  if(!is.null(seed)) set.seed(seed)

  # Mixture variables
  X = matrix(rep(NA_real_, n_runs*q), nrow = n_runs)

  for(i in 1:nrow(X)){
    rands = runif(q)

    # rows in X sum to 1:
    X[i,] = rands/sum(rands)
  }

  # Process variables
  if(n_pv > 0){
    if(is.null(pv_bounds)){
      warning("Number of process variables (n_pv = ", n_pv, ") provided but no information about the bounds. Using interval (0, 1) for all process variables.")
      pv_bounds = lapply(1:n_pv, function(.) return(c(0, 1)))
    } else{
      stopifnot(is.list(pv_bounds))
      for(k in 1:length(pv_bounds)){
        stopifnot(length(pv_bounds[[k]]) == 2)
      }
      stopifnot(length(pv_bounds) == n_pv)
    }

    # Create matrix and randomly fill them
    X_pv = matrix(rep(NA_real_, n_runs*n_pv), ncol = n_pv)
    for(k in 1:n_pv){
      X_pv[,k] = runif(n = n_runs, min = pv_bounds[[k]][1], max = pv_bounds[[k]][2])
    }
    X = cbind(X, X_pv)

  }

  return(X)
}




#' Coordinate exchange algorithm for a mixture model assuming Gaussian iid errors.
#'
#' \code{gaussian_mixture_coord_exch} Performs the coordinate exchange algorithm for a Scheffé mixture model. It can have many different random starts, or the coordinate exchange algorithm can be performed in a user-supplied matrix.
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
  n_cores = 1,
  n_pv = 0,
  pv_bounds = NULL
  ){



  #############################################
  ## Check that the order is okay
  #############################################
  if(order != 1 & order != 2 & order != 3 & order != 4){
    stop("Inadmissible value for order. Must be 1, 2, 3 or 4")
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
      des = gaussian_create_random_initial_design(n_runs = n_runs, q = q, seed = x, n_pv = n_pv, pv_bounds = pv_bounds)
      return(des)
    })
  } else{
    n_runs = nrow(X)
    q = ncol(X) - n_pv

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
    W = gaussian_create_moment_matrix(q = q, n_pv = n_pv, order = order)
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
        n_pv = n_pv
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
    n_pv = n_pv,
    seed = seeds_designs[which.min(optimality_values)]
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
  # res_alg: output of a call to gaussian_mixture_coord_exch() function

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











#' TODO: write doc
#' @export
gaussian_get_opt_crit_value = function(X, order = 1, opt_crit = 0, n_pv = 0){



  #############################################
  ## Check that the order is okay
  #############################################
  if(order != 1 & order != 2 & order != 3 & order != 4){
    stop("Inadmissible value for order. Must be 1, 2, 3 or 4")
  }


  q = dim(X)[2] - n_pv

  if(order == 4 & n_pv == 0) stop("If order == 4 then n_pv must be greater than 0")
  if(order %in% 1:3 & n_pv > 0) stop("If order is 1, 2, or 3 then n_pv must be 0")


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
    W = gaussian_create_moment_matrix(q, order)
  }

  return(getOptCritValueGaussian(X = X, order = order, q = q, opt_crit = opt_crit, W = W, n_pv = n_pv))
}










#' TODO: write doc
#' @export
gaussian_create_moment_matrix = function(q, n_pv = 0, order = 3){

  stopifnot(order %in% 1:4)

  if(order == 4 & n_pv == 0) stop("If order == 4 then n_pv must be greater than 0")
  if(order %in% 1:3 & n_pv > 0) stop("If order is 1, 2, or 3 then n_pv must be 0")

  if(order == 1){
    m = q
  } else{
    if(order == 2){
      m = q*(q-1)/2 + q
    } else{
      if(order == 3){
        m = (q^3+ 5*q)/6 # = q + q*(q-1)/2 + q*(q-1)*(q-2)/6
      } else{
        m = q + q*(q-1)/2 + q*n_pv + n_pv*(n_pv-1)/2 + n_pv
      }
    }
  }

  f_matrix = matrix(rep(0L, m*(q + n_pv)), ncol = q + n_pv)

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



  # Fill indicators of fourth part of the model expansion when n_pv > 0
  if(order == 4){

    for(i in 1:n_pv){
      for(k in 1:q){
        counter = counter + 1
        f_matrix[counter, q + i] = 1
        f_matrix[counter, k] = 1
      }
    }

    if(n_pv > 1){
      for(i in 1:(n_pv-1)){
        for(j in (i+1):n_pv){
          counter = counter + 1
          f_matrix[counter, q + i] = 1
          f_matrix[counter, q + j] = 1
        }
      }
    }


    for(i in 1:n_pv){
      counter = counter + 1
      f_matrix[counter, q + i] = 2
    }

  }


  W = matrix(rep(NA_real_, m*m), ncol = m)

  for(i in 1:m){
    for(j in 1:m){

      aux_ij_1 = f_matrix[i, 1:q] + f_matrix[j, 1:q]
      if(sum(aux_ij_1 > 0)){
        num_ij_1 = prod(factorial(aux_ij_1))
        denom_ij_1 = factorial(q - 1 + sum(aux_ij_1))
      } else{
        num_ij_1 = 1
        denom_ij_1 =  factorial(q - 1)
      }

      denom_ij_2 = 1
      if(n_pv > 0){
        aux_ij_2 = f_matrix[i, (q+1):(q+n_pv)] + f_matrix[j, (q+1):(q+n_pv)]
        if(sum(aux_ij_2 > 0)){
          denom_ij_2 = prod(1 + aux_ij_2)
        }
      }

      W[i,j] = (num_ij_1/denom_ij_1)*(1/denom_ij_2)
    }
  }

  return(W)
}
