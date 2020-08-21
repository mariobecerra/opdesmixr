
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




#' TODO: write doc
#' @export
mixture_coord_ex_gaussian = function(
  n_runs = NULL,
  q = NULL,
  n_random_starts = 100,
  X = NULL,
  order = 1,
  n_cox_points = 30,
  max_it = 10,
  plot_designs = F,
  verbose = 1,
  opt_crit = 0,
  seed = NULL,
  n_cores = 1){



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



  # If criterion is D-optimality, send a matrix with only one zero element as a moment matrix
  # Maybe a NULL value would be better. Gotta check.
  if(opt_crit == 0){
    # "D-optimality"
    W = matrix(0.0, nrow = 1)
  } else{
    # "I-optimality")
    W = create_moment_matrix_gaussian(q)
  }

  if(.Platform$OS.type != "unix") n_cores = 1

  # Apply the coordinate exchange algorithm to all the designs generated
  results = parallel::mclapply(seq_along(designs), function(i){
    X = designs[[i]]

    if(verbose > 0) cat("\nDesign", i, "\n")

    out = try(
      mixtureCoordinateExchangeGaussian(
        X, order, n_cox_points, max_it, verbose, opt_crit, W
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
      ggtern(aes(c1, c2, c3)) +
      geom_point(shape = "x", size = 4) +
      theme_minimal() +
      theme_nomask() +
      ggtitle(
        label = paste0("Criterion: ", res_alg$opt_crit),
        subtitle = paste0("Value = ", round(res_alg$opt_crit_value_orig, 3)))
    ,
    res_alg$X %>%
      dplyr::as_tibble() %>%
      purrr::set_names(c("c1", "c2", "c3")) %>%
      ggtern(aes(c1, c2, c3)) +
      geom_point(shape = "x", size = 4) +
      theme_minimal() +
      theme_nomask() +
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















#' Function that computes optimality criterion value for MNL model
#'
#' \code{get_opt_crit_value_MNL} computes optimality criterion value for MNL model
#' @param X 3 dimensional array with dimensions (q, J, S) where:
#'     q is the number of ingredient proportions,
#'     J is the number of alternatives within a choice set,
#'     S is the number of choice sets.
#' @param opt_crit optimality criterion: 0 is D-optimality and 1 is I-optimality
#' @return Returns the value of the optimality criterion for this particular design and this beta vector
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
