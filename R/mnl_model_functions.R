mnl_check_order = function(order, n_pv){
  # Checks whether the order is 1, 2, or 3. And that order is 2 in case n_pv > 0.
  if(!(order %in% 1:3)){
    stop("Inadmissible value for 'order' parameter. Must be 1, 2, or 3")
  }

  if(n_pv > 0 & order != 2) stop("If n_pv > 0, then 'order' must be 2.")
}





#' Get number of parameters.
#'
#' Returns an integer with the number of parameters depending on \code{q}, \code{n_pv}, and \code{order}.
#'
#' @param q integer specifying the number of mixture ingredients
#' @param n_pv Number of process variables
#' @param order integer corresponding to a Scheffé model order (1, 2, 3 if no process variables are used; 2 if process variables are used).
#'
#' @export
mnl_get_number_of_parameters = function(q, n_pv, order){

  # Check order
  mnl_check_order(order = order, n_pv = n_pv)

  if(order == 1){
    m = q
  } else{
    if(order == 2){
      m = q + (q-1)*q/2 + q*n_pv + n_pv*(n_pv-1)/2 + n_pv
    } else{
      if(order == 3){
        m = (q^3+ 5*q)/6 # = q + q*(q-1)/2 + q*(q-1)*(q-2)/6
      }  else{
        stop("Wrong order")
        # This condition won't be met because of the checking in the previous lines.
        # I leave it just for better legibility
      }
    }
  }

  return(m)
}




#' Creation of a random initial design using the MNL model
#'
#' Creates a random initial design using the MNL model of the specified dimensions.
#'
#' @param q integer specifying the number of mixture ingredients
#' @param J integer specifying the number of alternatives within each choice set
#' @param S integer specifying the number of choice sets
#' @param n_pv Number of process variables
#' @param pv_bounds Vector or matrix that specifies the bounds of the process variables.
#' @param seed integer used for reproducibility
#'
#' @return 3-dimensional array of size \code{(q, J, S)}.
#'
#' @examples
#' mnl_create_random_initial_design(q = 3, J = 2, S = 8, seed = 2020)
#'
#' @details
#' \code{pv_bounds} must be either: \itemize{
#'     \item NULL (for default values), or
#'     \item a vector with 2 scalars (the second element must be greater than the first), or
#'     \item a matrix with n_pv rows and 2 columns, i.e., each row represents a bound for each process variable.
#' }
#'
#' @export
mnl_create_random_initial_design = function(q, J, S, n_pv = 0,  pv_bounds = NULL, seed = NULL){
  X = array(rep(NA_real_, (q+n_pv)*J*S), dim = c((q+n_pv), J, S))

  if(!is.null(seed)) set.seed(seed)

  for(j in 1:J){
    for(s in 1:S){
      rands = runif(q)
      # mixture ingredients must sum up to 1
      X[1:q,j, s] = rands/sum(rands)
    }
  }

  if(n_pv > 0){
    if(is.null(pv_bounds)){
      warning("Number of process variables (n_pv = ", n_pv, ") provided but no information about the bounds. Using interval [-1, 1] for all process variables.")
      pv_bounds = t(sapply(1:n_pv, function(.) return(c(-1, 1))))
    } else{

      if(is.vector(pv_bounds) & !is.list(pv_bounds) & !is.matrix(pv_bounds)){ # If only one interval is given
        stopifnot(length(pv_bounds) == 2) # Make sure it's an interval of the form [a, b]
        stopifnot(pv_bounds[2] > pv_bounds[1])
        pv_bounds = t(sapply(1:n_pv, function(.) return(pv_bounds)))
      }

      stopifnot(is.matrix(pv_bounds))
      stopifnot(ncol(pv_bounds) == 2)
      stopifnot(nrow(pv_bounds) == n_pv)
      for(i in 1:nrow(pv_bounds)){
        stopifnot(pv_bounds[i, 2] > pv_bounds[i, 1])
      }
    }

    # Fill rest of array with process variables
    for(l in 1:n_pv){
      X[l + q, , ] = runif(n = J*S, min = pv_bounds[l, 1], max = pv_bounds[l, 2])
    }

  }

  return(X)
}



#' Computes optimality criterion value using the MNL model
#'
#' Computes optimality criterion value for a design array \code{X} and \code{beta} parameter.
#'
#' @param X 3 dimensional array of size \code{(q, J, S)} where: \itemize{
#'     \item q is the number of ingredient proportions,
#'     \item J is the number of alternatives within a choice set,
#'     \item S is the number of choice sets.
#' }
#' @param beta numeric vector containing the parameters or numeric matrix containing draws of the prior distribution of the parameters.
#' @param order integer corresponding to a Scheffé model order (1, 2, 3 if no process variables are used; 2 if process variables are used).
#' @param opt_crit optimality criterion: 0 or "D" is D-optimality; while 1 or "I" is I-optimality.
#' @param transform_beta boolean parameter. Should the beta vector/matrix be transformed by subtracting the \code{q}-th element?
#' @param n_pv Number of process variables.
#'
#' @return Returns the value of the optimality criterion for this particular design and this \code{beta} vector
#' @export
mnl_get_opt_crit_value = function(X, beta, order, opt_crit = "D", transform_beta = T, n_pv = 0){

  # Recode opt_crit
  if(opt_crit == "D") opt_crit = 0
  if(opt_crit == "I") opt_crit = 1

  q = dim(X)[1] - n_pv

  ## Get m (No need to check order because this function checks it).
  m = mnl_get_number_of_parameters(q = q, order = order, n_pv = n_pv)

  #############################################
  ## Check beta dimensions and transform if necessary
  #############################################

  if(is.vector(beta)) {
    if(m != length(beta) & transform_beta) stop("Incompatible size in beta and q: beta must be of length ", m)
    if(m-1 != length(beta) & !transform_beta) stop("Incompatible size in beta and q: beta must be of length ", m-1)
  }

  if(is.matrix(beta)) {
    if(m != ncol(beta) & transform_beta) stop("Incompatible size in beta and q: beta must have ", m,  " columns")
    if(m-1 != ncol(beta) & !transform_beta) stop("Incompatible size in beta and q: beta must have ", m-1,  " columns")
  }


  if(is.vector(beta)) beta_mat = matrix(beta, nrow = 1)
  else beta_mat = beta

  #############################################
  ## Get moment matrix
  #############################################

  if(opt_crit == 0){
    # "D-optimality"
    W = matrix(0.0, nrow = 1)
  } else{
    # "I-optimality")
    W = mnl_create_moment_matrix(q = q, order = order, n_pv = n_pv, pv_bounds = c(-1, 1))
  }

  #############################################
  ## Get optimality criterion value
  #############################################

  return(getOptCritValueMNL(X = X, beta = beta_mat, opt_crit = opt_crit, verbose = 0, W = W, order = order, transform_beta = transform_beta, n_pv = n_pv))

}



#' Creation of a random parameter vector using the MNL model
#'
#' Creates a random parameter vector with the appropriate dimensions as in \emph{Bayesian D-optimal choice designs for mixtures} by Ruseckaite, Goos & Fok (2017).
#'
#' @param q integer specifying the number of mixture ingredients
#' @param n_pv number of process variables
#' @param order integer corresponding to a Scheffé model order (1, 2, 3 if no process variables are used; 2 if process variables are used).
#' @param seed integer used for reproducibility
#'
#' @return Returns a list in which the first element of the list is
#'     a numerical vector with the parameters and the second element is a matrix with
#'     the indices as in \emph{Bayesian D-optimal choice designs for mixtures} by Ruseckaite, Goos & Fok (2017).
#'
#' @examples
#' mnl_create_random_beta(q = 3, order = 3, seed = 3)
#'
#' @export
mnl_create_random_beta = function(q, n_pv = 0, order = 3, seed = NULL){

  if(!all.equal(q, floor(q))) stop("q does not seem to an integer.")
  stopifnot(q > 0)
  if(!all.equal(n_pv, floor(n_pv))) stop("n_pv does not seem to an integer.")
  stopifnot(n_pv >= 0)

  ## Get m (No need to check order because this function checks it).
  m = mnl_get_number_of_parameters(q = q, order = order, n_pv = n_pv)

  if(!is.null(seed)) set.seed(seed)

  beta = rnorm(m)

  # Return a list for backwards compatibility, but eventually I'll have to change it to return only a vector
  return(list(beta = beta))
}


#' Coordinate exchange algorithm for a Multinomial Logit Scheffé model.
#'
#' Performs the coordinate exchange algorithm for a Multinomial Logit Scheffé model as described in
#' \emph{Bayesian D-optimal choice designs for mixtures} by Ruseckaite, Goos & Fok (2017).
#'
#' @param q number of mixture ingredient proportions.
#' @param J number of alternatives within a choice set.
#' @param S number of choice sets.
#' @param n_random_starts number or random starts. Defaults to 100.
#' @param X If an initial design is to be supplied, then it must be a 3-dimensional array of size \code{(q, J, S)}, with q, J, and S defined above.
#' @param beta Prior parameters. For a locally optimal design, it should be a numeric vector of length m = (q^3 + 5*q)/6. For a pseudo-Bayesian design, it must be a matrix with prior simulations of size (nxm) where m is previously defined and m is the number of prior draws, i.e., there is a prior draw per row.
#' @param transform_beta boolean parameter. Should the \code{beta} vector/matrix be transformed by subtracting the q-th element?
#' @param order integer corresponding to a Scheffé model order (1, 2, 3 if no process variables are used; 2 if process variables are used).
#' @param opt_method Optimization method in each step of the coordinate exchange algorithm.
#'      It can be "B" (Brent's algorithm) or "D" (discretization of Cox direction)
#' @param max_it integer for maximum number of iterations that the coordinate exchange algorithm will do
#' @param tol A positive error tolerance in Brent's method.
#' @param n_cox_points number of points to use in the discretization of Cox direction. Ignored if opt_method is Brent.
#' @param plot_designs boolean. If TRUE, shows a plot of the initial and the final design. Only works if q is 3 or 4.
#' @param verbose level of verbosity. See below for details.
#' @param opt_crit optimality criterion: D-optimality ("D" or 0) or I-optimality ("I" or 1).
#' @param seed Seed for reproducibility.
#' @param n_cores Number of cores for parallel processing.
#' @param save_all_designs Whether the function should return a list with all the designs created at random or only the best.
#' @param n_pv Number of process variables
#'
#' @return The function returns a list with 11 elements: \enumerate{
#'     \item \code{X_orig}: The original design. A 3-dimensional array of size \code{(q, J, S)}.
#'     \item \code{X}: The optimized design. A 3-dimensional array of size \code{(q, J, S)}.
#'     \item \code{beta}: The original \code{beta} vector or matrix.
#'     \item \code{opt_crit_value_orig}: efficiency of the original design.
#'     \item \code{opt_crit_value}: efficiency of the optimized design.
#'     \item \code{n_iter}: Number of iterations performed.
#'     \item \code{efficiency_value_per_iteration}: Efficiency value in each iteration of the algorithm.
#'     \item \code{opt_crit}: The optimality criterion used.
#'     \item \code{q}: Number of mixture ingredients.
#'     \item \code{n_pv}: Number of process variables.
#'     \item \code{seed}: seed used to generate the final design. If a design was used as input by the user, this will be NA.
#'  }
#'
#' @details
#'
#' Verbosity levels: each level prints the previous plus additional things:
#' \enumerate{
#'     \item Print the efficiency value in each iteration and a final summary
#'     \item Print the values of k, s, i, and efficiency value in each subiteration
#'     \item Print the resulting X after each iteration, i.e., after each complete pass on the data
#'     \item Print efficiency value for each point in the Cox direction discretization
#'     \item Print the resulting X and information matrix after each subiteration
#'     \item Print the resulting X or each point in the Cox direction discretization
#'  }
#'

#' @export
mnl_mixture_coord_exch = function(
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
  tol = 0.0001,
  n_cox_points = NULL,
  plot_designs = F,
  verbose = 1,
  opt_crit = 0,
  seed = NULL,
  n_cores = 1,
  save_all_designs = F,
  n_pv = 0
){
  # Assumes process variables are bwteen -1 and 1

  t1 = Sys.time()
  if(verbose >= 1) cat("Starts at", substr(as.character(t1), 12, 19), "\n")


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


  #############################################
  ## Create random initial designs or check that
  ## the provide design is okay. Also check betas.
  #############################################


  if(is.null(X)){

    if(!is.null(seed)) set.seed(seed)

    seeds_designs = sample.int(1e9, n_random_starts)

    designs = lapply(seeds_designs, function(x){
      des = mnl_create_random_initial_design(q, J, S, seed = x, n_pv = n_pv, pv_bounds = c(-1, 1))
      return(des)
    })

  } else{
    # Some input checks
    dim_X = dim(X)

    if(length(dim_X) != 3) stop("X must be a 3 dimensional array.")
    if(!(is.vector(beta) | is.matrix(beta))) stop("beta is not a vector or a matrix. It must be a numerical or integer vector or matrix.")

    seeds_designs = NA

    q = dim_X[1] - n_pv

    designs = list(X)

    n_random_starts = 0
  }



  if(!all.equal(q, floor(q))) stop("q does not seem to an integer.")
  stopifnot(q > 0)
  if(!all.equal(n_pv, floor(n_pv))) stop("n_pv does not seem to an integer.")
  stopifnot(n_pv >= 0)

  ## Get m (No need to check order because this function checks it).
  m = mnl_get_number_of_parameters(q = q, order = order, n_pv = n_pv)


  if(is.vector(beta)) {
    if(m != length(beta) & transform_beta) stop("Incompatible size in beta and q: beta must be of length ", m)
    if(m-1 != length(beta) & !transform_beta) stop("Incompatible size in beta and q: beta must be of length ", m-1)
  }

  if(is.matrix(beta)) {
    if(m != ncol(beta) & transform_beta) stop("Incompatible size in beta and q: beta must have ", m,  " columns")
    if(m-1 != ncol(beta) & !transform_beta) stop("Incompatible size in beta and q: beta must have ", m-1,  " columns")
  }



  if(is.vector(beta)) {
    beta_mat = matrix(beta, nrow = 1)
  } else {
    beta_mat = beta
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
    W = mnl_create_moment_matrix(q = q, order = order, n_pv = n_pv, pv_bounds = c(-1, 1))
  }

  #############################################
  ## Check operating system for parallel processing
  #############################################

  if(.Platform$OS.type != "unix") n_cores = 1

  #############################################
  ## Apply the coordinate exchange algorithm to all the created designs
  #############################################

  results = parallel::mclapply(seq_along(designs), function(i){
    X = designs[[i]]

    if(verbose > 0) cat("\nDesign", i, "\n")

    order_cpp = ifelse(n_pv > 0, 4, order)

    out = try(
      mixtureCoordinateExchangeMNL(
        X_orig = X,
        beta = beta_mat,
        order = order_cpp,
        max_it = max_it,
        verbose = verbose,
        opt_crit = opt_crit,
        W = W,
        opt_method = opt_method,
        lower = 0,
        upper = 1,
        tol = tol,
        n_cox_points = n_cox_points,
        transform_beta = transform_beta,
        n_pv = n_pv
      ), silent = T)

    # If there was an error with this design, return whatever
    if(class(out) == "try-error"){

      if(verbose > 0) cat("Design", i, "encountered the following error:", out, "\n")
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


  if(save_all_designs & n_random_starts > 1){
    # Indices of optimality values in increasing order
    optimality_values_indices = order(optimality_values, decreasing = F)

    out_list = lapply(optimality_values_indices, function(i){
      X_result_i = results[[i]]
      out = list(
        X_orig = X_result_i$X_orig,
        X = X_result_i$X,
        beta = beta,
        opt_crit_value_orig = X_result_i$opt_crit_value_orig,
        opt_crit_value = X_result_i$opt_crit_value,
        n_iter = X_result_i$n_iter,
        efficiency_value_per_iteration = X_result_i$efficiency_value_per_iteration,
        opt_crit = ifelse(opt_crit == 0, "D-optimality", "I-optimality"),
        q = q,
        n_pv = n_pv,
        seed = seeds_designs[i]
      )
    })

  } else{
    # Find the index of the result with the best optimality criterion
    best_index = which.min(optimality_values)

    # Keep only the result with the best optimality criterion
    X_result = results[[best_index]]

    out_list = list(
      X_orig = X_result$X_orig,
      X = X_result$X,
      beta = beta,
      opt_crit_value_orig = X_result$opt_crit_value_orig,
      opt_crit_value = X_result$opt_crit_value,
      n_iter = X_result$n_iter,
      efficiency_value_per_iteration = X_result$efficiency_value_per_iteration,
      opt_crit = ifelse(opt_crit == 0, "D-optimality", "I-optimality"),
      q = q,
      n_pv = n_pv,
      seed = seeds_designs[best_index]
    )

  }

  if(plot_designs) {
    if(q == 3 | q == 4){
      if(!save_all_designs){
        mnl_plot_result(out_list)
      } else{
        warning("Plotted best design")
        mnl_plot_result(out_list[[1]])
      }
    } else {
      warning("Could not plot results because q is not 3 or 4.")
    }
  }

  t2 = Sys.time()

  if(verbose >= 1) cat("Ends at", substr(as.character(t2), 12, 19), "\n")
  if(verbose >= 1) cat("Time:", round(as.numeric(difftime(t2, t1, units = "secs")), 1), "seconds.\n")

  return(out_list)
}








plot_simplex_4d = function(X, phi = 40, theta = 140, cex_points = 0.8, cex_names = 0.8, ...){

  # Compute tetrahedron coordinates according to https://mathoverflow.net/a/184585
  tetra <- qr.Q(qr(matrix(1, nrow = 4)) ,complete = TRUE)[,-1]

  # Convert barycentric coordinates (4D) to cartesian coordinates (3D)
  X_3D <- geometry::bary2cart(tetra, X)

  # Plot data

  plot3D::scatter3D(X_3D[,1], X_3D[,2], X_3D[,3],
                    xlim = range(tetra[,1]), ylim = range(tetra[,2]), zlim = range(tetra[,3]),
                    col = "blue", pch = 16, box = FALSE, theta = theta, phi = phi, cex = cex_points, ...)
  plot3D::lines3D(tetra[c(1,2,3,4,1,3,1,2,4),1],
                  tetra[c(1,2,3,4,1,3,1,2,4),2],
                  tetra[c(1,2,3,4,1,3,1,2,4),3],
                  col = "grey", add = TRUE)
  if(!is.null(colnames(X))){
    plot3D::text3D(tetra[,1], tetra[,2], tetra[,3],
                   colnames(X), add = TRUE, cex = cex_names)
  }

}




mnl_plot_result_3d = function(
  res_alg,
  design = "final",
  phi = 40, theta = 140, cex_points = 0.5, cex_names = 0.5, ...){
  # res_alg: output of a call to mnl_mixture_coord_exch() function.
  # design: Plot the original or the final design?
  # It must be a design of 4 ingredients.

  if(design == "final"){
    i = 2
  } else{
    if(design == "original"){
      i = 1
    } else{
      stop('Must write either "final" or "original".')
    }
  }

  dim_X = dim(res_alg$X_orig)
  q = res_alg$q
  S = dim_X[3]

  if(q != 4) stop("Design must be of 4 ingredients.")
  if(res_alg$n_pv > 0) warning("There are process variables that are being ignored for the plots")

  # Convert 3 dimensional arrays into matrices by vertically binding them
  X_mat = t(res_alg[[i]][1:q,,1])
  for(s in 2:S){
    X_mat = rbind(X_mat, t(res_alg[[i]][,,s]))
  }

  plot_simplex_4d(X_mat,  phi = phi, theta = theta, cex_points = cex_points, cex_names = cex_names, ...)

}





mnl_plot_result_2d = function(res_alg){
  # res_alg: output of a call to mnl_mixture_coord_exch() function.
  # It must be a design of 3 ingredients.

  dim_X = dim(res_alg$X_orig)
  q = res_alg$q
  S = dim_X[3]

  if(q != 3) stop("Design must be of 3 ingredients.")
  if(res_alg$n_pv > 0) warning("There are process variables that are being ignored for the plots")

  # Convert 3 dimensional arrays into matrices by vertically binding them
  X_orig_mat = t(res_alg$X_orig[1:q,,1])
  for(s in 2:S){
    X_orig_mat = rbind(X_orig_mat, t(res_alg$X_orig[1:q,,s]))
  }

  X_final_mat = t(res_alg$X[1:q,,1])
  for(s in 2:S){
    X_final_mat = rbind(X_final_mat, t(res_alg$X[1:q,,s]))
  }

  # Plot matrices
  ggtern::grid.arrange(
    X_orig_mat %>%
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
    X_final_mat %>%
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







#' Plots the result of the list that results from the \code{mnl_mixture_coord_exch()} function.
#'
#' @param res_alg List resulting from a call to the \code{mnl_mixture_coord_exch()} function.
#'
#' @details
#' Only works for \code{q = 3} and \code{q = 4}.
#'
#' @export
mnl_plot_result = function(res_alg, ...){
  dim_X = dim(res_alg$X_orig)
  q = res_alg$q

  if(!(q == 3 | q == 4)) stop("Design must be of 3 or 4 ingredients.")

  if(q == 3){
    mnl_plot_result_2d(res_alg = res_alg)
  } else{
    mnl_plot_result_3d(res_alg = res_alg, ...)
  }

}





#' Computes moment matrix using the MNL model.
#'
#' @param q Number of mixture ingredients.
#' @param n_pv Number of process variables.
#' @param order integer corresponding to a Scheffé model order (1, 2, 3 if no process variables are used; 2 if process variables are used).
#' @param pv_bounds Vector or matrix that specifies the bounds of the process variables.
#'
#' @details
#' \code{pv_bounds} must be either: \itemize{
#'     \item NULL (for default values), or
#'     \item a vector with 2 scalars (the second element must be greater than the first), or
#'     \item a matrix with n_pv rows and 2 columns, i.e., each row represents a bound for each process variable.
#' }
#'
#' @return Matrix of size \code{(m, m)}.
#'
#' @export
mnl_create_moment_matrix = function(q, n_pv = 0, order = 3, pv_bounds = NULL){

  ## Get m (No need to check order because this function checks it).
  m = mnl_get_number_of_parameters(q = q, order = order, n_pv = n_pv)

  if(n_pv > 0){
    # Input checks for process variables
    if(is.null(pv_bounds)){
      warning("Number of process variables (n_pv = ", n_pv, ") provided but no information about the bounds. Using interval [-1, 1] for all process variables.")
      pv_bounds = t(sapply(1:n_pv, function(.) return(c(-1, 1))))
    } else{

      if(is.vector(pv_bounds) & !is.list(pv_bounds) & !is.matrix(pv_bounds)){ # If only one interval is given
        stopifnot(length(pv_bounds) == 2) # Make sure it's an interval of the form [a, b]
        stopifnot(pv_bounds[2] > pv_bounds[1])
        pv_bounds = t(sapply(1:n_pv, function(.) return(pv_bounds)))
      }

      stopifnot(is.matrix(pv_bounds))
      stopifnot(ncol(pv_bounds) == 2)
      stopifnot(nrow(pv_bounds) == n_pv)
      for(i in 1:nrow(pv_bounds)){
        stopifnot(pv_bounds[i, 2] > pv_bounds[i, 1])
      }
    }
  }


  f_matrix = matrix(rep(0L, (m-1)*(q + n_pv)), ncol = q + n_pv)

  counter = 0
  # Fill indicators of first part of the model expansion
  for(i in 1:(q-1)){
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
  if(n_pv > 0){

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


  W = matrix(rep(NA_real_, (m-1)^2), ncol = m-1)

  for(i in 1:(m-1)){
    for(j in 1:(m-1)){



      aux_ij_1 = f_matrix[i, 1:q] + f_matrix[j, 1:q]
      num_ij_1 = prod(factorial(aux_ij_1))
      denom_ij_1 = factorial(q - 1 + sum(aux_ij_1))


      if(n_pv > 0){
        aux_ij_2 = f_matrix[i, (q+1):(q+n_pv)] + f_matrix[j, (q+1):(q+n_pv)]
        num_ij_2 = 1.0
        for(l in 1:n_pv){
          num_ij_2 = num_ij_2 * (pv_bounds[l,2]^(aux_ij_2[l] + 1) - pv_bounds[l,1]^(aux_ij_2[l] + 1))
        }

      } else{
        aux_ij_2 = 0.0
        num_ij_2 = 1.0
      }
      denom_ij_2 = prod(1 + aux_ij_2)

      # # Same as:
      # if(n_pv > 0){
      #   aux_ij_2 = f_matrix[i, (q+1):(q+n_pv)] + f_matrix[j, (q+1):(q+n_pv)]
      #   num_ij_2 = prod(pv_bounds[,2]^(aux_ij_2 + 1) - pv_bounds[,1]^(aux_ij_2 + 1))
      # } else{
      #   aux_ij_2 = 0
      #   num_ij_2 = 1
      # }
      # denom_ij_2 = prod(1 + aux_ij_2)
      # # This code is faster, but it may be harder to read. In this case we don't care about speed.



      W[i,j] = (num_ij_1 * num_ij_2)/(denom_ij_1 * denom_ij_2)

    }
  }

  return(W)

}











#' Returns a 3 dimensional array to be used in the functions in the package
#'
#' @param des_df Dataframe with \code{(k+1)} columns where the first k columns are the variables and the \code{(k+1)}-th column is the choice set.
#'
#' @return 3-dimensional array of size \code{(q, J, S)}, with \itemize{
#'    \item q number of mixture ingredient proportions.
#'    \item J number of alternatives within a choice set.
#'    \item S number of choice sets.
#' }
#'
#' @export
mnl_design_dataframe_to_array = function(des_df){
  q = ncol(des_df) - 1
  JS = nrow(des_df)
  S = length(unique(des_df$choice_set))
  J = JS/S

  des_array = array(rep(NA_real_, JS*q), dim = c(q, J, S))

  s = 1
  for(cs in unique(des_df$choice_set)){
    des_array[, , s] = t(as.matrix(des_df[which(des_df$choice_set == cs), 1:q]))
    s = s + 1
  }

  return(des_array)
}





#' Convert a design array to tibble.
#'
#' This function allows the user to convert a 3-dimensional design array into a tibble, which helps for readability and to plot.
#'
#' @param des_array must be a 3-dimensional design array of dimension (k, J, S) where k is the number of variables, J the number of alternatives in each choice set and S is the number of choice sets.
#' @param names Names of the \code{(k+1)} columns. If NULL, then the names of the columns are \code{c(paste0("c", 1:k), "choice_set")}.
#'
#' @return Returns a tibble with \code{(k+1)} columns where the first \code{k} columns are the variables and the \code{(k+1)}-th column is the choice set.
#'
#' @examples
#' design_array_to_dataframe(mnl_create_random_initial_design(q = 3, J = 2, S = 4)
#'
#' design_array_to_dataframe(mnl_create_random_initial_design(q = 3, J = 2, S = 4), names = c("v1", "v2", "v3", "cs"))
#'
#' @export
mnl_design_array_to_dataframe = function(des_array, names = NULL){
  dim_X = dim(des_array)
  k = dim_X[1]
  S = dim_X[3]

  if(is.null(names)) names = c(paste0("c", 1:k), "choice_set")

  X_final_tbl = lapply(1:S, function(s){
    t(des_array[,,s]) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(choice_set = as.character(s))
  }) %>%
    dplyr::bind_rows() %>%
    purrr::set_names(names)

  return(X_final_tbl)

}






#' Get the information matrix using the MNL model.
#'
#' Computes the information matrix for a given design \code{X} and parameter vector \code{beta}.
#'
#' @param X 3-dimensional design array
#' @param beta Parameter vector
#' @param order integer corresponding to a Scheffé model order (1, 2, 3 if no process variables are used; 2 if process variables are used).
#' @param transform_beta boolean parameter. Should the beta vector/matrix be transformed by subtracting the \code{q}-th element?
#' @param n_pv Number of process variables
#'
#' @return Information matrix
#'
#' @export
mnl_get_information_matrix = function(X, beta, order, transform_beta, n_pv = 0){

  mnl_check_order(order = order, n_pv = n_pv)

  IM = getInformationMatrixMNL(
    X = X,
    beta = beta,
    order = order,
    transform_beta = transform_beta,
    n_pv = n_pv)

  return(IM)
}


#' Compute design matrix of a specific choice set using the MNL model.
#'
#' Compute the design matrix of choice set \code{s} using the MNL model.
#'
#' @param X 3-dimensional design cube/array of size \code{(q, J, S)}.
#' @param s integer corresponding to a choice set in 1 to S.
#' @param order integer corresponding to a Scheffé model order (1, 2, 3 if no process variables are used; 2 if process variables are used).
#' @param n_pv number of process variables. Defaults to 0.
#'
#' @return Matrix of dimension (J, m-1) that corresponds to the design matrix of choice set s.
#'
#' @export
mnl_get_Xs = function(X, s, order, n_pv = 0){
  mnl_check_order(order = order, n_pv = n_pv)

  return(getXsMNL(X = X, s = s, order = order, n_pv = n_pv))
}




#' Compute choice probabilities of a specific choice set using the MNL model.
#' @param X 3-dimensional design cube/array of size \code{(q, J, S)}.
#' @param beta Parameter vector
#' @param s integer corresponding to a choice set in 1 to S.
#' @param order integer corresponding to a Scheffé model order (1, 2, 3 if no process variables are used; 2 if process variables are used).
#' @param transform_beta boolean parameter. Should the beta vector/matrix be transformed by subtracting the \code{q}-th element?
#' @param n_pv number of process variables.
#'
#' @return J-dimensional numerical vector.
#'
#' @export
mnl_get_Ps = function(X, beta, s, order = 3, transform_beta = T, n_pv = 0){

  mnl_check_order(order = order, n_pv = n_pv)

  Xs = mnl_get_Xs(X, s, order)
  return(as.vector(getPsMNL(X = X, beta = beta, s = s, Xs = Xs, transform_beta = transform_beta, n_pv = n_pv)))
}








#### Functions for FDS plots

#' Perform Monte Carlo simulations for FDS plots.
#'
#' This function creates a dataframe/tibble with Monte Carlo simulations for Fraction of Design Space plots.
#'
#' @param design_array 3-dimensional design array.
#' @param beta Either a vector with the value of the parameter, or a matrix with Halton draws of the prior distribution.
#' @param order integer corresponding to a Scheffé model order (1, 2, 3 if no process variables are used; 2 if process variables are used). This version of the package does not support process variables yet.
#' @param n_points_per_alternative Number of points to use to approximate the prediction variance in each alternative \code{j}, with \code{j} taking values from 1 to \code{J}; and \code{J} is the number of alternatives in each choice set.
#' @param transform_beta boolean parameter. Should the beta vector/matrix be transformed by subtracting the \code{q}-th element?
#' @param verbose Level of verbosity. 0 for no printing, anything bigger than 0 for printing the progress of the progress.
#' @param n_pv Number of process variables. Defaults to 0.
#'
#'
#' @return Tibble of size \code{J*n_points_per_alternative} where \code{J} is the number of alternatives per choice set, i.e., the second element in \code{dim(design_array)}.
#'
#'
#' @examples
#'
#' # Example 1:
#'
#' beta = mnl_create_random_beta(q = 3, order = 3, seed = 2)$beta # Create random beta
#'
#' # Find D-optimal design
#' mnl_design = mnl_mixture_coord_exch(q = 3, J = 2, S = 10, n_random_starts = 1, beta = get_halton_draws(beta, sd = 1, ndraws = 128), max_it = 5)
#'
#' # Get FDS simulations
#' fds_sims = mnl_get_fds_simulations(mnl_design$X, mnl_design$beta, order = 3, n_points_per_alternative = 100, transform_beta = T, verbose = 1)
#'
#' # Plot:
#' fds_sims %>%
#'   ggplot() +
#'   geom_line(aes(fraction, pred_var), size = 0.8)
#'
#'
#'
#'
#'
#'
#'
#'
#' # Example 2:
#'
#' library(dplyr)
#' library(ggplot2)
#'
#' beta = mnl_create_random_beta(q = 3, order = 3, seed = 4)$beta
#' beta_draws = get_halton_draws(beta, sd = 1, ndraws = 64)
#'
#' mnl_I_opt_design = mnl_mixture_coord_exch(
#'   q = 3, J = 2, S = 10,
#'   n_random_starts = 1,
#'   beta = beta_draws,
#'   max_it = 5,
#'   opt_crit = "I")
#'
#' mnl_D_opt_design = mnl_mixture_coord_exch(
#'   q = 3, J = 2, S = 10,
#'   n_random_starts = 1,
#'   beta = beta_draws,
#'   max_it = 5,
#'   opt_crit = "D")
#'
#'
#' fds_sims =  mnl_get_fds_simulations(
#'   design_array = mnl_I_opt_design$X,
#'   beta = mnl_I_opt_design$beta,
#'   order = 3,
#'   n_points_per_alternative = 200,
#'   transform_beta = T,
#'   verbose = 1) %>%
#'   mutate(Design = "I-optimal") %>%
#'   bind_rows(
#'     mnl_get_fds_simulations(
#'       design_array = mnl_D_opt_design$X,
#'       beta = mnl_D_opt_design$beta,
#'       order = 3,
#'       n_points_per_alternative = 200,
#'       transform_beta = T,
#'       verbose = 1) %>%
#'       mutate(Design = "D-optimal")
#'   )
#'
#'
#' # Plot designs
#' fds_sims %>%
#'   ggplot() +
#'   geom_vline(xintercept = 0.5, linetype = "dashed", size = 0.2) +
#'   geom_hline(yintercept = fds_sims %>%
#'                group_by(Design) %>%
#'                summarize(
#'                  med = median(pred_var),
#'                  mean = mean(pred_var)) %>%
#'                dplyr::pull(med),
#'              linetype = "dashed", size = 0.2) +
#'   geom_line(aes(fraction, pred_var, linetype = Design), size = 0.8) +
#'   xlab("Fraction of design space") +
#'   ylab("Prediction variance") +
#'   ggtitle("I-optimal vs D-optimal design") +
#'   theme_bw()
#'
#'
#'
#'
#'
#'
#'
#'
#' Example 3 (with process variables)
#'
#'
#' library(dplyr)
#' library(ggplot2)
#'
#' # A D-optimal utility neutral design with 3 mixture ingredients and 3 process variables
#' d_opt_des = mnl_mixture_coord_exch(q = 3, J = 2, S = 28, n_pv = 3, order = 4, beta = rep(0, 20), transform_beta = F, n_random_starts = 4, max_it = 5, n_cores = 1, opt_crit = "D")
#'
#' # An I-optimal utility neutral design with 3 mixture ingredients and 3 process variables
#' i_opt_des = mnl_mixture_coord_exch(q = 3, J = 2, S = 28, n_pv = 3, order = 4, beta = rep(0, 20), transform_beta = F, n_random_starts = 4, max_it = 5, n_cores = 1, opt_crit = "I")
#'
#'
#' fds_sims =  mnl_get_fds_simulations(
#'   design_array = i_opt_des$X,
#'   beta = i_opt_des$beta,
#'   order = 4,
#'   n_pv = 3,
#'   n_points_per_alternative = 500,
#'   transform_beta = F,
#'   verbose = 1) %>%
#'   mutate(Design = "I-optimal") %>%
#'   bind_rows(
#'     mnl_get_fds_simulations(
#'       design_array = d_opt_des$X,
#'       beta = d_opt_des$beta,
#'       order = 4,
#'       n_pv = 3,
#'       n_points_per_alternative = 500,
#'       transform_beta = F,
#'       verbose = 1) %>%
#'       mutate(Design = "D-optimal")
#'   )
#'
#'
#' # Plot designs
#' fds_sims %>%
#'   ggplot() +
#'   geom_vline(xintercept = 0.5, linetype = "dashed", size = 0.2) +
#'   geom_hline(yintercept = fds_sims %>%
#'                group_by(Design) %>%
#'                summarize(
#'                  med = median(pred_var),
#'                  mean = mean(pred_var)) %>%
#'                dplyr::pull(med),
#'              linetype = "dashed", size = 0.2) +
#'   geom_line(aes(fraction, pred_var, linetype = Design), size = 0.8) +
#'   xlab("Fraction of design space") +
#'   ylab("Prediction variance") +
#'   ggtitle("I-optimal vs D-optimal design") +
#'   theme_bw()
#'
#' @export
mnl_get_fds_simulations = function(design_array, beta, order, n_points_per_alternative = 500, transform_beta = F, verbose = 0, n_pv = 0){


  q = dim(design_array)[1] - n_pv
  J = dim(design_array)[2]
  S = dim(design_array)[3]

  ## Get m (No need to check order because this function checks it).
  m = mnl_get_number_of_parameters(q = q, order = order, n_pv = n_pv)

  #############################################
  ## Check beta dimensions and transform if necessary
  #############################################

  if(is.vector(beta)) {
    if(m != length(beta) & transform_beta) stop("Incompatible size in beta and q: beta must be of length ", m)
    if(m-1 != length(beta) & !transform_beta) stop("Incompatible size in beta and q: beta must be of length ", m-1)
  }

  if(is.matrix(beta)) {
    if(m != ncol(beta) & transform_beta) stop("Incompatible size in beta and q: beta must have ", m,  " columns")
    if(m-1 != ncol(beta) & !transform_beta) stop("Incompatible size in beta and q: beta must have ", m-1,  " columns")
  }


  if(is.vector(beta)) beta_mat = matrix(beta, nrow = 1)
  else beta_mat = beta




  #############################################
  ## Actual function
  #############################################

  pred_var = suppressWarnings(as.data.frame(matrix(rep(NA_real_, J*n_points_per_alternative), ncol = J))) %>%
    purrr::set_names(paste0("V", 1:J))



  progress_old = -1 # For printing progress

  for(k in 1:n_points_per_alternative){

    ## Print progress
    if(verbose > 0){
      progress = round(100*(k-1)/n_points_per_alternative, 0)
      if(progress - progress_old >= 1){
        cat('\r', "Progress: ", progress, "%", sep = "")
        flush.console()
      }
      progress_old = progress
    }

    # Random design
    des_k = suppressWarnings(mnl_create_random_initial_design(q, J, S, seed = k, n_pv = n_pv))

    # Vector to store variances
    vars_1 = rep(NA_real_, J)

    # Iterate over alternatives
    for(j in 1:J){
      f_x = mnl_get_Xs(des_k, 1, order = order, n_pv = n_pv)[j,]
      acc = 0

      # Iterate over prior draws of beta
      for(i in 1:nrow(beta_mat)){
        inf_mat = mnl_get_information_matrix(design_array, beta = beta_mat[i,], order = order, transform_beta = transform_beta, n_pv = n_pv)
        acc = acc + t(f_x) %*% solve(inf_mat, f_x)
      }
      vars_1[j] = acc/nrow(beta_mat)
    }

    pred_var[k,] = vars_1
  }

  pred_var = unlist(pred_var)

  out = dplyr::tibble(
    fraction = (0:length(pred_var))/length(pred_var)
  ) %>%
    dplyr::mutate(pred_var = quantile(pred_var, probs = fraction, type = 1))

  if(verbose > 0) cat("\nFinished\n\n")

  return(out)
}








#' Transform variance covariance matrix.
#'
#' This function transforms the variance covariance matrix of a multivariate normal vector to the identifiable space for a Scheffé model.
#'
#' @param Sigma Original variance covariance matrix of size \code{(m, m)}.
#' @param q Index of the vector that is subtracted to the first q-1 entries.
#'
#' @return Matrix of size \code{(m-1, m-1)}.
#'
#' @export
transform_varcov_matrix = function(Sigma, q){
  m = nrow(Sigma)
  if(m != ncol(Sigma)) stop("Sigma is not a square matrix")
  if(q > m) stop("q can't be greater than the dimension of Sigma")



  Sigma_prime = matrix(rep(0.0, (m-1)^2), m-1)
  for(i in 1:(m-1)){
    for(j in 1:(m-1)){
      if(i <= q-1 & j <= q-1) Sigma_prime[i, j] = Sigma[i, j] - Sigma[i, q] - Sigma[j, q] + Sigma[q, q]
      if(i >= q & j <= q-1) Sigma_prime[i, j] = Sigma[i+1, j] - Sigma[i+1, q]
      if(i <= q-1 & j >= q) Sigma_prime[i, j] = Sigma[i, j+1] - Sigma[q, j+1]
      if(i >= q & j >= q) Sigma_prime[i, j] = Sigma[i + 1, j + 1]
    }
  }
  Sigma_prime

  return(Sigma_prime)

}




# mnl_get_fds_simulations = function(design_array, beta, order, n_points_per_alternative = 500, transform_beta = F, verbose = 1){
#
#   q = dim(design_array)[1]
#   J = dim(design_array)[2]
#   S = dim(design_array)[3]
#
#   pred_var = suppressWarnings(as.data.frame(matrix(rep(NA_real_, J*n_points_per_alternative), ncol = J))) %>%
#     set_names(paste0("V", 1:J))
#   progress_old = -1
#   for(k in 1:n_points_per_alternative){
#
#     if(verbose > 0){
#       progress = round(100*(k-1)/n_points_per_alternative, 0)
#       if(progress - progress_old >= 1){
#         cat('\r', "Progress: ", progress, "%", sep = "")
#         flush.console()
#       }
#       progress_old = progress
#     }
#
#     des_k = mnl_create_random_initial_design(q, J, S, seed = k)
#     vars_1 = rep(NA_real_, J)
#     for(j in 1:J){
#       f_x = mnl_get_Xs(des_k, 1, order = order)[j,]
#       acc = 0
#       for(i in 1:nrow(beta)){
#         inf_mat = mnl_get_information_matrix(design_array, beta = beta[i,], order = order, transform_beta = transform_beta)
#         acc = acc + t(f_x) %*% solve(inf_mat, f_x)
#       }
#       vars_1[j] = acc/nrow(beta)
#     }
#
#     pred_var[k,] = vars_1
#   }
#
#   out = tibble(pred_var = sort(apply(pred_var, 1, sum))) %>%
#     mutate(fraction = 1:nrow(.)/nrow(.))
#
#   if(verbose > 0) cat("\nFinished\n\n")
#
#   return(out)
# }

