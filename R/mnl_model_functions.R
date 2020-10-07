
#' Creation of a random initial design for MNL model
#'
#' \code{create_random_initial_MNL_design} creates a random initial design for MNL model of the specified dimensions
#' @param q integer specifying the number of ingredients
#' @param J integer specifying the number of alternatives within each choice set
#' @param S integer specifying the number of choice sets
#' @param seed integer used for reproducibility
#' @return 3-dimensional array of dimensions (q, J, S)
#' @examples
#' create_random_initial_MNL_design(3, 5, 4, seed = 2020)
#' @export
create_random_initial_MNL_design = function(q, J, S, seed = NULL){
  X = array(rep(NA_real_, q*J*S), dim = c(q, J, S))

  if(!is.null(seed)) set.seed(seed)

  for(j in 1:J){
    for(s in 1:S){
      rands = runif(q)
      # ingredients must sum up to 1
      X[,j, s] = rands/sum(rands)
    }
  }

  return(X)
}



#' Function that computes optimality criterion value for MNL model
#'
#' \code{get_opt_crit_value_MNL} computes optimality criterion value for MNL model
#' @param X 3 dimensional array with dimensions (q, J, S) where:
#'     q is the number of ingredient proportions,
#'     J is the number of alternatives within a choice set,
#'     S is the number of choice sets.
#' @param beta numeric vector containing the parameters. Should be of length (q^3 + 5*q)/6. Must extend to Bayesian case and have the option of it to be a matrix.
#' @param opt_crit optimality criterion: 0 is D-optimality and 1 is I-optimality
#' @return Returns the value of the optimality criterion for this particular design and this beta vector
#' @export
get_opt_crit_value_MNL = function(X, beta, order, opt_crit = 0){

  q = dim(X)[1]

  if(opt_crit == 0){
    # "D-optimality"
    W = matrix(0.0, nrow = 1)
  } else{
    # "I-optimality")
    W = create_moment_matrix_MNL(q, order)
  }

  beta_mat = matrix(beta, byrow = T, nrow = 1)

  return(getOptCritValueMNL(X = X, beta = beta_mat, opt_crit = opt_crit, verbose = 0, W = W))
}



#' Creation of a random parameter vector for MNL model
#'
#' \code{create_random_beta} creates a random parameter vector
#' @param q integer specifying the number of ingredients
#' @return Returns a list in which the first element of the list is
#'     a numerical vector with the parameters and the second element is a matrix with
#'     the indices as in Ruseckaite, et al - Bayesian D-optimal choice designs for mixtures (2017)
#' @examples
#' create_random_beta(3)
#' @export
create_random_beta = function(q, order = 3, seed = NULL){

  stopifnot(order %in% 1:3)
  if(!all.equal(q, floor(q))) stop("q does not seem to an integer.")
  stopifnot(q > 0)

  m1 = q
  m2 = q*(q-1)/2
  m3 = q*(q-1)*(q-2)/6
  if(order == 1){
    m = m1
  } else{
    if(order == 2){
      m = m1 + m2
    } else{
      m = m1 + m2 + m3 # = q + q*(q-1)/2 + q*(q-1)*(q-2)/6
    }
  }

  if(!is.null(seed)) set.seed(seed)

  beta = rep(NA_real_, m)
  beta[1:m1] = rnorm(m1)

  if(order >= 2) beta[(m1+1):(m1+m2)] = rnorm(m2)

  if(order >= 3) beta[(m1+m2+1):m] = rnorm(m3)

  # Return a list for backwards compatibility, but eventually I'll have to change it to return only a vector
  return(list(beta = beta))
}


#' Coordinate exchange algorithm for a Multinomial Logit Scheffé model.
#'
#' \code{mixture_coord_ex_mnl} Performs the coordinate exchange algorithm for a Multinomial Logit Scheffé model.
#' @param q number of ingredient proportions.
#' @param J number of alternatives within a choice set.
#' @param S number of choice sets.
#' @param n_random_starts number or random starts. Defaults to 100.
#' @param X If an initial design is to be supplied, thenit must be a 3 dimensional array with dimensions (q, J, S), with q, J, and S are defined above.
#' @param beta Prior parameters. For a locally optimal design, it should be a numeric vector of length m = (q^3 + 5*q)/6. For a pseudo-Bayesian design, it must be a matrix with prior simulations of size (nxm) where m is previously defined and m is the number of prior draws, i.e., there is a prior draw per row.
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
#'
#' @return list with 7 elements. See below for details.
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
#' Return list has 7 elements:
#' \itemize{
#'     \item X_orig: The original design. Array with dimensions (q, J, S).
#'     \item X: The optimized design. Array with dimensions (q, J, S).
#'     \item beta: The original beta vector or matrix.
#'     \item opt_crit_value_orig: efficiency of the original design.
#'     \item opt_crit_value: efficiency of the optimized design.
#'     \item n_iter: Number of iterations performed.
#'     \item opt_crit: The optimality criterion used.
#'  }
#'
#' @export
mixture_coord_ex_mnl = function(
  q = NULL,
  J = NULL,
  S = NULL,
  n_random_starts = 100,
  X = NULL,
  beta,
  order = 3,
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

  t1 = Sys.time()
  cat("Starts at", substr(as.character(t1), 12, 19), "\n")


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
      des = create_random_initial_MNL_design(q, J, S, seed = x)
      return(des)
    })

  } else{
    # Some input checks
    dim_X = dim(X)

    if(length(dim_X) != 3) stop("X must be a 3 dimensional array.")
    if(!(is.vector(beta) | is.matrix(beta))) stop("beta is not a vector or a matrix. It must be a numerical or integer vector or matrix.")

    q = dim_X[1]

    designs = list(X)
  }


  # m = (q*q*q + 5*q)/6
  if(order == 1){
    m = q
  } else{
    if(order == 2){
      m = q*(q-1)/2 + q
    } else{
      m = (q^3+ 5*q)/6 # = q + q*(q-1)/2 + q*(q-1)*(q-2)/6
    }
  }

  if(is.vector(beta)) if(m != length(beta)) stop("Incompatible size in beta and q: beta must be of length (q^3 + 5*q)/6")

  if(is.matrix(beta)) if(m != ncol(beta)) stop("Incompatible size in beta and q: beta must have (q^3 + 5*q)/6 columns")


  if(is.vector(beta)) beta_mat = matrix(beta, nrow = 1)
  else beta_mat = beta



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
    W = create_moment_matrix_MNL(q, order = order)
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

    out = try(
      mixtureCoordinateExchangeMNL(
        X_orig = X,
        beta = beta_mat,
        order = order,
        max_it = max_it,
        verbose = verbose,
        opt_crit = opt_crit,
        W = W,
        opt_method = opt_method,
        lower = 0,
        upper = 1,
        tol = tol,
        n_cox_points = n_cox_points
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

  # Return the result with the best optimality criterion
  X_result = results[[which.min(optimality_values)]]


  out_list = list(
    X_orig = X_result$X_orig,
    X = X_result$X,
    beta = beta,
    opt_crit_value_orig = X_result$opt_crit_value_orig,
    opt_crit_value = X_result$opt_crit_value,
    n_iter = X_result$n_iter,
    opt_crit = ifelse(opt_crit == 0, "D-optimality", "I-optimality")
  )

  if(plot_designs) {
    if(q == 3 | q == 4) mnl_plot_result(out_list)
    else warning("Could not plot results because q is not 3 or 4.")
  }

  t2 = Sys.time()

  cat("Ends at", substr(as.character(t2), 12, 19), "\n")
  cat("Time:", round(as.numeric(difftime(t2, t1, units = "secs")), 1), "seconds.\n")

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


#' TODO: write doc
#' @export
mnl_plot_result_3d = function(
  res_alg,
  design = "final",
  phi = 40, theta = 140, cex_points = 0.5, cex_names = 0.5, ...){
  # res_alg: output of a call to mixture_coord_ex_mnl() function.
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
  q = dim_X[1]
  S = dim_X[3]

  if(q != 4) stop("Design must be of 4 ingredients.")

  # Convert 3 dimensional arrays into matrices by vertically binding them
  X_mat = t(res_alg[[i]][,,1])
  for(s in 2:S){
    X_mat = rbind(X_mat, t(res_alg[[i]][,,s]))
  }

  plot_simplex_4d(X_mat,  phi = phi, theta = theta, cex_points = cex_points, cex_names = cex_names, ...)

}





#' TODO: write doc
#' @export
mnl_plot_result_2d = function(res_alg){
  # res_alg: output of a call to mixture_coord_ex_mnl() function.
  # It must be a design of 3 ingredients.

  dim_X = dim(res_alg$X_orig)
  q = dim_X[1]
  S = dim_X[3]

  if(q != 3) stop("Design must be of 3 ingredients.")

  # Convert 3 dimensional arrays into matrices by vertically binding them
  X_orig_mat = t(res_alg$X_orig[,,1])
  for(s in 2:S){
    X_orig_mat = rbind(X_orig_mat, t(res_alg$X_orig[,,s]))
  }

  X_final_mat = t(res_alg$X[,,1])
  for(s in 2:S){
    X_final_mat = rbind(X_final_mat, t(res_alg$X[,,s]))
  }

  # PLot matrices
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







#' TODO: write doc
#' @export
mnl_plot_result = function(res_alg, ...){
  dim_X = dim(res_alg$X_orig)
  q = dim_X[1]

  if(!(q == 3 | q == 4)) stop("Design must be of 3 or 4 ingredients.")

  if(q == 3){
    mnl_plot_result_2d(res_alg = res_alg)
  } else{
    mnl_plot_result_3d(res_alg = res_alg, ...)
  }

}









#' TODO: write doc
#' @export
create_moment_matrix_MNL = function(q, order = 3){

  stopifnot(order %in% 1:3)

  if(order == 1){
    m = q
  } else{
    if(order == 2){
      m = q*(q-1)/2 + q
    } else{
      m = (q^3+ 5*q)/6 # = q + q*(q-1)/2 + q*(q-1)*(q-2)/6
    }
  }

  f = lapply(1:(m-1), function(x) rep(0, q))

  counter = 0
  # Fill indicators of first part of the model expansion
  for(i in 1:(q-1)){
    counter = counter + 1
    f[[counter]][i] = 1
  }

  # Fill indicators of second part of the model expansion
  if(order >= 2){
    for(i in 1:(q-1)){
      for(j in (i+1):q){
        counter = counter + 1
        f[[counter]][i] = 1
        f[[counter]][j] = 1
      }
    }
  }


  # Fill indicators of third part of the model expansion
  if(order >= 3){
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
  }


  W = matrix(rep(NA_real_, (m-1)^2), ncol = m-1)

  for(i in 1:(m-1)){
    for(j in 1:(m-1)){

      aux_ij = f[[i]] + f[[j]]
      num_ij = prod(factorial(aux_ij))
      denom_ij = factorial(2 + sum(aux_ij))
      W[i,j] = num_ij/denom_ij
    }
  }

  return(W)

}





