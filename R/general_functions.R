#' TODO: write doc
#' @export
compute_cox_direction = function(x, comp, n_points = 11){
  # Call C++ function
  cox_direction = computeCoxDirection(x, comp, n_points, verbose = 0)

  return(cox_direction)
}



#' TODO: write docs
#' @export
plot_cox_direction = function(x_in, comp = NULL, n_points = 3){
  # x_in: vector of length 3 that sums up to 1
  # comp: which component should be plotted. Can be a scalar or a vector

  if(length(x_in) != 3) stop("x_in must be of length 3")
  if(sum(x_in) != 1) stop("x_in must sum up to 1")

  # If component is null
  if(!is.null(comp) & length(comp) == 1){
    out = compute_cox_direction(x_in, comp, n_points) %>%
      dplyr::as_tibble() %>%
      purrr::set_names(c("c1", "c2", "c3")) %>%
      ggtern::ggtern(ggtern::aes(c1, c2, c3)) +
      ggplot2::geom_path(linetype = "dashed") +
      ggplot2::geom_point(data = tibble(c1 = x_in[1], c2 = x_in[2], c3 = x_in[3])) +
      ggtern::theme_minimal() +
      ggtern::theme_nomask()
  } else{
    if(is.null(comp)) comp = 1:length(x_in)

    cox_dirs = lapply(1:length(comp), function(i){
      compute_cox_direction(x_in, i, 3) %>%
        dplyr::as_tibble() %>%
        purrr::set_names(c("c1", "c2", "c3")) %>%
        dplyr::mutate(comp = i)
    }) %>%
      bind_rows()

    out = cox_dirs %>%
      ggtern::ggtern(ggtern::aes(c1, c2, c3)) +
      ggplot2::geom_path(linetype = "dashed", aes(group = comp)) +
      ggplot2::geom_point(data = tibble(c1 = x_in[1], c2 = x_in[2], c3 = x_in[3])) +
      ggtern::theme_minimal() +
      ggtern::theme_nomask()
  }

  return(out)
}


# Main function

#' TODO: write doc
#' @export
mixture_coord_ex = function(
  n_random_starts = 100,
  X = NULL,
  order = NULL,
  n_runs = NULL,
  q = NULL,
  J = NULL,
  S = NULL,
  beta = NULL,
  model = "Gaussian",
  opt_method = "B",
  max_it = 10,
  tol = 0.0001,
  n_cox_points = NULL,
  plot_designs = F,
  verbose = 1,
  opt_crit = 0,
  seed = NULL,
  n_cores = 1){

  available_models = c("Gaussian", "MNL")

  # For martial match
  pmatch_vec = sapply(available_models,
                      function(x) pmatch(model, x))

  model = names(which(!is.na(pmatch_vec)))

  if(all(is.null(pmatch_vec))){
    # if no model is recognized
    stop('Model unknown. Must be one of the following:\n\t',
         paste(available_models, collapse = ", "))
  }


  # Call Gaussian function
  if(model == "Gaussian"){

    if(is.null(order)) stop("Must supply order.")

    if(verbose >= 1) cat('Creating design for "', model, '" model of order ', order, '.\n', sep = "")

    if(!is.null(beta)) warning("beta was supplied but was ignored.")
    if(!is.null(J)) warning("J was supplied but was ignored.")
    if(!is.null(S)) warning("S was supplied but was ignored.")
    if(is.null(X) & (is.null(n_runs) | is.null(q))) stop("Must supply X, or both n_runs and q.")

    out = gaussian_mixture_coord_exch(
      n_runs = n_runs,
      q = q,
      n_random_starts = n_random_starts,
      X = X,
      order = order,
      opt_method = opt_method,
      max_it = max_it,
      tol = tol,
      n_cox_points = n_cox_points,
      plot_designs = plot_designs,
      verbose = verbose,
      opt_crit = opt_crit,
      seed = seed,
      n_cores = n_cores)
  }

  # Call MNL function
  if(model == "MNL"){

    if(verbose >= 1) cat('Creating design for "', model, '" model.\n', sep = "")

    if(is.null(beta)) stop("Must supply beta for MNL model.")
    if(is.null(X) & (is.null(q) | is.null(J) | is.null(S))) stop("Must supply X, or all of J, Q, and q.")

    if(!is.null(order)) warning("order was supplied but was ignored. MNL only works with cubic models")
    if(!is.null(n_runs)) warning("n_runs was supplied but was ignored")


    out = mnl_mixture_coord_exch(
      q = q,
      J = J,
      S = S,
      n_random_starts = n_random_starts,
      X = X,
      beta = beta,
      opt_method = opt_method,
      max_it = max_it,
      tol = tol,
      n_cox_points = n_cox_points,
      plot_designs = plot_designs,
      verbose = verbose,
      opt_crit = opt_crit,
      seed = seed,
      n_cores = n_cores)
  }

  return(out)

}






#' TODO: write doc
#' Returns a matrix of dimension m x ndraws where m is the length of the beta vector.
#' Example on how to use:
#'     beta_means = create_random_beta(q)
#'     beta_prior_draws = get_halton_draws(beta_means$beta)
#' @export
get_halton_draws = function(beta, sd = 1, ndraws = 120){
  draws_unif_beta = randtoolbox::halton(ndraws, dim = length(beta))

  beta_prior_draws = matrix(rep(NA_real_, length(draws_unif_beta)), ncol = ncol(draws_unif_beta))
  for(i in 1:ncol(beta_prior_draws)){
    beta_prior_draws[, i] = qnorm(draws_unif_beta[, i], mean = beta[i], sd = sd)
  }

  return(beta_prior_draws)
}












#' TODO: write doc
#' Returns a matrix of dimension m x ndraws where m is the length of the beta vector and sigma is the correlation matrix between the parameters.
#' Matrix sigma must be if size mxm.
#' Example on how to use:
#' beta = c(1.36, 1.57, 2.47, -0.43, 0.50, 1.09)
#'
#' sigma = matrix(
#'   c(6.14, 5.00, 2.74, -0.43, -2.81, -3.33,
#'     5.00, 6.76, 4.47, -1.79, -6.13, -3.51,
#'     2.74, 4.47, 3.45, -1.38, -4.71, -2.17,
#'     -0.43, -1.79, -1.38, 1.18, 2.39, 0.71,
#'     -2.81, -6.13, -4.71, 2.39, 7.43, 2.71,
#'     -3.33, -3.51, -2.17, 0.71, 2.71, 2.49),
#'   ncol = 6,
#'   byrow = T)
#' @export
get_correlated_halton_draws = function(beta, sigma, n_draws = 128){

  uncorrelated_draws = get_halton_draws(beta, sd = 1, ndraws = n_draws)
  uncorrelated_draws_mean = apply(uncorrelated_draws, 2, mean)

  chol_sigma = chol(sigma)

  correlated_draws = uncorrelated_draws %*% chol_sigma
  correlated_draws2 = correlated_draws - matrix(rep(apply(correlated_draws, 2, mean), n_draws), nrow = n_draws, byrow = T) + matrix(rep(uncorrelated_draws_mean, n_draws), nrow = n_draws, byrow = T)
  return(correlated_draws2)
}













#' Explicitly draw plot (\code{ggtern} and \code{ggplot2} compatible)
#'
#' Makes sure both \code{\link[ggtern:::print.ggplot]{ggtern}} and
#' \code{\link[ggplot2:::print.ggplot]{ggplot2}} objects
#' display correctly.
#'
#' @param x plot to display
#' @param newpage draw new (empty) page first?
#' @param vp viewport to draw plot in
#' @param ... other arguments not used by this method
#'
#' @return Invisibly returns the result of \code{ggplot_build()}, which is a
#' list with components that contain the plot itself, the data, information
#' about the scales, panels etc.
#' Created in https://github.com/skgallagher/ following discussion in https://bitbucket.org/nicholasehamilton/ggtern/issues/13/ggtern-breaks-ggplots-plot
#' @export
print.ggplot <- function(x, newpage = is.null(vp), vp = NULL, ...){
  if(inherits(x$coordinates, "CoordTern")){
    ggtern:::print.ggplot(x, newpage = newpage, vp = vp, ...)

  } else {
    ggplot2:::print.ggplot(x, newpage = newpage, vp = vp, ...)
  }
}





