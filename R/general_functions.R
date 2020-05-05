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
      ggtern(aes(c1, c2, c3)) +
      geom_path(linetype = "dashed") +
      theme_minimal() +
      theme_nomask() +
      geom_point(data = tibble(c1 = x_in[1], c2 = x_in[2], c3 = x_in[3]))
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
      ggtern(aes(c1, c2, c3)) +
      geom_path(linetype = "dashed", aes(group = comp)) +
      theme_minimal() +
      theme_nomask() +
      geom_point(data = tibble(c1 = x_in[1], c2 = x_in[2], c3 = x_in[3]))
  }

  return(out)
}


# Main function

#' TODO: write doc
#' @export
mixture_coord_ex = function(
  X,
  order = NULL,
  beta = NULL,
  model = "Gaussian",
  n_cox_points = 100,
  max_it = 50,
  plot_designs = F,
  opt_crit = 0, # 0 is D-optimality and 1 is I-optimality
  verbose = 1){

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

    cat('Creating design for "', model, '" model of order ', order, '.\n', sep = "")

    if(!is.null(beta)) warning("beta was supplied but was ignored.")
    out = mixture_coord_ex_gaussian(
      X = X,
      order = order,
      n_cox_points = n_cox_points,
      max_it = max_it,
      plot_designs = plot_designs,
      verbose = verbose,
      opt_crit = opt_crit)
  }

  # Call MNL function
  if(model == "MNL"){

    cat('Creating design for "', model, '" model.\n', sep = "")

    if(!is.null(order)) warning("order was supplied but was ignored.")

    out = mixture_coord_ex_mnl(
      X = X,
      beta = beta,
      n_cox_points = n_cox_points,
      max_it = max_it,
      plot_designs = plot_designs,
      verbose = verbose,
      opt_crit = opt_crit)
  }

  return(out)

}
