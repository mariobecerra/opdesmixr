create_2d_simplex_grid = function(n_out){
  # Returns a grid of points inside the 2-dimensional simplex.
  # The output is a dataframe of approximately n_out rows.

  n_points = floor(sqrt(2*n_out))

  simplex_grid = expand.grid(x1 = seq(0, 1, length.out = n_points),
                             x2 = seq(0, 1, length.out = n_points)) %>%
    filter(x1 + x2 <= 1) %>%
    as_tibble() %>%
    mutate(x3 = 1-x1-x2)

  return(simplex_grid)
}




scheffe_order3_q3_utilities = function(beta, n_out = 1000){

  out = create_2d_simplex_grid(n_out) %>%
    mutate(utility = beta[1]*x1 + beta[2]*x2 + beta[3]*x1*x2 + beta[4]*x1*x3 + beta[5]*x2*x3 + beta[6]*x1*x2*x3)

  return(out)
}





plot_choice_set_utility = function(
  design_object,
  utility_data,
  design_point_size = 2,
  utility_point_size = 0.3,
  utility_point_shape = "square",
  low_color_gradient = "light blue",
  high_color_gradient = "dark blue",
  legend.position = "bottom",
  legend.box = "vertical"){

  if(class(design_object) == "list"){
    # If object is the result of the algorithm
    dim_X = dim(design_object$X)
    q = dim_X[1]
    S = dim_X[3]

    if(q != 3) stop("Design must be of 3 ingredients.")

    X_final_tbl = mnl_design_array_to_dataframe(design_object$X)

  } else{
    # If it's just an array
    if(inherits(design_object, "array")){

      dim_X = dim(design_object)
      q = dim_X[1]
      S = dim_X[3]

      if(q != 3) stop("Design must be of 3 ingredients.")


      X_final_tbl = mnl_design_array_to_dataframe(design_object)

    } else{
      if(inherits(design_object, "data.frame")){
        q = ncol(design_object) - 1

        if(q != 3) stop("Design must be of 3 ingredients.")

        X_final_tbl = design_object %>%
          set_names(c("c1", "c2", "c3", "choice_set")) %>%
          mutate(choice_set = as.character(choice_set))

      } else{
        stop("Unknown type of design")
      }

    }
  }

  out_plot = ggtern::ggtern(data = utility_data) +
    ggplot2::geom_point(
      aes(x = x1, y = x2, z = x3, color = utility_avg),
      size = utility_point_size,
      shape = utility_point_shape) +
    ggplot2::scale_color_gradient(
      low = low_color_gradient, high = high_color_gradient) +
    ggplot2::geom_point(
      data = X_final_tbl,
      ggtern::aes(c1, c2, c3, shape = choice_set),
      # https://stackoverflow.com/questions/26223857/more-than-six-shapes-in-ggplot
      color = "black", size = design_point_size, stroke = 1.5, inherit.aes = F) +
    ggtern::theme_nomask() +
    scale_shape_manual(values = 1:length(unique(X_final_tbl$choice_set))) +
    theme(legend.position = legend.position, legend.box = legend.box)

  return(out_plot)

}










#### Functions for FDS plots


get_pred_var1 = function(design_array, beta, order, n_points_per_alternative = 500, transform_beta = F){

  q = dim(design_array)[1]
  J = dim(design_array)[2]
  S = dim(design_array)[3]

  pred_var = suppressWarnings(as.data.frame(matrix(rep(NA_real_, J*n_points_per_alternative), ncol = J))) %>%
    set_names(paste0("V", 1:J))
  for(k in 1:n_points_per_alternative){
    des_k = mnl_create_random_initial_design(q, J, S, seed = k)
    vars_1 = rep(NA_real_, J)
    for(j in 1:J){
      f_x = mnl_get_Xs(des_k, 1, order = order)[j,]
      acc = 0
      for(i in 1:nrow(beta)){
        inf_mat = mnl_get_information_matrix(design_array, beta = beta[i,], order = order, transform_beta = transform_beta)
        acc = acc + t(f_x) %*% solve(inf_mat, f_x)
      }
      vars_1[j] = acc/nrow(beta)
    }

    pred_var[k,] = vars_1
  }


  out = tibble(pred_var = sort(unlist(pred_var))) %>%
    mutate(fraction = 1:nrow(.)/nrow(.))

  return(out)
}




get_pred_var = function(design_array, beta, order, n_points_per_alternative = 500, transform_beta = F, verbose = 1){

  q = dim(design_array)[1]
  J = dim(design_array)[2]
  S = dim(design_array)[3]

  pred_var = suppressWarnings(as.data.frame(matrix(rep(NA_real_, J*n_points_per_alternative), ncol = J))) %>%
    set_names(paste0("V", 1:J))
  progress_old = -1
  for(k in 1:n_points_per_alternative){

    if(verbose > 0){
      progress = round(100*(k-1)/n_points_per_alternative, 0)
      if(progress - progress_old >= 1){
        cat('\r', "Progress: ", progress, "%", sep = "")
        flush.console()
      }
      progress_old = progress
    }

    des_k = mnl_create_random_initial_design(q, J, S, seed = k)
    vars_1 = rep(NA_real_, J)
    for(j in 1:J){
      f_x = mnl_get_Xs(des_k, 1, order = order)[j,]
      acc = 0
      for(i in 1:nrow(beta)){
        inf_mat = mnl_get_information_matrix(design_array, beta = beta[i,], order = order, transform_beta = transform_beta)
        acc = acc + t(f_x) %*% solve(inf_mat, f_x)
      }
      vars_1[j] = acc/nrow(beta)
    }

    pred_var[k,] = vars_1
  }

  out = tibble(pred_var = sort(apply(pred_var, 1, sum))) %>%
    mutate(fraction = 1:nrow(.)/nrow(.))

  if(verbose > 0) cat("\nFinished\n\n")

  return(out)
}





get_file_path = function(filepath){
  # This function is useful to use when developing the package and working with data in the inst/ folder.
  # If the package hasn't been installed, then system.file() won't work, then it goes to the folder.
  # It assumes that the working directory is the project directory.
  path_system_file = system.file(filepath, package = "opdesmixr", mustWork = F)
  if(path_system_file == ""){
    out_path = here::here(filepath)
  } else{
    out_path = path_system_file
  }
  return(out_path)
}



get_file_path_inst = function(filepath){
  # This function is useful to use when developing the package and working with data in the inst/ folder.
  # If the package hasn't been installed, then system.file() won't work, then it goes to the folder.
  # It assumes that the working directory is the project directory.
  path_system_file = system.file(filepath, package = "opdesmixr", mustWork = F)
  if(path_system_file == ""){
    out_path = here::here("inst/", filepath)
  } else{
    out_path = path_system_file
  }
  return(out_path)
}
