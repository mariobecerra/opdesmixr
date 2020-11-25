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


#### functions for building vignettes

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
