library(opdesmixr)
library(tidyverse)
library(xtable)
library(here)

# get_file_path = function(filepath){
#   # This function is useful to use when developing the package and working with data in the inst/ folder.
#   # If the package hasn't been installed, then system.file() won't work, then it goes to the folder.
#   # It assumes that the working directory is the project directory.
#   path_system_file = system.file(filepath, package = "opdesmixr", mustWork = F)
#   if(path_system_file == ""){
#     out_path = here::here("inst/", filepath)
#   } else{
#     out_path = path_system_file
#   }
#   return(out_path)
# }
#
# source(get_file_path("article_results/R_scripts/functions.R"))

# Assuming this script is run before installing package.
designs_folder = here("inst/misc_output/cocktail_cornell_designs/")

out_folder = here("inst/misc_output/cocktail_cornell_plots/")
dir.create(out_folder, showWarnings = F)


# Folder with the designs by Ruseckaite
ruseckaite_designs_folder = opdesmixr:::get_file_path_inst("exdata/")




##########################################################################################
#### Functions
##########################################################################################


mnl_get_choice_probabilities = function(X, beta, order, transform_beta = T){
  # Returns a tibble with the choice probabilities of each choice set.
  # The output is a tibble in which each row is a choice set and the first J columns are the probability of choosing that option.
  dim_X = dim(X)
  S = dim_X[3]

  probs_df = as.data.frame(matrix(rep(NA_real_, dim_X[2]*S), nrow = S)) %>%
    set_names(paste0("p_", 1:dim_X[2]))

  for(s in 1:S){
    probs_df[s,] = mnl_get_Ps(X, beta, s, 3, T)
  }

  probs_df %>%
    as_tibble() %>%
    mutate(choice_set = 1:nrow(.)) %>%
    return()

}




get_product_of_choice_probabilities = function(design_array, beta, order = 3, transform_beta = T){
  # Returns a tibble with the choice probabilities of each choice set. and its product.
  # The output is a tibble in which each row is a choice set and the first J columns are the probability of choosing that option.
  # The (J+1)-th column denotes the choice set, and the last column contains the product of the choice probabilities in each choice set.

  choice_probs = mnl_get_choice_probabilities(design_array, beta, order, transform_beta)

  prods = Reduce(`*`, select(choice_probs, all_of(1:(ncol(choice_probs)-1))))

  choice_probs %>%
    mutate(prods = prods) %>%
    return()

}




get_distances_within_choice_set = function(design_array){
  # Computes the Euclidean distances between the two alternatives within a choice set for a design.

  dim_array = dim(design_array)
  J = dim_array[2]
  if(J != 2) stop("This function only works for designs with two alternatives within each choice set.")

  S = dim_array[3]
  distances = tibble(dist = rep(NA_real_, S))

  for(s in 1:S){
    square_diff = (design_array[,1,s] - design_array[,2,s])^2
    distances[s,] = sqrt(sum(square_diff))
  }

  distances %>%
    mutate(choice_set = 1:nrow(.)) %>%
    return()

}





# Plotting function
plot_choice_set_utility_discrete_palette = function(
  design_object, utility_data, design_point_size = 2,
  utility_point_size = 0.3, utility_point_shape = "square",
  legend.position = "bottom", legend.box = "vertical",
  color_palette = NULL) {

  if (class(design_object) == "list") {
    dim_X = dim(design_object$X)
    q = dim_X[1]
    S = dim_X[3]
    if (q != 3)
      stop("Design must be of 3 ingredients.")
    X_final_tbl = mnl_design_array_to_dataframe(design_object$X)
  }
  else {
    if (inherits(design_object, "array")) {
      dim_X = dim(design_object)
      q = dim_X[1]
      S = dim_X[3]
      if (q != 3)
        stop("Design must be of 3 ingredients.")
      X_final_tbl = mnl_design_array_to_dataframe(design_object)
    }
    else {
      if (inherits(design_object, "data.frame")) {
        q = ncol(design_object) - 1
        if (q != 3)
          stop("Design must be of 3 ingredients.")
        X_final_tbl = design_object %>% set_names(c("c1",
                                                    "c2", "c3", "choice_set")) %>% mutate(choice_set = as.character(choice_set))
      }
      else {
        stop("Unknown type of design")
      }
    }
  }

  if(is.null(color_palette)){
    # Create palette that goes from yellow to red
    green_levels = 1 - seq(0, 1, length.out = length(unique(utilities_cocktail_bayesian_plot$utility_int)))

    color_palette = tibble(g = green_levels) %>%
      mutate(r = 1, b = 0) %>%
      mutate(color = rgb(r, g, b)) %>%
      pull(color)
  }


  out_plot = ggtern::ggtern(data = utility_data) +
    ggplot2::geom_point(
      aes(x = x1,
          y = x2, z = x3, color = factor(utility_int)),
      size = utility_point_size,
      shape = utility_point_shape) +
    ggplot2::scale_color_manual(values = color_palette) +
    ggplot2::geom_point(data = X_final_tbl, ggtern::aes(c1, c2, c3, shape = choice_set), color = "black", size = design_point_size, stroke = 1.5, inherit.aes = F) +
    ggtern::theme_nomask() +
    scale_shape_manual(values = 1:length(unique(X_final_tbl$choice_set))) +
    theme(legend.position = legend.position, legend.box = legend.box)
  return(out_plot)
}







min_dist_i = function(design_tibble, i, use_unique_points = F){
  # design_tibble must be a tibble where the first q columns are ingredient proportions and the q-th column is a choice set indicator.
  # If use_unique_points is TRUE, the function will only use unique points

  q = ncol(design_tibble) - 1

  xi = as.numeric(design_tibble[i, 1:q])

  if(use_unique_points){
    design_tibble = design_tibble %>%
      select(all_of(1:q)) %>%
      distinct()
  }

  min_dist = Inf
  for(j in setdiff(1:nrow(design_tibble), i)){

    xj = as.numeric(design_tibble[j, 1:q])

    dist_ij = sqrt(sum((xi - xj)^2))
    min_dist = min(dist_ij, min_dist)

  }

  return(min_dist)

}






min_dist = function(design_tibble, use_unique_points = F){
  # design_tibble must be a tibble where the first q columns are ingredient proportions and the q-th column is a choice set indicator.
  # If use_unique_points is TRUE, the function will only use unique points

  q = ncol(design_tibble) - 1

  if(use_unique_points){
    design_tibble = design_tibble %>%
      select(all_of(1:q)) %>%
      distinct()
  }

  min_dist = Inf
  for(i in 1:(nrow(design_tibble)-1)){

    min_dist = min(min_dist, min_dist_i(design_tibble = design_tibble, i = i, use_unique_points = use_unique_points))

  }

  return(min_dist)

}






mean_dist = function(design_tibble, use_unique_points = F){
  # design_tibble must be a tibble where the first q columns are ingredient proportions and the q-th column is a choice set indicator.
  # If use_unique_points is TRUE, the function will only use unique points

  q = ncol(design_tibble) - 1

  if(use_unique_points){
    design_tibble = design_tibble %>%
      select(all_of(1:q)) %>%
      distinct()
  }

  N = nrow(design_tibble)
  mean_dist_acc = 0.0
  for(i in 1:N){

    mean_dist_acc = mean_dist_acc + min_dist_i(design_tibble = design_tibble, i = i, use_unique_points = use_unique_points)

  }
  return(mean_dist_acc/N)

}











coverage = function(design_tibble, use_unique_points = F){

  if(use_unique_points){
    design_tibble = design_tibble %>%
      select(all_of(1:q)) %>%
      distinct()
  }

  N = nrow(design_tibble)
  min_distances = rep(NA_real_, N)
  for(i in seq_along(min_distances)){
    min_distances[i] = min_dist_i(design_tibble = design_tibble, i, use_unique_points = use_unique_points)
  }

  gamma_bar = mean(min_distances)

  return(sqrt(sum((min_distances - gamma_bar)^2)/N)/gamma_bar)
}






std_dev_of_min_dist = function(design_tibble, use_unique_points = F){

  if(use_unique_points){
    design_tibble = design_tibble %>%
      select(all_of(1:q)) %>%
      distinct()
  }

  N = nrow(design_tibble)
  min_distances = rep(NA_real_, N)
  for(i in seq_along(min_distances)){
    min_distances[i] = min_dist_i(design_tibble = design_tibble, i, use_unique_points = use_unique_points)
  }

  gamma_bar = mean(min_distances)

  return(sqrt(sum((min_distances - gamma_bar)^2)/N))

}





transform_from_pseudocomp_to_comp = function(x, l, sum_l){
  a = l + (1 - sum_l)*x
  return(a)
}


transform_tibble_from_pseudocomp_to_comp = function(design_tibble, lower_bounds, var_indices){
  # Returns a tibble with p + q columns and n rows with the transformation from pseudo-components to the original components in the mixture. Here, q is the number of ingredietns proportions, p is the number of original columns in design_tibble, and n is the number of rows in design_tibble.
  sum_lower_bounds = sum(lower_bounds)
  q = length(lower_bounds)

  if(q != length(var_indices)) stop("Incompatible sizes of lower_bounds and var_indices")

  transformed_df = as_tibble(matrix(rep(NA_real_, nrow(design_tibble)*q), ncol = q)) %>%
    set_names(paste0(names(design_tibble)[var_indices], "_original"))

  for(i in seq_along(var_indices)){
    transformed_df[[i]] = transform_from_pseudocomp_to_comp(design_tibble[[var_indices[i]]], lower_bounds[i], sum_lower_bounds)
  }

  return(bind_cols(design_tibble, transformed_df))
}



##########################################################################################
#### Load designs for cocktail experiment
##########################################################################################

beta0 = c(1.36, 1.57, 2.47, -0.43, 0.50, 1.09)

# x1, mango juice, 0.3
# x2, blackcurrant syrup, 0.15
# x3, lemon juice, 0.1
lower_bounds_cocktail = c(0.3, 0.15, 0.1)
sum_lower_bounds = sum(lower_bounds_cocktail)

cocktail_d_opt_filename = paste0(designs_folder, "cocktail_d_optimal.rds")
cocktail_i_opt_filename = paste0(designs_folder, "cocktail_i_optimal.rds")

cocktail_D_opt = readRDS(cocktail_d_opt_filename)
cocktail_I_opt = readRDS(cocktail_i_opt_filename)

ruseckaite_cocktail_designs = read_csv(paste0(ruseckaite_designs_folder, "/ruseckaite_cocktail_designs.csv"))

ruseckaite_cocktail_beta0_design = mnl_design_dataframe_to_array(
  ruseckaite_cocktail_designs %>%
    filter(prior == "beta_0 and sigma_0") %>%
    select(-prior)
)



##########################################################################################
#### Load designs for Cornell's experiment
##########################################################################################

beta_2 = c(0.86, 0.21, 0, 3.07, 2.34, 3.24, -20.59)
beta_2_prime = c(0.86, 0.21, 3.07, 2.34, 3.24, -20.59)
kappas = c(0.5, 5, 10, 30)

n_random_starts_2 = 128
max_it_bayes = 20

cornell_designs_basefilename_transf = paste0(designs_folder, "cornell_experiment_transformed_betas_rs", n_random_starts_2, "_maxit", max_it_bayes)

cornell_designs_basefilename_untransf = paste0(designs_folder, "cornell_experiment_untransformed_betas_rs", n_random_starts_2, "_maxit", max_it_bayes)

ruseckaite_cornell_designs = read_csv(paste0(ruseckaite_designs_folder, "/ruseckaite_cornell_designs.csv"))




# ## Read designs with transformed betas and save them in a list
# cornell_designs_transf = lapply(kappas, function(k){
#
#   cat("kappa =", k, "\n")
#
#   d_opt_filename = paste0(cornell_designs_basefilename_transf, "_kappa", k, "_Dopt.rds")
#   i_opt_filename = paste0(cornell_designs_basefilename_transf, "_kappa", k, "_Iopt.rds")
#
#
#
#   if(file.exists(d_opt_filename)){
#     cat("\tD optimal file exists. Loading.\n")
#     cornell_beta_2_bayesian_d_opt = readRDS(d_opt_filename)
#   }else{
#     stop("\tD optimal file does not exist.\n")
#
#   }
#
#   if(file.exists(i_opt_filename)){
#     cat("\tI optimal file exists. Loading.\n")
#     cornell_beta_2_bayesian_i_opt = readRDS(i_opt_filename)
#   }else{
#     stop("\tI optimal file does not exist.\n")
#   }
#
#   return(list(
#     d_opt = cornell_beta_2_bayesian_d_opt,
#     i_opt = cornell_beta_2_bayesian_i_opt,
#     kappa = k
#   ))
#
# })


## Read designs with untransformed betas and save them in a list
cornell_designs_untransf = lapply(kappas, function(k){

  cat("kappa =", k, "\n")

  d_opt_filename = paste0(cornell_designs_basefilename_untransf, "_kappa", k, "_Dopt.rds")
  i_opt_filename = paste0(cornell_designs_basefilename_untransf, "_kappa", k, "_Iopt.rds")



  if(file.exists(d_opt_filename)){
    cat("\tD optimal file exists. Loading.\n")
    cornell_beta_2_bayesian_d_opt = readRDS(d_opt_filename)
  }else{
    stop("\tD optimal file does not exist.\n")

  }

  if(file.exists(i_opt_filename)){
    cat("\tI optimal file exists. Loading.\n")
    cornell_beta_2_bayesian_i_opt = readRDS(i_opt_filename)
  }else{
    stop("\tI optimal file does not exist.\n")
  }

  return(list(
    d_opt = cornell_beta_2_bayesian_d_opt,
    i_opt = cornell_beta_2_bayesian_i_opt,
    kappa = k
  ))

})









##########################################################################################
#### Print designs for Latex
##########################################################################################


col_names_cocktail_table = c("Choice set", "x1", "x2", "x3", "a1", "a2", "a3")

cat("Writing table for cocktail D-optimal design...")
transform_tibble_from_pseudocomp_to_comp(mnl_design_array_to_dataframe(cocktail_D_opt$X), lower_bounds_cocktail, 1:3) %>%
  select(4, 1:3, 5:7) %>%
  set_names(col_names_cocktail_table) %>%
  xtable::xtable(.,
                 caption = "D-optimal design for cocktail experiment",
                 label = "tab:cocktail_exp_d_optimal_des") %>%
  print(., include.rownames = F, file = paste0(out_folder, "res_cocktail_table_d_opt_design.tex"))
cat("Done\n\n")



cat("Writing table for cocktail I-optimal design...")
transform_tibble_from_pseudocomp_to_comp(mnl_design_array_to_dataframe(cocktail_I_opt$X), lower_bounds_cocktail, 1:3) %>%
  select(4, 1:3, 5:7) %>%
  set_names(col_names_cocktail_table) %>%
  xtable::xtable(.,
                 caption = "I-optimal design for cocktail experiment",
                 label = "tab:cocktail_exp_i_optimal_des") %>%
  print(., include.rownames = F, file = paste0(out_folder, "res_cocktail_table_i_opt_design.tex"))
cat("Done\n\n")





col_names_cornell_table = c("Choice set", "x1", "x2", "x3")
cornell_tables_filename = paste0(out_folder, "res_cornell_table_designs.tex")
for(i in seq_along(cornell_designs_untransf)){

  kappa_i = cornell_designs_untransf[[i]]$kappa

  cat("%Cornell, kappa =", kappa_i, "\n\n")
  # Cornell's D-optimal design
  cat("%Cornell's D-optimal design\n")
  mnl_design_array_to_dataframe(cornell_designs_untransf[[i]]$d_opt$X) %>%
    select(4, 1:3) %>%
    set_names(col_names_cornell_table) %>%
    xtable::xtable(
      .,
      caption = paste0("D-optimal design for artificial sweetener experiment, $\\kappa = ", kappa_i, "$"),
      label = paste0("tab:cornell_exp_d_optimal_des_kappa_", kappa_i)
      ) %>%
    print(.,
          include.rownames = F,
          file = cornell_tables_filename,
          append = T)
  cat("\n\n")



  # Cornell's I-optimal design
  cat("%Cornell's I-optimal design\n")
  mnl_design_array_to_dataframe(cornell_designs_untransf[[i]]$i_opt$X) %>%
    select(4, 1:3) %>%
    set_names(col_names_cornell_table) %>%
    xtable::xtable(
      .,
      caption = paste0("I-optimal design for artificial sweetener experiment, $\\kappa = ", kappa_i, "$"),
      label = paste0("tab:cornell_exp_i_optimal_des_kappa_", kappa_i)
      ) %>%
    print(.,
          include.rownames = F,
          file = cornell_tables_filename,
          append = T)
  cat("\n\n\n\n\n\n\n")

}





##########################################################################################
#### Plot cocktail designs
##########################################################################################


### Utility balance
(
  get_product_of_choice_probabilities(cocktail_I_opt$X, cocktail_I_opt$beta, order = 3, transform_beta = T) %>%
    mutate(Design = "I-optimal") %>%
    bind_rows(
      get_product_of_choice_probabilities(cocktail_D_opt$X, cocktail_D_opt$beta, order = 3, transform_beta = T) %>%
        mutate(Design = "D-optimal")
    ) %>%
    ggplot() +
    geom_boxplot(aes(x = Design, y = prods)) +
    theme_bw() +
    ylab("Product of choice probabilities") +
    ylim(0, 0.25)
) %>%
  ggplot2::ggsave(
    filename = paste0(out_folder, "res_cocktail_choice_probs_plot.png"),
    plot = .,
    width = 8,
    height = 9,
    units = "cm"
  )



### Distances between  alternatives within each choice set
(
  get_distances_within_choice_set(cocktail_I_opt$X) %>%
    mutate(Design = "I-optimal") %>%
    bind_rows(
      get_distances_within_choice_set(cocktail_D_opt$X) %>%
        mutate(Design = "D-optimal")
    ) %>%
    ggplot() +
    geom_boxplot(aes(x = Design, y = dist)) +
    theme_bw() +
    ylab("Distance between alternatives")
) %>%
  ggplot2::ggsave(
    filename = paste0(out_folder, "res_cocktail_distances_within_choice_set.png"),
    plot = .,
    width = 8,
    height = 9,
    units = "cm"
  )




### Design plots

levels_bayesian = c(0, .5625, 1.125, 1.6875, 2.25)
utilities_cocktail_bayesian_plot = opdesmixr:::scheffe_order3_q3_utilities(beta0, 100000) %>%
  mutate(utility_fact = cut(utility, levels_bayesian, right = F)) %>%
  mutate(utility_int = as.integer(utility_fact)) %>%
  mutate(utility_avg = case_when(
    utility_int == 1 ~ (levels_bayesian[1] + levels_bayesian[2] - levels_bayesian[1])/2,
    utility_int == 2 ~ (levels_bayesian[2] + levels_bayesian[3] - levels_bayesian[2])/2,
    utility_int == 3 ~ (levels_bayesian[3] + levels_bayesian[4] - levels_bayesian[3])/2,
    utility_int == 4 ~ (levels_bayesian[4] + levels_bayesian[5] - levels_bayesian[4])/2
  ))



# Color palette that goes from yellow to red
color_palette_cocktail = tibble(cutoffs = unique(utilities_cocktail_bayesian_plot$utility_fact)) %>%
  mutate(g = 1 - seq(0, 1, length.out = nrow(.)), r = 1, b = 0) %>%
  mutate(color = rgb(r, g, b))

# # Plot cutoffs to create legend (not useful in the end because I did it in Latex)
# color_palette_cocktail %>%
#   ggplot() +
#   geom_rect(
#     aes(xmin = r,
#         xmax = r + resolution(r),
#         ymin = g,
#         ymax = g + resolution(g),
#         fill = rgb(r, g, b)), color = "white", size = 0.1) +
#   scale_fill_identity() +
#   geom_text(aes(x = r + resolution(r)/2, y = g + resolution(g)/2, label = cutoffs)) +
#   theme_void()

# Settings for plots in PNG files
width_cocktail_simplex = 10
height_cocktail_simplex = 10
utility_point_size_cocktail = 0.01
utility_point_shape_cocktail = "circle"

ggtern::ggsave(
  filename = paste0(out_folder, "res_cocktail_design_db_simplex.png"),
  plot = plot_choice_set_utility_discrete_palette(
    cocktail_D_opt,
    utilities_cocktail_bayesian_plot,
    utility_point_size = utility_point_size_cocktail,
    utility_point_shape = utility_point_shape_cocktail,
    color_palette = color_palette_cocktail$color
  ) +
    theme(legend.position = "none"),
  width = width_cocktail_simplex,
  height = height_cocktail_simplex,
  units = "cm"
)

ggtern::ggsave(
  filename = paste0(out_folder, "res_cocktail_design_ib_simplex.png"),
  plot = plot_choice_set_utility_discrete_palette(
    cocktail_I_opt,
    utilities_cocktail_bayesian_plot,
    utility_point_size = utility_point_size_cocktail,
    utility_point_shape = utility_point_shape_cocktail,
    color_palette = color_palette_cocktail$color
  ) +
    theme(legend.position = "none"),
  width = width_cocktail_simplex,
  height = height_cocktail_simplex,
  units = "cm"
)



# ggtern::grid.arrange(
#   opdesmixr:::plot_choice_set_utility(cocktail_D_opt, utilities_cocktail_bayesian_plot) +
#     ggtitle("Bayesian D-optimal for cocktail experiment")
#   ,
#   opdesmixr:::plot_choice_set_utility(cocktail_I_opt, utilities_cocktail_bayesian_plot) +
#     ggtitle("Bayesian I-optimal for cocktail experiment")
#   ,
#   ncol = 2
# )






### FDS plots

n_points_per_alternative_cocktail = 2000

pred_vars_cocktail = mnl_get_fds_simulations(
  design_array = cocktail_I_opt$X,
  beta = cocktail_I_opt$beta,
  order = 3,
  n_points_per_alternative = n_points_per_alternative_cocktail,
  transform_beta = F,
  verbose = 1) %>%
  mutate(Design = "I-optimal") %>%
  bind_rows(
    mnl_get_fds_simulations(
      design_array = cocktail_D_opt$X,
      beta = cocktail_D_opt$beta,
      order = 3,
      n_points_per_alternative = n_points_per_alternative_cocktail,
      transform_beta = F,
      verbose = 1) %>%
      mutate(Design = "D-optimal")
  ) %>%
  bind_rows(
    mnl_get_fds_simulations(
      design_array = ruseckaite_cocktail_beta0_design,
      beta = cocktail_D_opt$beta,
      order = 3,
      n_points_per_alternative = n_points_per_alternative_cocktail,
      transform_beta = F,
      verbose = 0) %>%
      mutate(Design = "Ruseckaite et al.")
  )



cocktail_fds_plots_untransf = pred_vars_cocktail %>%
  ggplot() +
  geom_vline(xintercept = 0.5, linetype = "dashed", size = 0.2) +
  geom_hline(yintercept = pred_vars_cocktail %>%
               group_by(Design) %>%
               summarize(
                 med = median(pred_var),
                 mean = mean(pred_var)) %>%
               pull(med),
             linetype = "dashed", size = 0.2) +
  geom_line(aes(fraction, pred_var, color = Design), size = 0.8) +
  xlab("Fraction of design space") +
  ylab("Prediction variance") +
  # ggtitle("Cocktail experiment") +
  theme_bw() +
  theme(legend.position = "right") +
  scale_color_manual(values = c("red", "blue", "dark green"))


width_cocktail_fds = 18
height_cocktail_fds = 10

ggplot2::ggsave(
  filename = paste0(out_folder, "res_cocktail_fds_db_vs_ib_plot.png"),
  plot = cocktail_fds_plots_untransf,
  width = width_cocktail_fds,
  height = height_cocktail_fds,
  units = "cm"
)


# mnl_get_fds_simulations_all = function(design_array, beta, order, n_points_per_alternative = 500, transform_beta = F, verbose = 0){
#
#   q = dim(design_array)[1]
#   J = dim(design_array)[2]
#   S = dim(design_array)[3]
#
#   pred_var = suppressWarnings(as.data.frame(matrix(rep(NA_real_, nrow(beta)*J*n_points_per_alternative), ncol = J))) %>%
#     purrr::set_names(paste0("V", 1:J)) %>%
#     mutate(beta_ix = rep(1:nrow(beta), n_points_per_alternative))
#
#   progress_old = -1
#   for(k in 1:n_points_per_alternative){
#     if(verbose > 0){
#       progress = round(100*(k-1)/n_points_per_alternative, 0)
#       if(progress - progress_old >= 1){
#         cat('\r', "Progress: ", progress, "%", sep = "")
#         flush.console()
#       }
#       progress_old = progress
#     }
#     des_k = mnl_create_random_initial_design(q, J, S, seed = k)
#     vars_1 = matrix(rep(NA_real_, J*nrow(beta)), ncol = J)
#     for(j in 1:J){
#       f_x = mnl_get_Xs(des_k, 1, order = order)[j,]
#       acc = 0
#       for(i in 1:nrow(beta)){
#         inf_mat = mnl_get_information_matrix(design_array, beta = beta[i,], order = order, transform_beta = transform_beta)
#         vars_1[i,j] = t(f_x) %*% solve(inf_mat, f_x)
#       }
#
#     }
#
#     pred_var[((k-1)*nrow(beta) + 1):(k*nrow(beta)),1:J] = vars_1
#   }
#
#   # out = tibble(pred_var = sort(unlist(pred_var))) %>%
#   # mutate(fraction = 1:nrow(.)/nrow(.))
#
#   out = pred_var %>%
#     pivot_longer(cols = -beta_ix, values_to = "pred_var") %>%
#     select(-name) %>%
#     group_by(beta_ix) %>%
#     # mutate(ix = 1:n()) %>%
#     # ungroup() %>%
#     # group_by(ix) %>%
#     mutate(fraction = (0:(n()-1))/n()) %>%
#     mutate(pred_var = quantile(pred_var, probs = fraction, type = 1)) %>%
#     ungroup()
#
#   # pred_var = unlist(pred_var)
#   #
#   # out = dplyr::tibble(
#   #   fraction = (0:length(pred_var))/length(pred_var)
#   # ) %>%
#   #   dplyr::mutate(pred_var = quantile(pred_var, probs = fraction, type = 1))
#
#
#   if(verbose > 0) cat("\nFinished\n\n")
#
#   return(out)
# }
#
#
# fds_test = mnl_get_fds_simulations_all(
#   design_array = cocktail_I_opt$X,
#   beta = cocktail_I_opt$beta,
#   order = 3,
#   n_points_per_alternative = 500,
#   transform_beta = F,
#   verbose = 1)
#
# fds_test %>%
#   ggplot() +
#   geom_line(aes(fraction, pred_var, group = beta_ix), size = 0.2) +
#   xlab("Fraction of design space") +
#   ylab("Prediction variance") +
#   ggtitle("Cocktail experiment") +
#   theme_bw() +
#   theme(legend.position = "right")
#
# fds_test %>%
#   ggplot() +
#   geom_line(aes(fraction, pred_var, group = beta_ix), size = 0.2) +
#   xlab("Fraction of design space") +
#   ylab("Prediction variance") +
#   ggtitle("Cocktail experiment") +
#   theme_bw() +
#   theme(legend.position = "right") +
#   ylim(0, 4)
#
# fds_test %>%
#   group_by(beta_ix) %>%
#   mutate(ix = 1:n()) %>%
#   ungroup() %>%
#   group_by(ix, fraction) %>%
#   summarize(pred_var = mean(pred_var)) %>%
#   ungroup() %>%
#   ggplot() +
#   geom_line(aes(fraction, pred_var), size = 0.8) +
#   xlab("Fraction of design space") +
#   ylab("Prediction variance") +
#   ggtitle("Cocktail experiment") +
#   theme_bw() +
#   theme(legend.position = "right") +
#   ylim(0, 4)
#
# aaa = mnl_get_fds_simulations(
#   design_array = cocktail_I_opt$X,
#   beta = cocktail_I_opt$beta,
#   order = 3,
#   n_points_per_alternative = 500,
#   transform_beta = F,
#   verbose = 1)
#
# aaa %>%
#   ggplot() +
#   geom_line(aes(fraction, pred_var), size = 0.8) +
#   xlab("Fraction of design space") +
#   ylab("Prediction variance") +
#   ggtitle("Cocktail experiment") +
#   theme_bw() +
#   theme(legend.position = "right") +
#   ylim(0, 4)













##########################################################################################
#### Plot designs for Cornell's experiment
##########################################################################################

cornell_untrans_choice_probs = lapply(seq_along(cornell_designs_untransf), function(i){

  kappa = cornell_designs_untransf[[i]]$kappa

  get_product_of_choice_probabilities(cornell_designs_untransf[[i]]$i_opt$X,
                                      cornell_designs_untransf[[i]]$i_opt$beta,
                                      order = 3, transform_beta = F) %>%
    mutate(Design = "I-optimal") %>%
    bind_rows(
      get_product_of_choice_probabilities(cornell_designs_untransf[[i]]$d_opt$X,
                                          cornell_designs_untransf[[i]]$d_opt$beta,
                                          order = 3, transform_beta = F) %>%
        mutate(Design = "D-optimal")
    ) %>%
    mutate(kappa = kappa)
}) %>%
  bind_rows()


cornell_untrans_choice_probs_plot = cornell_untrans_choice_probs %>%
  mutate(kappa2 = paste0("kappa = ", kappa)) %>%
  mutate(kappa2 = fct_reorder(kappa2, kappa)) %>%
  ggplot() +
  geom_boxplot(aes(x = Design, y = prods)) +
  facet_wrap(~kappa2, ncol = 4) +
  ylab("Product of choice probabilities") +
  ylim(0, 0.25) +
  theme_bw()


ggplot2::ggsave(
  filename = paste0(out_folder, "res_cornell_untransf_choice_probs_plot.png"),
  plot = cornell_untrans_choice_probs_plot,
  width = 21,
  height = 6,
  units = "cm"
)









cornell_untrans_distances_within_choice_set = lapply(seq_along(cornell_designs_untransf), function(i){

  kappa = cornell_designs_untransf[[i]]$kappa

  get_distances_within_choice_set(cornell_designs_untransf[[i]]$i_opt$X) %>%
    mutate(Design = "I-optimal") %>%
    bind_rows(
      get_distances_within_choice_set(cornell_designs_untransf[[i]]$d_opt$X) %>%
        mutate(Design = "D-optimal")
    )  %>%
    mutate(kappa = kappa)


}) %>%
  bind_rows()


cornell_untrans_distances_within_choice_set_plot = cornell_untrans_distances_within_choice_set %>%
  mutate(kappa2 = paste0("kappa = ", kappa)) %>%
  mutate(kappa2 = fct_reorder(kappa2, kappa)) %>%
  ggplot() +
  geom_boxplot(aes(x = Design, y = dist)) +
  theme_bw() +
  ylab("Distance between alternatives") +
  facet_wrap(~kappa2, ncol = 4)


ggplot2::ggsave(
  filename = paste0(out_folder, "res_cornell_untransf_distances_within_choice_set.png"),
  plot = cornell_untrans_distances_within_choice_set_plot,
  width = 21,
  height = 6,
  units = "cm"
)








# Parameters for saving plots in PNG
width_cornell_simplex = 10
height_cornell_simplex = 10
utility_point_size_cornell = 0.01
utility_point_shape_cornell = "circle"

width_cornell_fds = 20
height_cornell_fds = 12

# Number of sampled points for FDS plot
fds_n_points_per_alternative_cornell = 1000

# Create utility contours
levels_cornell_bayesian = c(0, 0.375, 0.75, 1.125, 1.34)
utilities_cornell_bayesian_plot = opdesmixr:::scheffe_order3_q3_utilities(beta_2_prime, 100000) %>%
  mutate(utility_fact = cut(utility, levels_cornell_bayesian, right = F)) %>%
  mutate(utility_int = as.integer(utility_fact)) %>%
  mutate(utility_avg = case_when(
    utility_int == 1 ~ (levels_cornell_bayesian[1] + levels_cornell_bayesian[2] - levels_cornell_bayesian[1])/2,
    utility_int == 2 ~ (levels_cornell_bayesian[2] + levels_cornell_bayesian[3] - levels_cornell_bayesian[2])/2,
    utility_int == 3 ~ (levels_cornell_bayesian[3] + levels_cornell_bayesian[4] - levels_cornell_bayesian[3])/2,
    utility_int == 4 ~ (levels_cornell_bayesian[4] + levels_cornell_bayesian[5] - levels_cornell_bayesian[4])/2
  ))


# Color palette that goes from yellow to red
color_palette_cornell = tibble(cutoffs = unique(utilities_cornell_bayesian_plot$utility_fact)) %>%
  mutate(g = 1 - seq(0, 1, length.out = nrow(.)), r = 1, b = 0) %>%
  mutate(color = rgb(r, g, b))




##### Untransformed betas


## Save PNGs of designs

for(i in seq_along(cornell_designs_untransf)){

  kappa = cornell_designs_untransf[[i]]$kappa
  cat("Doing kappa =", kappa, "\n")

  d_opt_plot_untransf_i = plot_choice_set_utility_discrete_palette(
    cornell_designs_untransf[[i]]$d_opt$X, utilities_cornell_bayesian_plot,
    utility_point_size = utility_point_size_cornell,
    utility_point_shape = utility_point_shape_cornell,
    legend.position = "none",
    color_palette = color_palette_cornell$color
  ) +
    theme(plot.margin = unit(c(-5000000, -0.5, -5000000, -0.5), "cm")) # top, right, bottom, and left margins

  i_opt_plot_untransf_i = plot_choice_set_utility_discrete_palette(
    cornell_designs_untransf[[i]]$i_opt$X, utilities_cornell_bayesian_plot,
    utility_point_size = utility_point_size_cornell,
    utility_point_shape = utility_point_shape_cornell,
    legend.position = "none",
    color_palette = color_palette_cornell$color
  ) +
    theme(plot.margin = unit(c(-5000000, -0.5, -5000000, -0.5), "cm")) # top, right, bottom, and left margins


  plot_filename_basename = paste0("res_cornell_untransf_simplex_kappa_", sprintf("%03d", kappa*10), "_")

  ggtern::ggsave(
    filename = paste0(out_folder, plot_filename_basename, "D_opt.png"),
    plot = d_opt_plot_untransf_i,
    width = width_cornell_simplex,
    height = height_cornell_simplex,
    units = "cm"
  )

  ggtern::ggsave(
    filename = paste0(out_folder, plot_filename_basename, "I_opt.png"),
    plot = i_opt_plot_untransf_i,
    width = width_cornell_simplex,
    height = height_cornell_simplex,
    units = "cm"
  )

}


##### FDS plots

pred_vars_cornell_untransf = lapply(
  seq_along(cornell_designs_untransf), function(i){

    kappa_i = cornell_designs_untransf[[i]]$kappa
    cat("\n\n\nkappa:", kappa_i, format(Sys.time(),'(%H:%M:%S)'), "\n")

    beta_prior_draws_cornell = get_halton_draws(beta_2_prime, sd = sqrt(kappa_i), ndraws = 100)

    pred_vars_cornell_untransf_i = mnl_get_fds_simulations(
      design_array = cornell_designs_untransf[[i]]$i_opt$X,
      beta = beta_prior_draws_cornell,
      order = 3,
      n_points_per_alternative = fds_n_points_per_alternative_cornell,
      transform_beta = F,
      verbose = 1) %>%
      mutate(Design = "I-optimal") %>%
      bind_rows(
        mnl_get_fds_simulations(
          design_array = cornell_designs_untransf[[i]]$d_opt$X,
          beta = beta_prior_draws_cornell,
          order = 3,
          n_points_per_alternative = fds_n_points_per_alternative_cornell,
          transform_beta = F,
          verbose = 1) %>%
          mutate(Design = "D-optimal")
      ) %>%
      bind_rows(
        mnl_get_fds_simulations(
          design_array = ruseckaite_cornell_designs %>%
            filter(k == kappa_i) %>%
            select(-k) %>%
            mnl_design_dataframe_to_array(),
          beta = beta_prior_draws_cornell,
          order = 3,
          n_points_per_alternative = fds_n_points_per_alternative_cornell,
          transform_beta = F,
          verbose = 1) %>%
          mutate(Design = "Ruseckaite et al.")
      )

    return(pred_vars_cornell_untransf_i %>%
             mutate(kappa = kappa_i))

  }) %>%
  bind_rows()



cornell_fds_plots_untransf = pred_vars_cornell_untransf %>%
  left_join(
    pred_vars_cornell_untransf %>%
      group_by(Design, kappa) %>%
      summarize(
        med = median(pred_var),
        mean = mean(pred_var))
  ) %>%
  mutate(kappa2 = paste0("kappa = ", kappa)) %>%
  mutate(kappa2 = fct_reorder(kappa2, kappa)) %>%
  ggplot() +
  geom_vline(xintercept = 0.5, linetype = "dashed", size = 0.2) +
  geom_hline(aes(yintercept = med), linetype = "dashed", size = 0.2) +
  geom_line(aes(fraction, pred_var, color = Design), size = 0.8) +
  xlab("Fraction of design space") +
  ylab("Prediction variance") +
  # ggtitle("Cornell's experiment") +
  theme_bw() +
  theme(legend.position = "bottom") +
  facet_wrap(~kappa2, scales = "free_y") +
  scale_color_manual(values = c("red", "blue", "dark green"))





ggplot2::ggsave(
  filename = paste0(out_folder, "res_cornell_untransf_fds_db_vs_ib_plot.png"),
  plot = cornell_fds_plots_untransf,
  width = width_cornell_fds,
  height = height_cornell_fds,
  units = "cm"
)



# print(cornell_fds_plots_untransf)

#
# ## Create a list with the plots of the optimal designs
# cornell_optimal_designs_plots_untransf = lapply(seq_along(cornell_designs_untransf), function(i){
#
#   out_list = list(
#     d_opt_plot = opdesmixr:::plot_choice_set_utility(
#       cornell_designs_untransf[[i]]$d_opt$X, utilities_cornell_bayesian_plot,
#       legend.position = "none"
#     ) +
#       ggtitle("", subtitle = paste0("D-optimal, kappa = ", cornell_designs_untransf[[i]]$kappa))
#     ,
#     i_opt_plot = opdesmixr:::plot_choice_set_utility(
#       cornell_designs_untransf[[i]]$i_opt$X, utilities_cornell_bayesian_plot,
#       legend.position = "none"
#     ) +
#       ggtitle("", subtitle = paste0("I-optimal, kappa = ", cornell_designs_untransf[[i]]$kappa))
#
#     ,
#     kappa = cornell_designs_untransf[[i]]$kappa
#   )
#
#   return(out_list)
#
#
# })
# # Couldn't get it to print with the list only
#
# ggtern::grid.arrange(
#   cornell_optimal_designs_plots_untransf[[1]][[1]], cornell_optimal_designs_plots_untransf[[1]][[2]],
#   cornell_optimal_designs_plots_untransf[[2]][[1]], cornell_optimal_designs_plots_untransf[[2]][[2]],
#   cornell_optimal_designs_plots_untransf[[3]][[1]], cornell_optimal_designs_plots_untransf[[3]][[2]],
#   cornell_optimal_designs_plots_untransf[[4]][[1]], cornell_optimal_designs_plots_untransf[[4]][[2]],
#   top = grid::textGrob(paste0("kappa = 0.5, 5, 10, 30"), gp = grid::gpar(fontsize = 15)),
#   ncol = 2)





# ##### Transformed betas
#
# ## Save PNGs of designs
#
# for(i in seq_along(cornell_designs_transf)){
#
#   kappa = cornell_designs_transf[[i]]$kappa
#   cat("Doing kappa =", kappa, "\n")
#
#   d_opt_plot_transf_i = plot_choice_set_utility_discrete_palette(
#     cornell_designs_transf[[i]]$d_opt$X, utilities_cornell_bayesian_plot,
#     utility_point_size = utility_point_size_cornell,
#     utility_point_shape = utility_point_shape_cornell,
#     legend.position = "none",
#     color_palette = color_palette_cornell$color
#   ) +
#     theme(plot.margin = unit(c(-5000000, -0.5, -5000000, -0.5), "cm")) # top, right, bottom, and left margins
#
#   i_opt_plot_transf_i = plot_choice_set_utility_discrete_palette(
#     cornell_designs_transf[[i]]$i_opt$X, utilities_cornell_bayesian_plot,
#     utility_point_size = utility_point_size_cornell,
#     utility_point_shape = utility_point_shape_cornell,
#     legend.position = "none",
#     color_palette = color_palette_cornell$color
#   ) +
#     theme(plot.margin = unit(c(-5000000, -0.5, -5000000, -0.5), "cm")) # top, right, bottom, and left margins
#
#
#   plot_filename_basename = paste0("res_cornell_transf_simplex_kappa_", sprintf("%03d", kappa*10), "_")
#
#   ggtern::ggsave(
#     filename = paste0(out_folder, plot_filename_basename, "D_opt.png"),
#     plot = d_opt_plot_transf_i,
#     width = width_cornell_simplex,
#     height = height_cornell_simplex,
#     units = "cm"
#   )
#
#   ggtern::ggsave(
#     filename = paste0(out_folder, plot_filename_basename, "I_opt.png"),
#     plot = i_opt_plot_transf_i,
#     width = width_cornell_simplex,
#     height = height_cornell_simplex,
#     units = "cm"
#   )
# }
#
#
# ##### FDS plots
#
# pred_vars_cornell_transf = lapply(
#   seq_along(cornell_designs_transf), function(i){
#
#     kappa_i = cornell_designs_transf[[i]]$kappa
#     cat("\n\n\nkappa:", kappa_i, format(Sys.time(),'(%H:%M:%S)'), "\n")
#
#     pred_vars_cornell_transf_i = mnl_get_fds_simulations(
#       design_array = cornell_designs_transf[[i]]$i_opt$X,
#       beta = cornell_designs_transf[[i]]$i_opt$beta,
#       order = 3,
#       n_points_per_alternative = fds_n_points_per_alternative_cornell,
#       transform_beta = T,
#       verbose = 1) %>%
#       mutate(Design = "I-optimal") %>%
#       bind_rows(
#         mnl_get_fds_simulations(
#           design_array = cornell_designs_transf[[i]]$d_opt$X,
#           beta = cornell_designs_transf[[i]]$d_opt$beta,
#           order = 3,
#           n_points_per_alternative = fds_n_points_per_alternative_cornell,
#           transform_beta = T,
#           verbose = 1) %>%
#           mutate(Design = "D-optimal")
#       )
#
#     return(pred_vars_cornell_transf_i %>%
#              mutate(kappa = kappa_i))
#
#   }) %>%
#   bind_rows()
#
#
#
# cornell_fds_plots_transf = pred_vars_cornell_transf %>%
#   left_join(
#     pred_vars_cornell_transf %>%
#       group_by(Design, kappa) %>%
#       summarize(
#         med = median(pred_var),
#         mean = mean(pred_var))
#   ) %>%
#   mutate(kappa2 = paste0("kappa = ", kappa)) %>%
#   mutate(kappa2 = fct_reorder(kappa2, kappa)) %>%
#   ggplot() +
#   geom_vline(xintercept = 0.5, linetype = "dashed", size = 0.2) +
#   geom_hline(aes(yintercept = med), linetype = "dashed", size = 0.2) +
#   geom_line(aes(fraction, pred_var, color = Design), size = 0.8) +
#   xlab("Fraction of design space") +
#   ylab("Prediction variance") +
#   # ggtitle("Cornell's experiment") +
#   theme_bw() +
#   theme(legend.position = "bottom") +
#   facet_wrap(~kappa2, scales = "free_y") +
#   scale_color_manual(values = c("red", "blue", "dark green"))
#
#
#
#
#
# ggplot2::ggsave(
#   filename = paste0(out_folder, "res_cornell_transf_fds_db_vs_ib_plot.png"),
#   plot = cornell_fds_plots_transf,
#   width = width_cornell_fds,
#   height = height_cornell_fds,
#   units = "cm"
# )
#
#
#
#
#
#
#
#
#
#
