
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
  design_object, utility_data, design_point_size = 2, design_point_alpha = 1,
  utility_point_size = 0.3, utility_point_shape = "square", utility_point_alpha = 1.0,
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
      shape = utility_point_shape,
      alpha = utility_point_alpha
    ) +
    ggplot2::scale_color_manual(values = color_palette) +
    ggplot2::geom_point(data = X_final_tbl, ggtern::aes(c1, c2, c3, shape = choice_set), color = "black", size = design_point_size, stroke = 1.5, inherit.aes = F, alpha = design_point_alpha) +
    ggtern::theme_nomask() +
    scale_shape_manual(values = 1:length(unique(X_final_tbl$choice_set))) +
    theme(legend.position = legend.position, legend.box = legend.box)
  return(out_plot)
}






# Transforms from pseudocomponent x to original component a.
# x is a scalar or a vector, l is a scalar denoting the lower bound.
# sum_l is a scalar denoting the sum of the lower bounds.
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
#### Setup
##########################################################################################

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

out_folder = here("inst/misc_output/cocktail_cornell_plots_and_tables/")
dir.create(out_folder, showWarnings = F)


# Folder with the designs by Ruseckaite
ruseckaite_designs_folder = opdesmixr:::get_file_path_inst("exdata/")





























# Cocktail experiment -----------------------------------------------------
############################################################################
############################################################################
############################################################################




#### Load designs for cocktail experiment

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


###############################
#### Print designs for Latex
###############################

col_names_cocktail_table = c("Choice set", "x1", "x2", "x3", "a1", "a2", "a3")

cat("Writing table for cocktail Bayesian D-optimal design...")
transform_tibble_from_pseudocomp_to_comp(mnl_design_array_to_dataframe(cocktail_D_opt$X), lower_bounds_cocktail, 1:3) %>%
  select(4, 1:3, 5:7) %>%
  set_names(col_names_cocktail_table) %>%
  xtable::xtable(.,
                 caption = "Bayesian D-optimal design for cocktail experiment",
                 label = "tab:cocktail_exp_d_optimal_des") %>%
  print(., include.rownames = F, file = paste0(out_folder, "res_cocktail_table_d_opt_design.tex"))
cat("Done\n\n")



cat("Writing table for cocktail Bayesian I-optimal design...")
transform_tibble_from_pseudocomp_to_comp(mnl_design_array_to_dataframe(cocktail_I_opt$X), lower_bounds_cocktail, 1:3) %>%
  select(4, 1:3, 5:7) %>%
  set_names(col_names_cocktail_table) %>%
  xtable::xtable(.,
                 caption = "Bayesian I-optimal design for cocktail experiment",
                 label = "tab:cocktail_exp_i_optimal_des") %>%
  print(., include.rownames = F, file = paste0(out_folder, "res_cocktail_table_i_opt_design.tex"))
cat("Done\n\n")



###############################
#### Plot utility balance
###############################

(
  get_product_of_choice_probabilities(cocktail_I_opt$X, cocktail_I_opt$beta, order = 3, transform_beta = T) %>%
    mutate(Design = "Bayesian\nI-optimal") %>%
    bind_rows(
      get_product_of_choice_probabilities(cocktail_D_opt$X, cocktail_D_opt$beta, order = 3, transform_beta = T) %>%
        mutate(Design = "Bayesian\nD-optimal")
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



###############################
#### Plot distances between alternatives within each choice set
###############################

(
  get_distances_within_choice_set(cocktail_I_opt$X) %>%
    mutate(Design = "Bayesian\nI-optimal") %>%
    bind_rows(
      get_distances_within_choice_set(cocktail_D_opt$X) %>%
        mutate(Design = "Bayesian\nD-optimal")
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




###############################
### Cocktail design plots
###############################

levels_cocktail = c(0, .5625, 1.125, 1.6875, 2.25)
utilities_cocktail_bayesian_plot = opdesmixr:::scheffe_order3_q3_utilities(beta0, 100000) %>%
  mutate(utility_fact = cut(utility, levels_cocktail, right = F)) %>%
  mutate(utility_int = as.integer(utility_fact)) %>%
  mutate(utility_avg = case_when(
    utility_int == 1 ~ (levels_cocktail[1] + levels_cocktail[2] - levels_cocktail[1])/2,
    utility_int == 2 ~ (levels_cocktail[2] + levels_cocktail[3] - levels_cocktail[2])/2,
    utility_int == 3 ~ (levels_cocktail[3] + levels_cocktail[4] - levels_cocktail[3])/2,
    utility_int == 4 ~ (levels_cocktail[4] + levels_cocktail[5] - levels_cocktail[4])/2
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
    theme(legend.position = "none") +
    labs(x = "x_1", y = "x_2", z = "x_3") +
    ggtern::theme_latex(TRUE),
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
    theme(legend.position = "none") +
    labs(x = "x_1", y = "x_2", z = "x_3") +
    ggtern::theme_latex(TRUE),
  width = width_cocktail_simplex,
  height = height_cocktail_simplex,
  units = "cm"
)


###############################
### Cocktail FDS plots
###############################

n_points_per_alternative_cocktail = 2000

pred_vars_cocktail = mnl_get_fds_simulations(
  design_array = cocktail_I_opt$X,
  beta = cocktail_I_opt$beta,
  order = 3,
  n_points_per_alternative = n_points_per_alternative_cocktail,
  transform_beta = F,
  verbose = 1) %>%
  mutate(Design = "Bayesian I-optimal") %>%
  bind_rows(
    mnl_get_fds_simulations(
      design_array = cocktail_D_opt$X,
      beta = cocktail_D_opt$beta,
      order = 3,
      n_points_per_alternative = n_points_per_alternative_cocktail,
      transform_beta = F,
      verbose = 1) %>%
      mutate(Design = "Bayesian D-optimal")
  ) %>%
  bind_rows(
    mnl_get_fds_simulations(
      design_array = ruseckaite_cocktail_beta0_design,
      beta = cocktail_D_opt$beta,
      order = 3,
      n_points_per_alternative = n_points_per_alternative_cocktail,
      transform_beta = F,
      verbose = 1) %>%
      mutate(Design = "Ruseckaite et al.")
  )



cocktail_fds_plots = pred_vars_cocktail %>%
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
  plot = cocktail_fds_plots,
  width = width_cocktail_fds,
  height = height_cocktail_fds,
  units = "cm"
)
























# Cornell's experiment ----------------------------------------------------
############################################################################
############################################################################
############################################################################

ruseckaite_cornell_designs = read_csv(paste0(ruseckaite_designs_folder, "/ruseckaite_cornell_designs.csv"))
kappas = c(0.5, 5, 10, 30)

beta_2 = c(0.86, 0.21, 0, 3.07, 2.34, 3.24, -20.59)
beta_2_prime = c(0.86, 0.21, 3.07, 2.34, 3.24, -20.59)

cornell_designs_basefilename_transf = paste0(designs_folder, "cornell_experiment_transformed_betas_maxit20")
cornell_designs_basefilename_untransf = paste0(designs_folder, "cornell_experiment_untransformed_betas_maxit20")
cornell_designs_basefilename_analytic_transf = paste0(designs_folder, "cornell_experiment_analytic_transformed_betas_maxit20")









##############################################################
##############################################################
##### Analytically transformed betas
##############################################################
##############################################################


## Read designs with analytic_transformed betas and save them in a list
cornell_designs_analytic_transf = lapply(kappas, function(k){

  cat("kappa =", k, "\n")

  d_opt_filename = paste0(cornell_designs_basefilename_analytic_transf, "_kappa", k, "_Dopt.rds")
  i_opt_filename = paste0(cornell_designs_basefilename_analytic_transf, "_kappa", k, "_Iopt.rds")



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






###############################
### D-optimalities of our design and Ruseckaite's
###############################

cornell_d_optimalities_analytic_transf = lapply(
  seq_along(cornell_designs_analytic_transf), function(i){

    kappa_i = cornell_designs_analytic_transf[[i]]$kappa
    cat("\n\n\nkappa:", kappa_i, format(Sys.time(),'(%H:%M:%S)'), "\n")

    # beta_prior_draws_cornell = get_halton_draws(beta_2_prime, sd = sqrt(kappa_i), ndraws = 100)
    Sigma_prime = transform_varcov_matrix(kappa_i*diag(7), 3)
    beta_prior_draws_cornell = get_correlated_halton_draws(beta = beta_2_prime, sigma = Sigma_prime, n_draws = 512)

    our_d_opt_design = cornell_designs_analytic_transf[[i]]$d_opt$X

    ruseckaite_d_opt_design = ruseckaite_cornell_designs %>%
      filter(k == kappa_i) %>%
      select(-k) %>%
      mnl_design_dataframe_to_array()

    our_d_optimality = mnl_get_opt_crit_value(X = our_d_opt_design, beta = beta_prior_draws_cornell, order = 3, opt_crit = "D", transform_beta = F)
    rus_d_optimality = mnl_get_opt_crit_value(X = ruseckaite_d_opt_design, beta = beta_prior_draws_cornell, order = 3, opt_crit = "D", transform_beta = F)

    return(
      tibble(
        our_d_optimality = our_d_optimality,
        rus_d_optimality = rus_d_optimality,
        kappa = kappa_i
      )
    )

  }) %>%
  bind_rows()


write_csv(cornell_d_optimalities_analytic_transf, file = paste0(out_folder, "res_cornell_d_optimalities_analytic_transf.csv"))



###############################
#### Print designs for Latex
###############################

col_names_cornell_analytic_transf_table = c("Choice set", "x1", "x2", "x3")
cornell_analytic_transf_tables_filename = paste0(out_folder, "res_cornell_analytic_transf_table_designs.tex")
if(file.exists(cornell_analytic_transf_tables_filename)){
  stop("Tex file already exists!")
} else{
  for(i in seq_along(cornell_designs_analytic_transf)){

    kappa_i = cornell_designs_analytic_transf[[i]]$kappa

    cat("% Saving optimal design tables for Cornell's experiment, kappa =", kappa_i, "\n")


    # # Cornell's Bayesian D-optimal design
    # cat("%Cornell's Bayesian D-optimal design\n")
    mnl_design_array_to_dataframe(cornell_designs_analytic_transf[[i]]$d_opt$X) %>%
      select(4, 1:3) %>%
      set_names(col_names_cornell_analytic_transf_table) %>%
      xtable::xtable(
        .,
        caption = paste0("Bayesian D-optimal design for artificial sweetener experiment, $\\kappa = ", kappa_i, "$"),
        label = paste0("tab:cornell_exp_d_optimal_des_kappa_", kappa_i)
      ) %>%
      print(.,
            include.rownames = F,
            file = cornell_analytic_transf_tables_filename,
            append = T)
    # cat("\n\n")



    # # Cornell's Bayesian I-optimal design
    # cat("%Cornell's Bayesian I-optimal design\n")
    mnl_design_array_to_dataframe(cornell_designs_analytic_transf[[i]]$i_opt$X) %>%
      select(4, 1:3) %>%
      set_names(col_names_cornell_analytic_transf_table) %>%
      xtable::xtable(
        .,
        caption = paste0("Bayesian I-optimal design for artificial sweetener experiment, $\\kappa = ", kappa_i, "$"),
        label = paste0("tab:cornell_exp_i_optimal_des_kappa_", kappa_i)
      ) %>%
      print(.,
            include.rownames = F,
            file = cornell_analytic_transf_tables_filename,
            append = T)
    # cat("\n\n\n\n\n\n\n")

  }
}






###############################
## Choice probability plots
###############################

cornell_analytic_trans_choice_probs = lapply(seq_along(cornell_designs_analytic_transf), function(i){

  kappa = cornell_designs_analytic_transf[[i]]$kappa

  get_product_of_choice_probabilities(cornell_designs_analytic_transf[[i]]$i_opt$X,
                                      cornell_designs_analytic_transf[[i]]$i_opt$beta,
                                      order = 3, transform_beta = F) %>%
    mutate(Design = "Bayesian\nI-optimal") %>%
    bind_rows(
      get_product_of_choice_probabilities(cornell_designs_analytic_transf[[i]]$d_opt$X,
                                          cornell_designs_analytic_transf[[i]]$d_opt$beta,
                                          order = 3, transform_beta = F) %>%
        mutate(Design = "Bayesian\nD-optimal")
    ) %>%
    mutate(kappa = kappa)
}) %>%
  bind_rows()


cornell_analytic_trans_choice_probs_plot = cornell_analytic_trans_choice_probs %>%
  mutate(kappa2 = paste0("kappa == ", kappa)) %>%
  mutate(kappa2 = fct_reorder(kappa2, kappa)) %>%
  ggplot() +
  geom_boxplot(aes(x = Design, y = prods)) +
  facet_wrap(~kappa2, ncol = 4, labeller = label_parsed) +
  ylab("Product of choice probabilities") +
  ylim(0, 0.25) +
  theme_bw()




ggplot2::ggsave(
  filename = paste0(out_folder, "res_cornell_analytic_transf_choice_probs_plot.png"),
  plot = cornell_analytic_trans_choice_probs_plot,
  width = 21,
  height = 6,
  units = "cm"
)







###############################
## Distance plots
###############################

cornell_analytic_trans_distances_within_choice_set = lapply(seq_along(cornell_designs_analytic_transf), function(i){

  kappa = cornell_designs_analytic_transf[[i]]$kappa

  get_distances_within_choice_set(cornell_designs_analytic_transf[[i]]$i_opt$X) %>%
    mutate(Design = "Bayesian\nI-optimal") %>%
    bind_rows(
      get_distances_within_choice_set(cornell_designs_analytic_transf[[i]]$d_opt$X) %>%
        mutate(Design = "Bayesian\nD-optimal")
    )  %>%
    mutate(kappa = kappa)


}) %>%
  bind_rows()


cornell_analytic_trans_distances_within_choice_set_plot = cornell_analytic_trans_distances_within_choice_set %>%
  mutate(kappa2 = paste0("kappa == ", kappa)) %>%
  mutate(kappa2 = fct_reorder(kappa2, kappa)) %>%
  ggplot() +
  geom_boxplot(aes(x = Design, y = dist)) +
  theme_bw() +
  ylab("Distance between alternatives") +
  facet_wrap(~kappa2, ncol = 4, labeller = label_parsed)


ggplot2::ggsave(
  filename = paste0(out_folder, "res_cornell_analytic_transf_distances_within_choice_set.png"),
  plot = cornell_analytic_trans_distances_within_choice_set_plot,
  width = 21,
  height = 6,
  units = "cm"
)





###############################
#### Plot designs for Cornell's experiment
###############################


# Parameters for saving plots in PNG
width_cornell_simplex = 10
height_cornell_simplex = 10
utility_point_size_cornell = 0.01
utility_point_shape_cornell = "circle"

width_cornell_fds = 20
height_cornell_fds = 12

# Number of sampled points for FDS plot
# fds_n_points_per_alternative_cornell_analytic_transf = 100
fds_n_points_per_alternative_cornell_analytic_transf = 1000

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



for(i in seq_along(cornell_designs_analytic_transf)){

  kappa = cornell_designs_analytic_transf[[i]]$kappa
  cat("Doing kappa =", kappa, "\n")

  d_opt_plot_analytic_transf_i = plot_choice_set_utility_discrete_palette(
    cornell_designs_analytic_transf[[i]]$d_opt$X, utilities_cornell_bayesian_plot,
    utility_point_size = utility_point_size_cornell,
    utility_point_shape = utility_point_shape_cornell,
    legend.position = "none",
    color_palette = color_palette_cornell$color
  ) +
    labs(x = "x_1", y = "x_2", z = "x_3") +
    ggtern::theme_latex(TRUE) +
    theme(plot.margin = unit(c(-5000000, -0.5, -5000000, -0.5), "cm")) # top, right, bottom, and left margins

  i_opt_plot_analytic_transf_i = plot_choice_set_utility_discrete_palette(
    cornell_designs_analytic_transf[[i]]$i_opt$X, utilities_cornell_bayesian_plot,
    utility_point_size = utility_point_size_cornell,
    utility_point_shape = utility_point_shape_cornell,
    legend.position = "none",
    color_palette = color_palette_cornell$color
  ) +
    labs(x = "x_1", y = "x_2", z = "x_3") +
    ggtern::theme_latex(TRUE) +
    theme(plot.margin = unit(c(-5000000, -0.5, -5000000, -0.5), "cm")) # top, right, bottom, and left margins

  plot_filename_basename = paste0("res_cornell_analytic_transf_simplex_kappa_", sprintf("%03d", kappa*10), "_")

  ggtern::ggsave(
    filename = paste0(out_folder, plot_filename_basename, "D_opt.png"),
    plot = d_opt_plot_analytic_transf_i,
    width = width_cornell_simplex,
    height = height_cornell_simplex,
    units = "cm"
  )

  ggtern::ggsave(
    filename = paste0(out_folder, plot_filename_basename, "I_opt.png"),
    plot = i_opt_plot_analytic_transf_i,
    width = width_cornell_simplex,
    height = height_cornell_simplex,
    units = "cm"
  )
}



###############################
##### FDS plots
###############################



pred_vars_cornell_analytic_transf = lapply(
  seq_along(cornell_designs_analytic_transf), function(i){

    kappa_i = cornell_designs_analytic_transf[[i]]$kappa
    cat("\n\n\nkappa:", kappa_i, format(Sys.time(),'(%H:%M:%S)'), "\n")

    # beta_prior_draws_cornell = get_halton_draws(beta_2_prime, sd = sqrt(kappa_i), ndraws = 128)
    Sigma_prime = transform_varcov_matrix(kappa_i*diag(7), 3)
    beta_prior_draws_cornell = get_correlated_halton_draws(beta = beta_2_prime, sigma = Sigma_prime, n_draws = 128)

    pred_vars_cornell_analytic_transf_i = mnl_get_fds_simulations(
      design_array = cornell_designs_analytic_transf[[i]]$i_opt$X,
      beta = beta_prior_draws_cornell,
      order = 3,
      n_points_per_alternative = fds_n_points_per_alternative_cornell_analytic_transf,
      transform_beta = F,
      verbose = 1) %>%
      mutate(Design = "Bayesian I-optimal") %>%
      bind_rows(
        mnl_get_fds_simulations(
          design_array = cornell_designs_analytic_transf[[i]]$d_opt$X,
          beta = beta_prior_draws_cornell,
          order = 3,
          n_points_per_alternative = fds_n_points_per_alternative_cornell_analytic_transf,
          transform_beta = F,
          verbose = 1) %>%
          mutate(Design = "Bayesian D-optimal")
      ) %>%
      bind_rows(
        mnl_get_fds_simulations(
          design_array = ruseckaite_cornell_designs %>%
            filter(k == kappa_i) %>%
            select(-k) %>%
            mnl_design_dataframe_to_array(),
          beta = beta_prior_draws_cornell,
          order = 3,
          n_points_per_alternative = fds_n_points_per_alternative_cornell_analytic_transf,
          transform_beta = F,
          verbose = 1) %>%
          mutate(Design = "Ruseckaite et al.")
      )

    return(pred_vars_cornell_analytic_transf_i %>%
             mutate(kappa = kappa_i))

  }) %>%
  bind_rows()



cornell_fds_plots_analytic_transf = pred_vars_cornell_analytic_transf %>%
  left_join(
    pred_vars_cornell_analytic_transf %>%
      group_by(Design, kappa) %>%
      summarize(
        med = median(pred_var),
        mean = mean(pred_var))
  ) %>%
  mutate(kappa2 = paste0("kappa == ", kappa)) %>%
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
  facet_wrap(~kappa2, ncol = 2, scales = "free_y", labeller = label_parsed) +
  scale_color_manual(values = c("red", "blue", "dark green"))




ggplot2::ggsave(
  filename = paste0(out_folder, "res_cornell_analytic_transf_fds_db_vs_ib_plot.png"),
  plot = cornell_fds_plots_analytic_transf,
  width = width_cornell_fds,
  height = height_cornell_fds,
  units = "cm"
)




















