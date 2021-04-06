library(opdesmixr)
library(tidyverse)
library(here)



plot_array = function(X){
  dim_X = dim(X)
  q = dim_X[1]
  S = dim_X[3]

  if(q != 3) stop("Design must be of 3 ingredients.")

  # Convert 3 dimensional arrays into matrices by vertically binding them
  design_df = mnl_design_array_to_dataframe(X)

  out_plot = design_df %>%
    ggtern::ggtern(ggtern::aes(c1, c2, c3)) +
    geom_point(shape = "x", size = 4) +
    theme_minimal() +
    ggtern::theme_nomask()

  return(out_plot)
}




designs_folder = here("inst/misc_output/cocktail_cornell_designs/")
dir.create(designs_folder, showWarnings = F)


q = 3
J = 2
S = 480
n_draws = 5000

beta0 = c(1.36, 1.57, 2.47, -0.43, 0.50, 1.09)

sigma0 = matrix(
  c(6.14, 5.00, 2.74, -0.43, -2.81, -3.33,
    5.00, 6.76, 4.47, -1.79, -6.13, -3.51,
    2.74, 4.47, 3.45, -1.38, -4.71, -2.17,
    -0.43, -1.79, -1.38, 1.18, 2.39, 0.71,
    -2.81, -6.13, -4.71, 2.39, 7.43, 2.71,
    -3.33, -3.51, -2.17, 0.71, 2.71, 2.49),
  ncol = 6,
  byrow = T)

beta_correlated_draws_cocktail = get_correlated_halton_draws(beta0, sigma0, n_draws)


n_rand_starts = 12
# max_it = 15
max_it = 10
seed = 2020


cocktail_d_opt_filename_small = paste0(designs_folder, "cocktail_d_optimal.rds")
cocktail_i_opt_filename_small = paste0(designs_folder, "cocktail_i_optimal.rds")
cocktail_d_opt_filename_480 = paste0(designs_folder, "cocktail_d_optimal_480_choice_sets.rds")
cocktail_i_opt_filename_480 = paste0(designs_folder, "cocktail_i_optimal_480_choice_sets.rds")


cocktail_d_opt_small = readRDS(cocktail_d_opt_filename_small)
cocktail_i_opt_small = readRDS(cocktail_i_opt_filename_small)
cocktail_d_opt_480 = readRDS(cocktail_d_opt_filename_480)
cocktail_i_opt_480 = readRDS(cocktail_i_opt_filename_480)

cocktail_d_opt_480_X = cocktail_d_opt_480$X
cocktail_i_opt_480_X = cocktail_i_opt_480$X

dim_480 = c(3, 2, 480)
cocktail_d_opt_small_replicated = array(0.0, dim = dim_480)
cocktail_i_opt_small_replicated = array(0.0, dim = dim_480)

n = 16
for(i in 1:(480/16)){
  lo = (i-1)*n + 1
  hi = (i)*n
  cocktail_d_opt_small_replicated[,,lo:hi] = cocktail_d_opt_small$X
  cocktail_i_opt_small_replicated[,,lo:hi] = cocktail_i_opt_small$X
}


mnl_get_opt_crit_value(X = cocktail_d_opt_480_X, beta = beta_correlated_draws_cocktail, order = 3, transform_beta = F, opt_crit = "D")
mnl_get_opt_crit_value(X = cocktail_d_opt_small_replicated, beta = beta_correlated_draws_cocktail, order = 3, transform_beta = F, opt_crit = "D")
mnl_get_opt_crit_value(X = cocktail_d_opt_480_X, beta = beta_correlated_draws_cocktail, order = 3, transform_beta = F, opt_crit = "I")
mnl_get_opt_crit_value(X = cocktail_d_opt_small_replicated, beta = beta_correlated_draws_cocktail, order = 3, transform_beta = F, opt_crit = "I")

mnl_get_opt_crit_value(X = cocktail_i_opt_480_X, beta = beta_correlated_draws_cocktail, order = 3, transform_beta = F, opt_crit = "D")
mnl_get_opt_crit_value(X = cocktail_i_opt_small_replicated, beta = beta_correlated_draws_cocktail, order = 3, transform_beta = F, opt_crit = "D")
mnl_get_opt_crit_value(X = cocktail_i_opt_480_X, beta = beta_correlated_draws_cocktail, order = 3, transform_beta = F, opt_crit = "I")
mnl_get_opt_crit_value(X = cocktail_i_opt_small_replicated, beta = beta_correlated_draws_cocktail, order = 3, transform_beta = F, opt_crit = "I")





ggtern::grid.arrange(

  plot_array(cocktail_d_opt_480_X) +
    ggtitle(
      "D-optimal for 480 choice sets",
      subtitle = paste0(
        "D-optimality = ",
        round(mnl_get_opt_crit_value(X = cocktail_d_opt_480_X, beta = beta_correlated_draws_cocktail, order = 3, transform_beta = F, opt_crit = "D"), 2),
        ", ",
        "I-optimality = ",
        round(mnl_get_opt_crit_value(X = cocktail_d_opt_480_X, beta = beta_correlated_draws_cocktail, order = 3, transform_beta = F, opt_crit = "I"), 2)
        )
      )
  ,
  plot_array(cocktail_d_opt_small_replicated) +
    ggtitle(
      "D-optimal for 16 choice sets replicated 30 times",
      subtitle = paste0(
        "D-optimality = ",
        round(mnl_get_opt_crit_value(X = cocktail_d_opt_small_replicated, beta = beta_correlated_draws_cocktail, order = 3, transform_beta = F, opt_crit = "D"), 2),
        ", ",
        "I-optimality = ",
        round(mnl_get_opt_crit_value(X = cocktail_d_opt_small_replicated, beta = beta_correlated_draws_cocktail, order = 3, transform_beta = F, opt_crit = "I"), 2)
      )
    )

  ,

  plot_array(cocktail_i_opt_480_X) +
    ggtitle(
      "I-optimal for 480 choice sets",
      subtitle = paste0(
        "D-optimality = ",
        round(mnl_get_opt_crit_value(X = cocktail_i_opt_480_X, beta = beta_correlated_draws_cocktail, order = 3, transform_beta = F, opt_crit = "D"), 2),
        ", ",
        "I-optimality = ",
        round(mnl_get_opt_crit_value(X = cocktail_i_opt_480_X, beta = beta_correlated_draws_cocktail, order = 3, transform_beta = F, opt_crit = "I"), 2)
      )
    )
  ,
  plot_array(cocktail_i_opt_small_replicated) +
    ggtitle(
      "I-optimal for 16 choice sets replicated 30 times",
      subtitle = paste0(
        "D-optimality = ",
        round(mnl_get_opt_crit_value(X = cocktail_i_opt_small_replicated, beta = beta_correlated_draws_cocktail, order = 3, transform_beta = F, opt_crit = "D"), 2),
        ", ",
        "I-optimality = ",
        round(mnl_get_opt_crit_value(X = cocktail_i_opt_small_replicated, beta = beta_correlated_draws_cocktail, order = 3, transform_beta = F, opt_crit = "I"), 2)
      )
    )

  ,
  ncol = 2)

