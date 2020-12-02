library(opdesmixr)
library(tidyverse)
library(here)


designs_folder = here("inst/misc_output/cocktail_cornell_designs/")
dir.create(designs_folder, showWarnings = F)


n_cores = parallel::detectCores()

##########################################################################################
#### Cocktail experiment
##########################################################################################

q = 3
J = 2
S = 480
n_draws = 128

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


cocktail_d_opt_filename = paste0(designs_folder, "cocktail_d_optimal_480_choice_sets.rds")
cocktail_i_opt_filename = paste0(designs_folder, "cocktail_i_optimal_480_choice_sets.rds")

##########################################################################################
#### Load or create the designs
##########################################################################################

# 14 iterations took around 13 hours with a 12-core PC
if(file.exists(cocktail_d_opt_filename)){
  cat("D_B optimal design already exists.\n")
} else{
  # 10 mins
  (t1D = Sys.time())
  cocktail_D_opt = mnl_mixture_coord_exch(
    n_random_starts = n_rand_starts,
    q = q,
    J = J,
    S = S,
    beta = beta_correlated_draws_cocktail,
    transform_beta = F,
    opt_method = "B",
    opt_crit = "D",
    max_it = max_it,
    verbose = 1,
    plot_designs = F,
    seed = seed,
    n_cores = n_cores
  )
  (t2D = Sys.time())
  t2D - t1D

  saveRDS(cocktail_D_opt, cocktail_d_opt_filename)
}


# 10 iterations took around 10.5 hours with a 12-core PC
if(file.exists(cocktail_i_opt_filename)){
  cat("I_B optimal design already exists.\n")
} else{
  # 14 mins
  (t1I = Sys.time())
  cocktail_I_opt =  mnl_mixture_coord_exch(
    n_random_starts = n_rand_starts,
    q = q,
    J = J,
    S = S,
    beta = beta_correlated_draws_cocktail,
    transform_beta = F,
    opt_method = "B",
    opt_crit = "I",
    max_it = max_it,
    verbose = 1,
    plot_designs = F,
    seed = seed,
    n_cores = n_cores
  )
  (t2I = Sys.time())
  t2I - t1I

  saveRDS(cocktail_I_opt, cocktail_i_opt_filename)
}


