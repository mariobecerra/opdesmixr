library(opdesmixr)
library(tidyverse)
library(here)

# Assuming this sript is run before installing.
designs_folder = here("inst/misc_output/cocktail_cornell_designs/")
dir.create(designs_folder, showWarnings = F)


n_cores = parallel::detectCores()









##########################################################################################
#### Cocktail experiment
##########################################################################################

q = 3
J = 2
S = 16
n_draws_1 = 128

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

beta_correlated_draws_cocktail = get_correlated_halton_draws(beta0, sigma0, n_draws_1)


n_random_initial_starts_1 = 80
max_it_cocktail = 10
seed = 2020


cocktail_d_opt_filename = paste0(designs_folder, "cocktail_d_optimal.rds")
cocktail_i_opt_filename = paste0(designs_folder, "cocktail_i_optimal.rds")

##########################################################################################
#### Load or create the designs
##########################################################################################


if(file.exists(cocktail_d_opt_filename)){
  cat("D_B optimal design already exists.\n")
} else{
  # 15 mins with 64 random starts and 4 cores
  # 25 mins with 80 random starts and 4 cores
  cat("Doing D_B optimal design for cocktail experiment.\n")
  (t1D = Sys.time())
  cocktail_D_opt = mnl_mixture_coord_exch(
    n_random_starts = n_random_initial_starts_1,
    q = q,
    J = J,
    S = S,
    beta = beta_correlated_draws_cocktail,
    transform_beta = F,
    opt_method = "B",
    opt_crit = "D",
    max_it = max_it_cocktail,
    verbose = 1,
    plot_designs = F,
    seed = seed,
    n_cores = n_cores,
    save_all_designs = F
  )

  (t2D = Sys.time())
  t2D - t1D

  saveRDS(cocktail_D_opt, cocktail_d_opt_filename)
}



if(file.exists(cocktail_i_opt_filename)){
  cat("I_B optimal design already exists.\n")
} else{
  # 20 mins with 64 random starts and 4 cores
  cat("Doing I_B optimal design for cocktail experiment.\n")
  (t1I = Sys.time())
  cocktail_I_opt =  mnl_mixture_coord_exch(
    n_random_starts = n_random_initial_starts_1,
    q = q,
    J = J,
    S = S,
    beta = beta_correlated_draws_cocktail,
    transform_beta = F,
    opt_method = "B",
    opt_crit = "I",
    max_it = max_it_cocktail,
    verbose = 1,
    plot_designs = F,
    seed = seed,
    n_cores = n_cores,
    save_all_designs = F
  )
  (t2I = Sys.time())
  t2I - t1I

  saveRDS(cocktail_I_opt, cocktail_i_opt_filename)
}













##########################################################################################
#### Cornell's experiment
##########################################################################################

kappas = c(0.5, 5, 10, 30)
n_draws_2 = 128
n_random_initial_starts_2 = 80
max_it_cornell = 20


beta_2 = c(0.86, 0.21, 0, 3.07, 2.34, 3.24, -20.59)
beta_2_prime = c(0.86, 0.21, 3.07, 2.34, 3.24, -20.59)






###### Analytic transformation of betas


cornell_designs_basefilename_analytic_transf = paste0(designs_folder, "cornell_experiment_analytic_transformed_betas_maxit", max_it_cornell)


## Creating designs if they don't already exist
(start_time = Sys.time())

for(k in kappas){
  # Each design takes around 53 and 87 seconds with 16 initial random designs, 128 halton draws, 10 max iterations and 4 cores.
  # Each design takes around 99 and 160 seconds with 32 initial random designs, 128 halton draws, 10 max iterations and 4 cores.
  # Each design takes around 13 and 20 minutes with 128 initial random designs, 128 halton draws, 20 max iterations and 4 cores.

  cat("kappa =", k, "\n")


  Sigma_prime = transform_varcov_matrix(k*diag(7), 3)

  beta_2_prior_draws = get_correlated_halton_draws(beta_2_prime, Sigma_prime, n_draws_2)
  # beta_2_prior_draws = get_halton_draws(beta_2, sd = sqrt(k), ndraws = n_draws_2)

  d_opt_filename_analytic_transf = paste0(cornell_designs_basefilename_analytic_transf, "_kappa", k, "_Dopt.rds")

  if(file.exists(d_opt_filename_analytic_transf)){
    cat("\tD optimal file exists.\n")
  }else{
    cat("\tD optimal file does not exist. Creating design.\n")
    cornell_beta_2_analytic_transf_pseudo_bayesian_d_opt = mnl_mixture_coord_exch(
      q = 3,
      J = 2,
      S = 7,
      verbose = 1,
      n_random_starts = n_random_initial_starts_2,
      beta = beta_2_prior_draws,
      transform_beta = F,
      max_it = max_it_cornell,
      n_cores = n_cores,
      opt_crit = "D",
      plot_designs = F,
      save_all_designs = F)

    saveRDS(cornell_beta_2_analytic_transf_pseudo_bayesian_d_opt, d_opt_filename_analytic_transf)
  }


  i_opt_filename_analytic_transf = paste0(cornell_designs_basefilename_analytic_transf, "_kappa", k, "_Iopt.rds")

  if(file.exists(i_opt_filename_analytic_transf)){
    cat("\tI optimal file exists.\n")
  }else{
    cat("\tI optimal file does not exist. Creating design.\n")
    cornell_beta_2_analytic_transf_pseudo_bayesian_i_opt = mnl_mixture_coord_exch(
      q = 3,
      J = 2,
      S = 7,
      verbose = 1,
      n_random_starts = n_random_initial_starts_2,
      beta = beta_2_prior_draws,
      transform_beta = F,
      max_it = max_it_cornell,
      n_cores = n_cores,
      opt_crit = "I",
      plot_designs = F,
      save_all_designs = F)

    saveRDS(cornell_beta_2_analytic_transf_pseudo_bayesian_i_opt, i_opt_filename_analytic_transf)
  }

}










