library(tidyverse)
library(here)

devtools::load_all(".")

ggplot2::theme_set(ggplot2::theme_bw())

# Rcpp::sourceCpp("src/utils.cpp")
# source("R/mnl_model_functions.R")
# source("R/general_functions.R")


out_dir_1 = "tests/i_optimality_vs_d_optimality/out/"
out_dir_2 = paste0(out_dir_1, "rds_files/")
out_dir = paste0(out_dir_2, "many_random_designs/")

file_names = list.files(here(out_dir))


efficiencies = map_df(seq_along(file_names), function(i){

  filename_i = file_names[i]

  q = str_extract(pattern = "[0-9]+", str_extract(pattern = "_.*?_", string = filename_i))
  J = str_extract(pattern = "[0-9]+", str_extract(pattern = "J.*?_", string = filename_i))
  S = str_extract(pattern = "[0-9]+", str_extract(pattern = "S.*?_", string = filename_i))

  designs = readRDS(here(out_dir, filename_i))
  D_opt = designs$D_opt
  I_opt = designs$I_opt
  beta = designs$beta

  d_eff_of_d_optimal_design = try(get_opt_crit_value_MNL(D_opt$X, beta = beta, opt_crit = 0), silent = T)
  d_eff_of_i_optimal_design = try(get_opt_crit_value_MNL(I_opt$X, beta = beta, opt_crit = 0), silent = T)
  i_eff_of_d_optimal_design = try(get_opt_crit_value_MNL(D_opt$X, beta = beta, opt_crit = 1), silent = T)
  i_eff_of_i_optimal_design = try(get_opt_crit_value_MNL(I_opt$X, beta = beta, opt_crit = 1), silent = T)


  if(
    class(d_eff_of_d_optimal_design) == 'try-error' |
    class(d_eff_of_i_optimal_design) == 'try-error' |
    class(i_eff_of_d_optimal_design) == 'try-error' |
    class(i_eff_of_i_optimal_design) == 'try-error'
  ){
    d_eff_of_d_optimal_design = NA_real_
    d_eff_of_i_optimal_design = NA_real_
    i_eff_of_d_optimal_design = NA_real_
    i_eff_of_i_optimal_design = NA_real_
  }

  out = tibble(
    d_eff_of_d_optimal_design = d_eff_of_d_optimal_design,
    d_eff_of_i_optimal_design = d_eff_of_i_optimal_design,
    i_eff_of_d_optimal_design = i_eff_of_d_optimal_design,
    i_eff_of_i_optimal_design = i_eff_of_i_optimal_design
  ) %>%
    mutate(q = as.numeric(q), J = J, S = S)

  return(out)

}) %>%
  mutate(n_params = (q^3 + 5*q)/6-1) %>%
  mutate(
    ratio_d_eff = (exp(d_eff_of_i_optimal_design - d_eff_of_d_optimal_design))^(1/n_params),
    ratio_i_eff = (exp(i_eff_of_d_optimal_design - i_eff_of_i_optimal_design))
  )


saveRDS(efficiencies, here(out_dir_2, "efficiencies_all.rds"))

