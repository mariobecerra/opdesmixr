library(opdesmixr)
library(here)

out_folder = "tests/testthat/rds_mnl_initial_designs_do_not_edit_by_hand/"
dir.create(here(out_folder))

seed = 10
n_draws = 32
variance = 0.5

for(q in c(3:5)){
  J = 2
  S = 10*q

  initial_design = mnl_create_random_initial_design(q = q, J = J, S = S, seed = seed)
  design_filename = paste0("mnl_initial_design_q", q, "_J", J, "_S", S, ".rds")
  saveRDS(initial_design, here(out_folder, design_filename))

  for(order in 1:3){

    beta_vec = mnl_create_random_beta(q, order = order, seed = seed)$beta
    beta_prior_draws = get_halton_draws(beta_vec, sd = sqrt(variance), ndraws = n_draws)

    vector_filename = paste0("mnl_beta_q", q, "_J", J, "_S", S, "_order", order, ".rds")
    matrix_filename = paste0("mnl_beta_draws_q", q, "_J", J, "_S", S, "_order", order, ".rds")

    saveRDS(beta_vec, here(out_folder, vector_filename))
    saveRDS(beta_prior_draws, here(out_folder, matrix_filename))
  }
}

