library(testthat)
library(opdesmixr)

out_folder = here::here("tests/testthat/rds_mnl_pseudo_bayesian_optimal_designs/")

seed = 10
n_random_starts = 2
n_cox_points = 10
n_cores = parallel::detectCores()

n_draws = 32
variance = 0.5



get_size = function(q, order){

  m1 = q
  m2 = q*(q-1)/2
  m3 = q*(q-1)*(q-2)/6
  if(order == 1){
    m = m1
  } else{
    if(order == 2){
      m = m1 + m2
    } else{
      m = m1 + m2 + m3 # = q + q*(q-1)/2 + q*(q-1)*(q-2)/6
    }
  }
  return(m)
}




test_that("rds_mnl_pseudo_bayesian_optimal_designs",{

  for(optimization_method in c("B", "D")){
    for(optimality_criterion in c("D", "I")){
      for(order in 1:3){
        for(q in c(4, 5)){
          J = 2
          S = 5*q

          m = get_size(q, order)

          beta_vec = rep(1, m)
          beta_prior_draws = get_halton_draws(beta_vec, sd = sqrt(variance), ndraws = n_draws)


            res_alg = suppressWarnings(
              mnl_mixture_coord_exch(
              q = q,
              J = J,
              S = S,
              beta = beta_prior_draws,
              n_random_starts = n_random_starts,
              order = order,
              opt_crit = optimality_criterion,
              opt_method = optimization_method,
              plot_designs = F,
              verbose = 0,
              n_cores = n_cores,
              max_it = 3,
              n_cox_points = n_cox_points,
              seed = seed)
          )

          base_filename = paste0(
            "mnl_", optimality_criterion, "_q", q, "_J", J, "_S", S, "_rs", n_random_starts, "_order", order, "_optmeth", optimization_method
          )

          # expect_equal_to_reference(
          #   object = res_alg$X_orig,
          #   file = paste0(out_folder, base_filename, "_X_orig.rds")
          # )

          # expect_equal_to_reference(
          #   object = res_alg$X,
          #   file = paste0(out_folder, base_filename, "_X.rds")
          # )

          expect_equal_to_reference(
            object = res_alg$opt_crit_value_orig,
            file = paste0(out_folder, base_filename, "_opt_crit_value_orig.rds")
          )

          expect_equal_to_reference(
            object = res_alg$opt_crit_value,
            file = paste0(out_folder, base_filename, "_opt_crit_value_final.rds")
          )

          expect_equal_to_reference(
            object = res_alg$beta,
            file = paste0(out_folder, base_filename, "_beta.rds")
          )

        }

      }
    }
  }


})
