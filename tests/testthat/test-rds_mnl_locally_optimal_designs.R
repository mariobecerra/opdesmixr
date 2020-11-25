library(testthat)
library(opdesmixr)

out_folder = here::here("tests/testthat/rds_mnl_locally_optimal_designs/")

seed = 10
n_random_starts = 4
n_cores = parallel::detectCores()


test_that("mnl_locally_optimal",{

  for(optimization_method in c("B", "D")){
    for(optimality_criterion in c("D", "I")){
      for(order in 1:3){
        for(q in c(3, 4, 5)){
          J = 2
          S = 5*q

          beta_vec = mnl_create_random_beta(q, order = order, seed = seed)$beta

          res_alg = mnl_mixture_coord_exch(
            q = q,
            J = J,
            S = S,
            beta = beta_vec,
            n_random_starts = n_random_starts,
            order = order,
            opt_crit = optimality_criterion,
            opt_method = optimization_method,
            plot_designs = F,
            verbose = 0,
            n_cores = n_cores,
            max_it = 3,
            seed = seed)

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
