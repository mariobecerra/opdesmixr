library(testthat)
library(opdesmixr)


# This test file checks that the final I-optimality and D-optimality values are equal to the previous versions.

out_folder = here::here("tests/testthat/rds_gaussian_designs/")


n_runs = 50
n_random_starts = 4
n_cores = parallel::detectCores()

test_that("gaussian_rds",{

  for(optimization_method in c("B", "D")){
    for(optimality_criterion in c("D", "I")){
      for(order in 1:3){
        for(q in c(3, 4, 5)){

          res_alg = gaussian_mixture_coord_exch(
            q = q,
            n_runs = n_runs,
            n_random_starts = n_random_starts,
            order = order,
            opt_crit = optimality_criterion,
            opt_method = optimization_method,
            n_cox_points = 50,
            plot_designs = F,
            verbose = 0,
            n_cores = n_cores,
            max_it = 3,
            seed = 10)

          base_filename = paste0(
            "gaussian_", optimality_criterion, "opt_q", q, "_nruns", n_runs, "_rs", n_random_starts, "_order", order, "_optmeth", optimization_method
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


        }

      }
    }
  }


})










