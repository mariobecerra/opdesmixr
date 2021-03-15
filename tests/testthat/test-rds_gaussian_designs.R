library(testthat)
library(opdesmixr)


# This test file checks that the final I-optimality and D-optimality values are equal to the previous versions.

out_folder = here::here("tests/testthat/rds_gaussian_designs/")


n_runs = 50
n_random_starts = 4
n_cores = parallel::detectCores()

test_that("gaussian_rds_no_pv",{

  for(optimization_method in c("B", "D")){
    for(optimality_criterion in c("D", "I")){
      for(order in 1:3){
        for(q in c(3, 4, 5)){

          res_alg = gaussian_mixture_coord_exch(
            q = q,
            n_pv = 0,
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







test_that("gaussian_rds_with_pv_pvboundsminus1to1",{

  optimization_method = "B"
  for(optimality_criterion in c("D", "I")){
    for(q in c(3, 4, 5)){
      for(n_pv in 1:5){

        res_alg = gaussian_mixture_coord_exch(
          q = q,
          n_pv = n_pv,
          n_runs = n_runs,
          n_random_starts = n_random_starts,
          order = 4,
          pv_bounds = c(-1, 1),
          opt_crit = optimality_criterion,
          opt_method = optimization_method,
          n_cox_points = NULL,
          plot_designs = F,
          verbose = 0,
          n_cores = n_cores,
          max_it = 3,
          seed = 10)

        base_filename = paste0(
          "gaussian_", optimality_criterion, "opt_q", q, "_npv", n_pv, "_nruns", n_runs, "_rs", n_random_starts, "_order4_optmeth", optimization_method, "_pvboundsminus1to1"
        )

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


})










test_that("gaussian_rds_with_pv_pvbounds0to1",{

  optimization_method = "B"
  for(optimality_criterion in c("D", "I")){
    for(q in c(3, 4, 5)){
      for(n_pv in 1:5){

        res_alg = gaussian_mixture_coord_exch(
          q = q,
          n_pv = n_pv,
          n_runs = n_runs,
          n_random_starts = n_random_starts,
          order = 4,
          pv_bounds = c(0, 1),
          opt_crit = optimality_criterion,
          opt_method = optimization_method,
          n_cox_points = NULL,
          plot_designs = F,
          verbose = 0,
          n_cores = n_cores,
          max_it = 3,
          seed = 10)

        base_filename = paste0(
          "gaussian_", optimality_criterion, "opt_q", q, "_npv", n_pv, "_nruns", n_runs, "_rs", n_random_starts, "_order4_optmeth", optimization_method, "_pvbounds0to1"
        )

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


})









