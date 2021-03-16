library(testthat)
library(opdesmixr)
library(here)

# The initial designs in the folder defined by the 'design_and_vector_folder' variable were created programmatically with an external script.

design_and_vector_folder = here("tests/testthat/rds_mnl_initial_designs_do_not_edit_by_hand/")
out_folder = here("tests/testthat/rds_mnl_locally_optimal_designs_initial_designs_output/")
dir.create(out_folder)

seed = 10



test_that("mnl_locally_optimal_with_initial_designs_beta_0",{

  for(q in c(3:5)){
    J = 2
    S = 10*q

    design_filename = paste0("mnl_initial_design_q", q, "_J", J, "_S", S, ".rds")
    initial_design = readRDS(paste0(design_and_vector_folder, design_filename))

    for(order in 1:3){

      m = length(mnl_create_random_beta(q = q, order = order)$beta)
      beta_vec = rep(0, m)

      for(optimization_method in c("B", "D")){
        for(optimality_criterion in c("D", "I")){

          if(optimization_method == "D") n_cox_points = 30
          else n_cox_points = NULL

          res_alg = mnl_mixture_coord_exch(
            q = q,
            J = J,
            S = S,
            beta = beta_vec,
            X = initial_design,
            order = order,
            opt_crit = optimality_criterion,
            opt_method = optimization_method,
            n_cox_points = n_cox_points,
            plot_designs = F,
            verbose = 0,
            n_cores = 1,
            max_it = 3,
            seed = seed)

          base_filename = paste0(
            "mnl_", optimality_criterion, "opt_q", q, "_J", J, "_S", S, "_order", order, "_optmeth", optimization_method, "_beta0"
            )


          expect_equal_to_reference(
            object = res_alg$X_orig,
            file = paste0(out_folder, base_filename, "_X_orig.rds")
          )

          expect_equal_to_reference(
            object = res_alg$X,
            file = paste0(out_folder, base_filename, "_X.rds")
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
  }

})






test_that("mnl_locally_optimal_with_initial_designs_beta_1",{

  for(q in c(3:5)){
    J = 2
    S = 10*q

    design_filename = paste0("mnl_initial_design_q", q, "_J", J, "_S", S, ".rds")
    initial_design = readRDS(paste0(design_and_vector_folder, design_filename))

    for(order in 1:3){

      m = length(mnl_create_random_beta(q = q, order = order)$beta)
      beta_vec = rep(1, m)

      for(optimization_method in c("B", "D")){
        for(optimality_criterion in c("D", "I")){

          if(optimization_method == "D") n_cox_points = 30
          else n_cox_points = NULL

          res_alg = mnl_mixture_coord_exch(
            q = q,
            J = J,
            S = S,
            beta = beta_vec,
            X = initial_design,
            order = order,
            opt_crit = optimality_criterion,
            opt_method = optimization_method,
            n_cox_points = n_cox_points,
            plot_designs = F,
            verbose = 0,
            n_cores = 1,
            max_it = 3,
            seed = seed)

          base_filename = paste0(
            "mnl_", optimality_criterion, "opt_q", q, "_J", J, "_S", S, "_order", order, "_optmeth", optimization_method, "_beta1"
          )


          expect_equal_to_reference(
            object = res_alg$X_orig,
            file = paste0(out_folder, base_filename, "_X_orig.rds")
          )

          expect_equal_to_reference(
            object = res_alg$X,
            file = paste0(out_folder, base_filename, "_X.rds")
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
  }

})





test_that("mnl_locally_optimal_with_initial_designs_random_beta",{

  for(q in c(3:5)){
    J = 2
    S = 10*q

    design_filename = paste0("mnl_initial_design_q", q, "_J", J, "_S", S, ".rds")
    initial_design = readRDS(paste0(design_and_vector_folder, design_filename))

    for(order in 1:3){

      vector_filename = paste0("mnl_beta_q", q, "_J", J, "_S", S, "_order", order, ".rds")
      beta_vec = readRDS(paste0(design_and_vector_folder, vector_filename))

      for(optimization_method in c("B", "D")){
        for(optimality_criterion in c("D", "I")){

          if(optimization_method == "D") n_cox_points = 30
          else n_cox_points = NULL

          res_alg = mnl_mixture_coord_exch(
            q = q,
            J = J,
            S = S,
            beta = beta_vec,
            X = initial_design,
            order = order,
            opt_crit = optimality_criterion,
            opt_method = optimization_method,
            n_cox_points = n_cox_points,
            plot_designs = F,
            verbose = 0,
            n_cores = 1,
            max_it = 3,
            seed = seed)

          base_filename = paste0(
            "mnl_", optimality_criterion, "opt_q", q, "_J", J, "_S", S, "_order", order, "_optmeth", optimization_method, "_randombeta"
          )


          expect_equal_to_reference(
            object = res_alg$X_orig,
            file = paste0(out_folder, base_filename, "_X_orig.rds")
          )

          expect_equal_to_reference(
            object = res_alg$X,
            file = paste0(out_folder, base_filename, "_X.rds")
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
  }

})
