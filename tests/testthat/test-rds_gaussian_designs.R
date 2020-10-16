library(testthat)
library(opdesmixr)


# This test file checks that the final design matrices are equal to the previous versions.
# I only check the X matrices because I wrote the tests by hand.
# Should have done it programatically but I didn't know how to do it at the time.


out_folder = here::here("tests/testthat/rds_gaussian_designs/")

# D-optimality, 30 runs, 3 ingredients

test_that("gaussian_D_q3_nruns30_rs100_o1_seed10", {
  expect_equal_to_reference(
    object = gaussian_mixture_coord_exch(
      n_runs = 30,
      q = 3,
      n_random_starts = 100,
      order = 1,
      opt_crit = "D",
      opt_method = "B",
      plot_designs = F,
      verbose = 0,
      n_cores = 1,
      max_it = 5,
      seed = 10)$X,
    file = paste0(out_folder, "gaussian_D_q3_nruns30_rs100_o1_seed10_brent_X.rds")
  )

  expect_equal_to_reference(
    object = gaussian_mixture_coord_exch(
      n_runs = 30,
      q = 3,
      n_random_starts = 100,
      order = 1,
      opt_crit = "D",
      opt_method = "D",
      n_cox_points = 50,
      plot_designs = F,
      verbose = 0,
      n_cores = 1,
      max_it = 5,
      seed = 10)$X,
    file = paste0(out_folder, "gaussian_D_q3_nruns30_rs100_o1_seed10_discrete_X.rds"))
})


test_that("gaussian_D_q3_nruns30_rs100_o2_seed10", {
  expect_equal_to_reference(
    object = gaussian_mixture_coord_exch(
      n_runs = 30,
      q = 3,
      n_random_starts = 100,
      order = 2,
      opt_crit = "D",
      opt_method = "B",
      plot_designs = F,
      verbose = 0,
      n_cores = 1,
      max_it = 5,
      seed = 10)$X,
    file = paste0(out_folder, "gaussian_D_q3_nruns30_rs100_o2_seed10_brent_X.rds"))

  expect_equal_to_reference(
    object = gaussian_mixture_coord_exch(
      n_runs = 30,
      q = 3,
      n_random_starts = 100,
      order = 2,
      opt_crit = "D",
      opt_method = "D",
      n_cox_points = 50,
      plot_designs = F,
      verbose = 0,
      n_cores = 1,
      max_it = 5,
      seed = 10)$X,
    file = paste0(out_folder, "gaussian_D_q3_nruns30_rs100_o2_seed10_discrete_X.rds"))
})



test_that("gaussian_D_q3_nruns30_rs50_o3_seed10", {
  expect_equal_to_reference(
    object = gaussian_mixture_coord_exch(
      n_runs = 30,
      q = 3,
      n_random_starts = 50,
      order = 2,
      opt_crit = "D",
      opt_method = "B",
      plot_designs = F,
      verbose = 0,
      n_cores = 1,
      max_it = 5,
      seed = 10)$X,
    file = paste0(out_folder, "gaussian_D_q3_nruns30_rs50_o3_seed10_brent_X.rds")
  )

  expect_equal_to_reference(
    object = gaussian_mixture_coord_exch(
      n_runs = 30,
      q = 3,
      n_random_starts = 50,
      order = 2,
      opt_crit = "D",
      opt_method = "D",
      n_cox_points = 50,
      plot_designs = F,
      verbose = 0,
      n_cores = 1,
      max_it = 5,
      seed = 10)$X,
    file = paste0(out_folder, "gaussian_D_q3_nruns30_rs50_o3_seed10_discrete_X.rds")
  )
})





# I-optimality, 30 runs, 3 ingredients

test_that("gaussian_I_q3_nruns30_rs50_o3_seed10", {
  expect_equal_to_reference(
    object = gaussian_mixture_coord_exch(
      n_runs = 30,
      q = 3,
      n_random_starts = 50,
      order = 3,
      opt_crit = "I",
      opt_method = "B",
      plot_designs = F,
      verbose = 0,
      n_cores = 1,
      max_it = 5,
      seed = 10)$X,
    file = paste0(out_folder, "gaussian_I_q3_nruns30_rs50_o3_seed10_brent_X.rds")
  )

  expect_equal_to_reference(
    object = gaussian_mixture_coord_exch(
      n_runs = 30,
      q = 3,
      n_random_starts = 50,
      order = 3,
      opt_crit = "I",
      opt_method = "D",
      n_cox_points = 30,
      plot_designs = F,
      verbose = 0,
      n_cores = 1,
      max_it = 5,
      seed = 10)$X,
    file = paste0(out_folder, "gaussian_I_q3_nruns30_rs50_o3_seed10_discrete_X.rds")
  )
})



# D-optimality, 60 runs, 5 ingredients

test_that("gaussian_D_q5_nruns60_rs50_o1_seed10", {
  expect_equal_to_reference(
    object = gaussian_mixture_coord_exch(
      n_runs = 60,
      q = 5,
      n_random_starts = 50,
      order = 1,
      opt_crit = "D",
      opt_method = "B",
      plot_designs = F,
      verbose = 0,
      n_cores = 1,
      max_it = 5,
      seed = 10)$X,
    file = paste0(out_folder, "gaussian_D_q5_nruns60_rs50_o1_seed10_brent_X.rds")
  )

  expect_equal_to_reference(
    object = gaussian_mixture_coord_exch(
      n_runs = 60,
      q = 5,
      n_random_starts = 50,
      order = 1,
      opt_crit = "D",
      opt_method = "D",
      n_cox_points = 30,
      plot_designs = F,
      verbose = 0,
      n_cores = 1,
      max_it = 5,
      seed = 10)$X,
    file = paste0(out_folder, "gaussian_D_q5_nruns60_rs50_o1_seed10_discrete_X.rds")
  )
})



test_that("gaussian_D_q5_nruns60_o2_seed10", {
  expect_equal_to_reference(
    object = gaussian_mixture_coord_exch(
      n_runs = 60,
      q = 5,
      n_random_starts = 10,
      order = 2,
      opt_crit = "D",
      opt_method = "B",
      plot_designs = F,
      verbose = 0,
      n_cores = 1,
      max_it = 3,
      seed = 10)$X,
    file = paste0(out_folder, "gaussian_D_q5_nruns60_o2_seed10_brent_X.rds")
  )

  expect_equal_to_reference(
    object = gaussian_mixture_coord_exch(
      n_runs = 60,
      q = 5,
      n_random_starts = 10,
      order = 2,
      opt_crit = "D",
      opt_method = "D",
      n_cox_points = 30,
      plot_designs = F,
      verbose = 0,
      n_cores = 1,
      max_it = 3,
      seed = 10)$X,
    file = paste0(out_folder, "gaussian_D_q5_nruns60_o2_seed10_discrete_X.rds")
  )
})


test_that("gaussian_D_q5_nruns60_o3_seed10", {
  expect_equal_to_reference(
    object = gaussian_mixture_coord_exch(
      n_runs = 60,
      q = 5,
      n_random_starts = 10,
      order = 3,
      opt_crit = "D",
      opt_method = "B",
      plot_designs = F,
      verbose = 0,
      n_cores = 1,
      max_it = 3,
      seed = 10)$X,
    file = paste0(out_folder, "gaussian_D_q5_nruns60_o3_seed10_brent_X.rds")
  )

  expect_equal_to_reference(
    object = gaussian_mixture_coord_exch(
      n_runs = 60,
      q = 5,
      n_random_starts = 10,
      order = 3,
      opt_crit = "D",
      opt_method = "D",
      n_cox_points = 30,
      plot_designs = F,
      verbose = 0,
      n_cores = 1,
      max_it = 3,
      seed = 10)$X,
    file = paste0(out_folder, "gaussian_D_q5_nruns60_o3_seed10_discrete_X.rds")
  )
})


# I-optimality, 60 runs, 5 ingredients

test_that("gaussian_I_q5_nruns60_rs50_o3_seed10", {
  expect_equal_to_reference(
    object = gaussian_mixture_coord_exch(
      n_runs = 60,
      q = 5,
      n_random_starts = 10,
      order = 3,
      opt_crit = "I",
      opt_method = "B",
      plot_designs = F,
      verbose = 0,
      n_cores = 1,
      max_it = 2,
      seed = 10)$X,
    file = paste0(out_folder, "gaussian_I_q5_nruns60_rs50_o3_seed10_brent_X.rds")
  )

  expect_equal_to_reference(
    object = gaussian_mixture_coord_exch(
      n_runs = 60,
      q = 5,
      n_random_starts = 10,
      order = 3,
      opt_crit = "I",
      opt_method = "D",
      n_cox_points = 30,
      plot_designs = F,
      verbose = 0,
      n_cores = 1,
      max_it = 2,
      seed = 10)$X,
    file = paste0(out_folder, "gaussian_I_q5_nruns60_rs50_o3_seed10_discrete_X.rds")
  )
})





