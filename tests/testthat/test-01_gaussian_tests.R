library(testthat)
library(digest)
library(opdesmixr)

# D-optimality, 30 runs, 3 ingredients

test_that("gaussian_D_q3_nruns30_rs100_o1_seed10", {
  expect_identical(
    digest::digest(
      gaussian_mixture_coord_exch(
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
      seed = 10)
      ),
    "1459dc65acc87933d1829ac4c4cd02e8"
  )

  expect_identical(
    digest::digest(
      gaussian_mixture_coord_exch(
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
        seed = 10)
    ),
    "f3ec31fe144b6576702ccda0277d51f1"
  )
})


test_that("gaussian_D_q3_nruns30_rs100_o2_seed10", {
  expect_identical(
    digest::digest(
      gaussian_mixture_coord_exch(
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
        seed = 10)
    ),
    "b7b4d8f9d995a17239199ef4d7952364"
  )

  expect_identical(
    digest::digest(
      gaussian_mixture_coord_exch(
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
        seed = 10)
    ),
    "049c48c0e71d58a97e9a8ab7c70cfb8a"
  )
})



test_that("gaussian_D_q3_nruns30_rs50_o3_seed10", {
  expect_identical(
    digest::digest(
      gaussian_mixture_coord_exch(
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
        seed = 10)
    ),
    "9783f6eabd99123dc5ff04117464f759"
  )

  expect_identical(
    digest::digest(
      gaussian_mixture_coord_exch(
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
        seed = 10)
    ),
    "049c48c0e71d58a97e9a8ab7c70cfb8a"
  )
})





# I-optimality, 30 runs, 3 ingredients

test_that("gaussian_I_q3_nruns30_rs50_o3_seed10", {
  expect_identical(
    digest::digest(
      gaussian_mixture_coord_exch(
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
        seed = 10)
    ),
    "814a95707267c69373f3422a49f45159"
  )

  expect_identical(
    digest::digest(
      gaussian_mixture_coord_exch(
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
        seed = 10)
    ),
    "aa19ab2a4f1cafeba8ed0c57d7cc4691"
  )
})



# D-optimality, 60 runs, 5 ingredients

test_that("gaussian_D_q5_nruns60_rs50_o1_seed10", {
  expect_identical(
    digest::digest(
      gaussian_mixture_coord_exch(
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
        seed = 10)
    ),
    "85a7cfe0b422c0ca300ccb928f1ef7e9"
  )

  expect_identical(
    digest::digest(
      gaussian_mixture_coord_exch(
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
        seed = 10)
    ),
    "4b2d39f62bd53d77ad2d41f826728b70"
  )
})



test_that("gaussian_D_q5_nruns60_rs10_o2_seed10", {
  expect_identical(
    digest::digest(
      gaussian_mixture_coord_exch(
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
        seed = 10)
    ),
    "dea7b9d9e3297539399c6fe9ad2b1caa"
  )

  expect_identical(
    digest::digest(
      gaussian_mixture_coord_exch(
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
        seed = 10)
    ),
    "d9ac19cd257f457f7df7356553314f89"
  )
})


test_that("gaussian_D_q5_nruns60_rs10_o3_seed10", {
  expect_identical(
    digest::digest(
      gaussian_mixture_coord_exch(
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
        seed = 10)
    ),
    "48d78f55c907886fe3baef1eca27b277"
  )

  expect_identical(
    digest::digest(
      gaussian_mixture_coord_exch(
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
        seed = 10)
    ),
    "c01b20f95541261f825ce76ae4736979"
  )

})


# I-optimality, 60 runs, 5 ingredients

test_that("gaussian_I_q5_nruns60_rs50_o3_seed10", {
  expect_identical(
    digest::digest(
      gaussian_mixture_coord_exch(
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
        seed = 10)
    ),
    "f9f391607cc21fd1844bac607afe4af2"
  )

  expect_identical(
    digest::digest(
      gaussian_mixture_coord_exch(
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
        seed = 10)
    ),
    "834388be0ba978078aedeb1143a3794e"
  )

})





