library(testthat)
library(opdesmixr)


out_folder = here::here("tests/testthat/rds_mnl_moments_matrix/")



test_that("mnl_moments_matrix_no_pv",{
  for(order in 1:3){
    for(q in 3:15){
      base_filename = paste0(
        "mnl_moments_matrix_q", q, "_order", order, "_npv0"
      )

      mm = mnl_create_moment_matrix(q = q, n_pv = 0, order = order, pv_bounds = NULL)

      expect_equal_to_reference(
        object = mm,
        file = paste0(out_folder, base_filename, ".rds")
      )
    }
  }
})



test_that("mnl_moments_matrix_with_pv_pvboundsminus1to1",{
  for(q in c(3:10)){
    for(n_pv in 1:10){
      base_filename = paste0(
        "mnl_moments_matrix_q", q, "_npv", n_pv, "_order4", "_pvboundsminus1to1"
      )

      mm = mnl_create_moment_matrix(q = q, n_pv = n_pv, order = 4, pv_bounds = c(-1, 1))

      expect_equal_to_reference(
        object = mm,
        file = paste0(out_folder, base_filename, ".rds")
      )
    }
  }
})




test_that("mnl_moments_matrix_with_pv_pvbounds0to1",{
  for(q in c(3:10)){
    for(n_pv in 1:10){
      base_filename = paste0(
        "mnl_moments_matrix_q", q, "_npv", n_pv, "_order4", "_pvbounds0to1"
      )

      mm = mnl_create_moment_matrix(q = q, n_pv = n_pv, order = 4, pv_bounds = c(0, 1))

      expect_equal_to_reference(
        object = mm,
        file = paste0(out_folder, base_filename, ".rds")
      )
    }
  }
})








