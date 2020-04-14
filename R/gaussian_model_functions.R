
#' TODO: write doc
#' @export
create_random_initial_design_gaussian = function(n_runs, q, seed = NULL){
  X = matrix(rep(NA_real_, n_runs*q), nrow = n_runs)

  if(!is.null(seed)) set.seed(seed)

  for(i in 1:nrow(X)){
    rands = runif(q)

    # rows in X sum to 1:
    X[i,] = rands/sum(rands)
  }

  return(X)
}




#' TODO: write doc
#' @export
mixture_coord_ex_gaussian = function(
  X,
  order = 1,
  n_cox_points = 100,
  max_it = 50,
  plot_designs = F,
  verbose = 1,
  opt_crit = 0){



  n_runs = nrow(X)
  q = ncol(X)

  # Check that rows in X sum to 1
  row_sums = apply(X, 1, sum)

  if(sum(abs(row_sums - rep(1, n_runs)) < 1e-10) != n_runs){
    stop("Rows in X must sum 1")
  }


  # Coordinate exchanges:
  X_result = mixtureCoordinateExchangeGaussian(X, order, n_cox_points, max_it, verbose, opt_crit)

  out_list = list(
    X_orig = X_result$X_orig,
    X = X_result$X,
    d_eff_orig = X_result$d_eff_orig,
    d_eff = X_result$d_eff,
    n_iter = X_result$n_iter
  )

  if(plot_designs) {
    if(q == 3) gaussian_plot_result(out_list)
    else warning("Could not plot results because q != 3")
  }
  return(out_list)

}




#' TODO: write doc
#' @export
gaussian_plot_result = function(res_alg){
  # res_alg: output of a call to mixture_coord_ex_gaussian() function

  ggtern::grid.arrange(
    res_alg$X_orig %>%
      dplyr::as_tibble() %>%
      purrr::set_names(c("c1", "c2", "c3")) %>%
      ggtern(aes(c1, c2, c3)) +
      geom_point(shape = "x", size = 4) +
      theme_minimal() +
      ggtitle(label = "",
              subtitle = paste0("log D-efficiency = ", round(res_alg$d_eff_orig, 3)))
    ,
    res_alg$X %>%
      dplyr::as_tibble() %>%
      purrr::set_names(c("c1", "c2", "c3")) %>%
      ggtern(aes(c1, c2, c3)) +
      geom_point(shape = "x", size = 4) +
      theme_minimal() +
      ggtitle(label = "",
              subtitle = paste0("log D-efficiency = ", round(res_alg$d_eff, 3)))
    ,
    ncol = 2
  )


}

