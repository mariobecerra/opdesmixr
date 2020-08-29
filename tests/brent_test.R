library(tidyverse)
# library(opdesmixr)
devtools::load_all(".")

X = opdesmixr::create_random_initial_design_gaussian(10, 3, seed = 4)
get_opt_crit_value_Gaussian(X, order = 1, opt_crit = 0)
efficiency_cox_scheffe_gaussian(theta = 0.65932440, X = X, j = 0, i = 0, order = 1, opt_crit = 0)
efficiency_cox_scheffe_gaussian(theta = 1, X = X, j = 0, i = 0, order = 1, opt_crit = 0)


X2 = opdesmixr::create_random_initial_design_gaussian(30, 3, seed = 4)
get_opt_crit_value_Gaussian(X2, order = 3, opt_crit = 1)
efficiency_cox_scheffe_gaussian(theta = 0.65932440, X = X2, j = 0, i = 0, order = 3, opt_crit = 1)
efficiency_cox_scheffe_gaussian(theta = 0.1, X = X2, j = 0, i = 0, order = 3, opt_crit = 1)
efficiency_cox_scheffe_gaussian(theta = 1, X = X2, j = 0, i = 0, order = 3, opt_crit = 1)


# We want to minimize this
tibble(theta = 0:50/50) %>%
  mutate(u = map_dbl(theta, function(x) efficiency_cox_scheffe_gaussian(theta = x, X = X2, j = 0, i = 0, order = 3, opt_crit = 1))) %>%
  ggplot() +
  geom_line(aes(theta, u))



# We want to minimize this
tibble(theta = 0:50/50) %>%
  mutate(u = map_dbl(theta, function(x) efficiency_cox_scheffe_gaussian(theta = x, X = X2, j = 0, i = 0, order = 3, opt_crit = 0))) %>%
  ggplot() +
  geom_line(aes(theta, u))





fn_optimize = function(X, j, i, order, opt_crit){
  out_fn = function(theta){
    # efficiency_cox_scheffe_gaussian(theta = x, X = X2, j = 0, i = 0, order = 3, opt_crit = 1)
    return(efficiency_cox_scheffe_gaussian(theta = theta, X = X, j = j, i = i, order = order, opt_crit = opt_crit))
  }
  return(out_fn)
}


f_i = fn_optimize(X = X2, j = 0, i = 0, order = 3, opt_crit = 0)

optim(0.5, f_i, method = "Brent", lower = 0, upper = 1)
optimize(f_i, c(0, 1))
optim(0, f_i, method = "Brent", lower = 0, upper = 1)
BrentCoxScheffeGaussian(X = X2, j = 0, i = 0, order = 3, opt_crit = 0, W = matrix(0.0, nrow = 1))

brent_cox_scheffe_gaussian(
  X = X2, j = 0, i = 0, order = 3, opt_crit = 0,
  lower = 0, upper = 1, tol = 0.0001)

brent_global_cox_scheffe_gaussian(
  X = X2, j = 0, i = 0, order = 3, opt_crit = 0,
  lower = 0,
  upper = 1,
  initial_guess = 0.5,
  hessian_bound = 1e5,
  abs_err_tol = 0.0001,
  tol = 0.0001)




microbenchmark::microbenchmark(
  brent_cox_scheffe_gaussian(
    X = X2, j = 0, i = 0, order = 3, opt_crit = 0,
    lower = 0, upper = 1, tol = 0.001)

  ,

  brent_global_cox_scheffe_gaussian(
    X = X2, j = 0, i = 0, order = 3, opt_crit = 0,
    lower = 0,
    upper = 1,
    initial_guess = 0.5,
    hessian_bound = 100,
    abs_err_tol = 0.001,
    tol = 0.0001)

  ,


  unit = "relative"
)









# Third degree
# res_alg_order_3_2 = mixture_coord_ex(
#   X = X2,
#   order = 3,
#   model = "Gaussian",
#   plot_designs = F,
#   n_cox_points = 100,
#   max_it = 2)

# res_alg_order_3_2 = mixture_coord_ex_gaussian(
#   X = X2,
#   order = 3,
#   opt_crit = 0,
#   plot_designs = F,
#   n_cox_points = 100,
#   max_it = 2,
#   verbose = 0,
#   opt_method = 1)

res_alg_order_3_2 = mixtureCoordinateExchangeGaussian(
  X_orig = X2,
  order = 3, n_cox_points = 100, max_it = 2, verbose = 0, opt_crit = 0,
  W = matrix(0.0, nrow = 1),
  opt_method = 1, lower = 0, upper = 1, tol = 0.0001)

X3 = res_alg_order_3_2$X


res_alg_order_3_3 = mixtureCoordinateExchangeGaussian(
  X_orig = X2,
  order = 3, n_cox_points = 100, max_it = 2, verbose = 0, opt_crit = 0,
  W = matrix(0.0, nrow = 1),
  opt_method = 0, lower = 0, upper = 1, tol = 0.0001)

X3_2 = res_alg_order_3_3$X





tibble(theta = 0:50/50) %>%
  mutate(u = map_dbl(theta, function(x) efficiency_cox_scheffe_gaussian(theta = x, X = X3, j = 0, i = 0, order = 3, opt_crit = 0))) %>%
  ggplot() +
  geom_line(aes(theta, u))



brent_cox_scheffe_gaussian(
  X = X3, j = 0, i = 0, order = 3, opt_crit = 0,
  lower = 0, upper = 1, tol = 0.0001)

brent_global_cox_scheffe_gaussian(
  X = X3, j = 0, i = 0, order = 3, opt_crit = 0,
  lower = 0,
  upper = 1,
  initial_guess = 0.5,
  hessian_bound = 1e5,
  abs_err_tol = 0.0001,
  tol = 0.0001)






