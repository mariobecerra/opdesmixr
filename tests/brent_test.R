library(tidyverse)
# library(opdesmixr)
devtools::load_all(".")


efficiency_cox_scheffe_gaussian = function(theta, X, j, i, order, opt_crit){
  # // Computes efficiency criterion of a design matrix X but where the j-th ingredient in i-th observation is changed to theta.
  # // theta must be between 0 and 1 because it's an ingredient proportion.
  # // j and i are 0-indexed.
  # // We want to minimize this.
  #
  # // Create new matrix Y that is identical to the one pointed by X.
  # // Note: This is the easiest way to do it because we have to modify a row in this matrix.
  # // A more computationally effective way would be to only store the new modified vector since
  # // we don't need a copy of the whole matrix. But to do that I would have to either modify some
  # // existing functions, or create some new ones, or both. IDK if the gain in performance is worth it.


  q = dim(X)[2]

  if(opt_crit == 0){
    # "D-optimality"
    W = matrix(0.0, nrow = 1)
  } else{
    # "I-optimality")
    W = create_moment_matrix_gaussian(q)
  }

  return(efficiencyCoxScheffeGaussian(theta, X, j, i, order, opt_crit, W))

}

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











mixture_coord_ex_gaussian(
  n_runs = 30,
  q = 3,
  n_random_starts = 100,
  X = NULL,
  order = 3,
  opt_method = "B",
  max_it = 10,
  tol = 0.0001,
  n_cox_points = NULL,
  plot_designs = T,
  verbose = 0,
  opt_crit = "D",
  seed = 10,
  n_cores = 8)

mixture_coord_ex_gaussian(
  n_runs = 30,
  q = 3,
  n_random_starts = 100,
  X = NULL,
  order = 3,
  opt_method = "B",
  max_it = 10,
  tol = 0.0001,
  n_cox_points = NULL,
  plot_designs = T,
  verbose = 0,
  opt_crit = "I",
  seed = 10,
  n_cores = 8)












X_mnl_1 = opdesmixr::create_random_initial_MNL_design(q = 3, J = 2, S = 20, seed = 4)





efficiency_cox_scheffe_mnl = function(theta, X, beta, i, j, s, opt_crit){
  # // Computes efficiency criterion of a design matrix X but where the j-th ingredient in i-th observation is changed to theta.
  # // theta must be between 0 and 1 because it's an ingredient proportion.
  # // j and i are 0-indexed.
  # // We want to minimize this.
  #
  # // Create new matrix Y that is identical to the one pointed by X.
  # // Note: This is the easiest way to do it because we have to modify a row in this matrix.
  # // A more computationally effective way would be to only store the new modified vector since
  # // we don't need a copy of the whole matrix. But to do that I would have to either modify some
  # // existing functions, or create some new ones, or both. IDK if the gain in performance is worth it.


  q = dim(X)[1]

  if(opt_crit == 0){
    # "D-optimality"
    W = matrix(0.0, nrow = 1)
  } else{
    # "I-optimality")
    W = create_moment_matrix_gaussian(q)
  }

  beta_mat = matrix(beta, byrow = T, nrow = 1)

  return(efficiencyCoxScheffeMNL(theta, X, beta_mat, i, j, s, opt_crit, W))

}

q = 3
beta_1 = rep(0, (q*q*q + 5*q)/6)

get_opt_crit_value_MNL(X_mnl_1, beta_1, opt_crit = 0)
efficiency_cox_scheffe_mnl(theta = 0.6593244, X = X_mnl_1, beta = beta_1, i = 0, j = 0, s = 0, opt_crit = 0)

efficiency_cox_scheffe_mnl(theta = 0.1, X = X_mnl_1, beta = beta_1, i = 0, j = 0, s = 0, opt_crit = 0)
efficiency_cox_scheffe_mnl(theta = 1, X = X_mnl_1, beta = beta_1, i = 0, j = 0, s = 0, opt_crit = 0)


# We want to minimize this
tibble(theta = 0:50/50) %>%
  mutate(u = map_dbl(theta, function(x) efficiency_cox_scheffe_mnl(theta = x, X = X_mnl_1, beta = beta_1, i = 0, j = 0, s = 0, opt_crit = 0))) %>%
  ggplot() +
  geom_line(aes(theta, u))










q = 3
J = 3
S = 10
X_mnl_2 = create_random_initial_MNL_design(q, J, S, seed = 3)
beta4_2 = create_random_beta(q)

X4_2_opt_disc = mixture_coord_ex(
  X = X_mnl_2,
  beta = beta4_2$beta,
  model = "MNL",
  n_cox_points = 100,
  max_it = 10,
  verbose = 1,
  plot_designs = F,
  opt_crit = "D",
  opt_method = "D"
)



X4_2_opt_brent = mixture_coord_ex(
  X = X_mnl_2,
  beta = beta4_2$beta,
  model = "MNL",
  max_it = 15,
  verbose = 1,
  plot_designs = F,
  opt_crit = "D",
  opt_method = "B"
)




# We want to minimize this
tibble(theta = 0:50/50) %>%
  mutate(u = map_dbl(theta, function(x) efficiency_cox_scheffe_mnl(theta = x, X = X_mnl_2, beta = beta4_2$beta, i = 0, j = 0, s = 0, opt_crit = 0))) %>%
  ggplot() +
  geom_line(aes(theta, u))




# We want to minimize this
tibble(theta = 0:50/50) %>%
  mutate(u = map_dbl(theta, function(x) efficiency_cox_scheffe_mnl(theta = x, X = X4_2_opt_brent$X, beta = beta4_2$beta, i = 0, j = 0, s = 0, opt_crit = 0))) %>%
  ggplot() +
  geom_line(aes(theta, u))











