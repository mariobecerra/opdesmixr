library(tidyverse)
library(Rcpp)
# library(RcppBrent)

# sourceCpp("tests/brent_test.cpp")
# devtools::load_all(".")
library(opdesmixr)

devtools::unload("ggtern")

R.methodsS3::setMethodS3("print", "ggplot", ggplot2:::print.ggplot)
R.methodsS3::setMethodS3("plot", "ggplot", ggplot2:::plot.ggplot)
R.methodsS3::setMethodS3("grid.draw", "ggplot", ggplot2:::grid.draw.ggplot)

# John Denker's implementation: https://people.sc.fsu.edu/~jburkardt/cpp_src/brent/brent.html
# brent.cpp: https://people.sc.fsu.edu/~jburkardt/cpp_src/brent/brent.cpp
# brent.hpp: https://people.sc.fsu.edu/~jburkardt/cpp_src/brent/brent.hpp





tibble(x = seq(-1, 1, length.out = 200)) %>%
  mutate(y = x^2) %>%
  ggplot() +
  geom_line(aes(x, y))


tibble(x = seq(-7, 15, length.out = 200)) %>%
  mutate(y = map_dbl(x, RcppBrent::f)) %>%
  ggplot() +
  geom_line(aes(x, y))


RcppBrent::min_f()







tibble(x = seq(-10, 10, length.out = 200)) %>%
  mutate(y = map_dbl(x, f2)) %>%
  ggplot() +
  geom_line(aes(x, y))

tibble(x = seq(-1, 5, length.out = 200)) %>%
  mutate(y = map_dbl(x, f3)) %>%
  ggplot() +
  geom_line(aes(x, y))



min_f2()
min_f2_2()
min_f3()










## Banana function


tibble(x = seq(-10, 10, length.out = 200)) %>%
  mutate(y = map_dbl(x, banana_x_y1)) %>%
  ggplot() +
  geom_line(aes(x, y))



# Local minimum at x = -1
tibble(x = seq(-1.05, -0.95, length.out = 200)) %>%
  mutate(y = map_dbl(x, banana_x_y1)) %>%
  ggplot() +
  geom_line(aes(x, y))
banana_x_y1(-1)


# Global minimum at x = 1
tibble(x = seq(0.95, 1.05, length.out = 200)) %>%
  mutate(y = map_dbl(x, banana_x_y1)) %>%
  ggplot() +
  geom_line(aes(x, y))
banana_x_y1(1)




min_banana_x_y1()
min_banana_x_y1(-10, 10)
min_banana_x_y1(-1.02, 1.02)




tibble(x = seq(-10, 10, length.out = 200)) %>%
  mutate(y = map_dbl(x, function(.) banana_xy(., 2))) %>%
  ggplot() +
  geom_line(aes(x, y))

minimize_banana_fixed_y(1)


tibble(x = seq(1.35, 1.45, length.out = 200)) %>%
  mutate(y = map_dbl(x, function(.) banana_xy(., 2))) %>%
  ggplot() +
  geom_line(aes(x, y))

minimize_banana_fixed_y(2)






library(tidyverse)
library(opdesmixr)
# devtools::load_all(".")

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
res_alg_order_3_2 = mixture_coord_ex(
  X = X2,
  order = 3,
  model = "Gaussian",
  plot_designs = F,
  n_cox_points = 100,
  max_it = 2)

X3 = res_alg_order_3_2$X


tibble(theta = 0:50/50) %>%
  mutate(u = map_dbl(theta, function(x) efficiency_cox_scheffe_gaussian(theta = x, X = X3, j = 0, i = 0, order = 3, opt_crit = 0))) %>%
  ggplot() +
  geom_line(aes(theta, u))








