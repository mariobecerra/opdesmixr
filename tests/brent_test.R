library(tidyverse)
library(Rcpp)
# library(RcppBrent)

sourceCpp("src/brent_test.cpp")

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
  mutate(y = map_dbl(x, banana)) %>%
  ggplot() +
  geom_line(aes(x, y))



# Local minimum at x = -1
tibble(x = seq(-1.05, -0.95, length.out = 200)) %>%
  mutate(y = map_dbl(x, banana)) %>%
  ggplot() +
  geom_line(aes(x, y))
banana(-1)


# Global minimum at x = 1
tibble(x = seq(0.95, 1.05, length.out = 200)) %>%
  mutate(y = map_dbl(x, banana)) %>%
  ggplot() +
  geom_line(aes(x, y))
banana(1)




min_banana()
min_banana(-10, 10)
min_banana(-1.02, 1.02)







# // Iterate over ingredient proportions
# for(int j_aux = 0; j_aux < setDiff.n_elem; j_aux++){
#   j = setDiff(j_aux);
#   if(abs(1 - x(i)) < 1e-16) {
#     // In case x(i) is numerically 1
#     result = (1 - cox_direction(n, i))/(q-1);
#   } else{
#     // In case x(i) is not numerically 1
#     result = x(j) - deltas(n)*x(j)/(1 - x(i));
#   }
#   cox_direction(n, j) = result;
#   j++;
# } // end for j







fn_cox_scheffe = function(theta, X, j, i, order){

  x_row = X[j,]

  delta = theta - x_row[i]
  # recompute proportions in cox dir
  for(k in setdiff(1:q, i)){
    # In case it's a corner case, i.e., x[i] = 1
    # This may be wrong. Check later.
    if(abs(1 - x_row[i]) < 1e-16) res = (1 - x_row[i])/(q-1)
    else{
      res = x_row[k] - delta*x_row[k]/(1 - x_row[i])
    }
    x_row[k] = res
  } # end k
  x_row[i] = theta

  Y = X
  Y[j,] = x_row

  utility_funct = get_scheffe_log_D_efficiency(Y, order = order)
  return(-as.numeric(utility_funct))
}









X = opdesmixr::create_random_initial_design_gaussian(10, 3, seed = 4)
get_opt_crit_value_Gaussian(X, order = 1, opt_crit = 0)
efficiency_cox_scheffe(theta = 0.65932440, X = X, j = 0, i = 0, order = 1, opt_crit = 0, W = matrix(0.0, nrow = 1))
efficiency_cox_scheffe(theta = 1, X = X, j = 0, i = 0, order = 1, opt_crit = 0, W = matrix(0.0, nrow = 1))


X2 = opdesmixr::create_random_initial_design_gaussian(30, 3, seed = 4)
get_opt_crit_value_Gaussian(X2, order = 3, opt_crit = 1)
efficiency_cox_scheffe(theta = 0.65932440, X = X2, j = 0, i = 0, order = 3, opt_crit = 1, W = create_moment_matrix_gaussian(3))
efficiency_cox_scheffe(theta = 0.1, X = X2, j = 0, i = 0, order = 3, opt_crit = 1, W = create_moment_matrix_gaussian(3))
efficiency_cox_scheffe(theta = 1, X = X2, j = 0, i = 0, order = 3, opt_crit = 1, W = create_moment_matrix_gaussian(3))



tibble(theta = 0:50/50) %>%
  mutate(u = -map_dbl(theta, function(x) efficiency_cox_scheffe(theta = x, X = X2, j = 0, i = 0, order = 3, opt_crit = 1, W = create_moment_matrix_gaussian(3)))) %>%
  ggplot() +
  geom_line(aes(theta, u))




tibble(theta = 0:50/50) %>%
  mutate(u = -map_dbl(theta, function(x) efficiency_cox_scheffe(theta = x, X = X2, j = 0, i = 0, order = 3, opt_crit = 0, W = matrix(0.0, nrow = 1)))) %>%
  ggplot() +
  geom_line(aes(theta, u))







