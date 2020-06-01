library(opdesmixr)
library(tidyverse)


q = 3
J = 5
S = 3
X2 = create_random_initial_MNL_design(q, J, S, seed = 4)
beta2 = rep(0, (q*q*q + 5*q)/6)



X2_opt2 = mixture_coord_ex_mnl(
  X = X2,
  beta = beta2,
  n_cox_points = 5,
  max_it = 2,
  verbose = 5,
  plot_designs = T,
  n_cores = 1
)









q = 3
J = 5
S = 4
X4 = create_random_initial_MNL_design(q, J, S, seed = 3)
beta4_2 = create_random_beta(q)

(t1 = Sys.time())
X_q3_J5_s4_D_2 = mixture_coord_ex_mnl(
  n_random_starts = 500,
  q = q,
  J = J,
  S = S,
  beta = beta4_2$beta,
  n_cox_points = 50,
  max_it = 5,
  verbose = 1,
  plot_designs = T,
  n_cores = 4,
  seed = 10,
  opt_crit = 0
)
(t2 = Sys.time())
t2 - t1















