library(tidyverse)
library(opdesmixr)
library(randtoolbox)
library(here)


q = 3
J = 2
S = 16

m = (q^3 + 5*q)/6

beta = rep(0, m)
kappa = 5



### Bayesian
draws_unif_beta = halton(120, dim = length(beta))

beta_prior_draws = matrix(rep(NA_real_, length(draws_unif_beta)), ncol = ncol(draws_unif_beta))
for(i in 1:ncol(beta_prior_draws)){
  beta_prior_draws[, i] = qnorm(draws_unif_beta[, i], mean = beta[i], sd = sqrt(kappa))
}


beta_prior_draws %>%
  as.data.frame() %>%
  pivot_longer(starts_with("V")) %>%
  ggplot() +
  geom_density(aes(value)) +
  facet_wrap(~name)






n_cores = parallel::detectCores()

n_cox_points = 10
n_rand_starts = 50
max_it = 5
seed = 2020


# 4 mins
(t1D = Sys.time())
D_opt = mixture_coord_ex(
  n_random_starts = n_rand_starts,
  q = q,
  J = J,
  S = S,
  beta = beta_prior_draws,
  model = "MNL",
  n_cox_points = n_cox_points,
  max_it = max_it,
  verbose = 0,
  plot_designs = T,
  opt_crit = 0,
  seed = seed,
  n_cores = n_cores
)
(t2D = Sys.time())
t2D - t1D

mnl_plot_result(D_opt)





# 5.5 mins
(t1I = Sys.time())
I_opt = mixture_coord_ex(
  n_random_starts = n_rand_starts,
  q = q,
  J = J,
  S = S,
  beta = beta_prior_draws,
  model = "MNL",
  n_cox_points = n_cox_points,
  max_it = max_it,
  verbose = 0,
  plot_designs = F,
  opt_crit = 1,
  seed = seed,
  n_cores = n_cores
)
(t2I = Sys.time())
t2I - t1I

mnl_plot_result(I_opt)




