library(randtoolbox)
library(tidyverse)
devtools::load_all(".")

# Monte Carlo integration

joint.sam <- cbind(runif(10000, min = -1, max = 1),
                   runif(10000, min = -1, max = 1))

g.indA <- function(point){
  if ((point[1]^2 + point[2]^2) <= 1) return (1.0)
  else return(0)
}

theta.mc <- mean(apply(joint.sam, 1, g.indA))
4*theta.mc








monte_carlo0 = function(n1, n2){
  mat = matrix(rep(NA_real_, n1*n2), nrow = n1)
  for(i in 1:n1){
    for(j in 1:n2){
      X <- runif(1, 0, 1)
      Y <- runif(1, 0, 1)

      mat[i, j] = X*Y
    }
  }
  return(mean(mat))
}


monte_carlo0(500, 500)







monte_carlo = function(fun, n){
  mat = matrix(runif(2*n, 0, 1), ncol = 2)
  acc = 0.0
  for(i in 1:n){
    for(j in 1:n){
      acc = acc + fun(mat[i, 1], mat[j, 2])
    }
  }

  return(acc/n^2)
}


monte_carlo(function(x, y) x*y, 1000) # 0.25
monte_carlo(function(x, y) x*y^2, 1000) # 0.1666667









# Halton sequences




sims_unif = halton(100)

qplot(qnorm(sims_unif))




qplot(qnorm(halton(100)))
qplot(qnorm(ghalton(100)))
qplot(rnorm(100))



qplot(qnorm(halton(1000)))
qplot(qnorm(ghalton(1000)))
qplot(rnorm(1000))


qplot(qnorm(halton(1000), 2, 3))






sims_unif_6 = halton(120, dim = 6)
sims_norm_6 = matrix(rep(NA_real_, length(sims_unif_6)), ncol = ncol(sims_unif_6))
for(i in 1:ncol(sims_norm_6)){
  sims_norm_6[, i] = qnorm(sims_unif_6[, i], 2, 3)
}

sims_norm_6 %>%
  as.data.frame() %>%
  pivot_longer(V1:V6) %>%
  ggplot() +
  geom_histogram(aes(value)) +
  facet_wrap(~name)


sims_norm_6 %>%
  as.data.frame() %>%
  pivot_longer(V1:V6) %>%
  ggplot() +
  geom_density(aes(value)) +
  facet_wrap(~name)










n1 = 100
results_100 = microbenchmark::microbenchmark(
  a = qnorm(halton(n1)),
  b = qnorm(ghalton(n1)),
  c = rnorm(n2),
  times = 10000
)

results_100



n2 = 1000
results_1000 = microbenchmark::microbenchmark(
  a = qnorm(halton(n2)),
  b = qnorm(ghalton(n2)),
  c = rnorm(n2),
  times = 1000
)

results_1000








n3 = 10000
results_10000 = microbenchmark::microbenchmark(
  a = qnorm(halton(n3)),
  b = qnorm(ghalton(n3)),
  c = rnorm(n3),
  times = 1000
)

results_10000
















monte_carlo_halton = function(fun, n){
  mat = halton(n, 2)
  acc = 0.0
  for(i in 1:n){
    for(j in 1:n){
      acc = acc + fun(mat[i, 1], mat[j, 2])
    }
  }

  return(acc/n^2)
}


monte_carlo_halton(function(x, y) x*y, 100) # 0.25
monte_carlo(function(x, y) x*y, 100) # 0.25

monte_carlo_halton(function(x, y) x*y^2, 100) # 0.1666667
monte_carlo(function(x, y) x*y^2, 100) # 0.1666667



monte_carlo_halton(function(x, y) x*y, 1000) # 0.25
monte_carlo_halton(function(x, y) x*y^2, 1000) # 0.1666667













q = 3
J = 5
S = 4
X4 = create_random_initial_MNL_design(q, J, S, seed = 3)
beta4_2 = create_random_beta(q)





(t1 = Sys.time())
X_q3_J5_s4_D_2 = mixture_coord_ex_mnl(
  n_random_starts = 100,
  q = q,
  J = J,
  S = S,
  beta = beta4_2$beta,
  n_cox_points = 20,
  max_it = 5,
  verbose = 1,
  plot_designs = T,
  n_cores = 4,
  seed = 10,
  opt_crit = 0
)
(t2 = Sys.time())
t2 - t1








### Bayesian
draws_unif_beta4_2 = halton(120, dim = length(beta4_2$beta))

beta4_2_prior_draws = matrix(rep(NA_real_, length(draws_unif_beta4_2)), ncol = ncol(draws_unif_beta4_2))
for(i in 1:ncol(beta4_2_prior_draws)){
  beta4_2_prior_draws[, i] = qnorm(draws_unif_beta4_2[, i], beta4_2$beta[i], 3)
}


beta4_2_prior_draws %>%
  as.data.frame() %>%
  pivot_longer(starts_with("V")) %>%
  ggplot() +
  geom_density(aes(value)) +
  facet_wrap(~name)





# 1.7 mins
(t1 = Sys.time())
X_q3_J5_s4_D_bayes = mixture_coord_ex_mnl(
  n_random_starts = 8,
  q = q,
  J = J,
  S = S,
  beta = beta4_2_prior_draws,
  n_cox_points = 20,
  max_it = 4,
  verbose = 1,
  plot_designs = T,
  n_cores = 4,
  seed = 10,
  opt_crit = 0
)
(t2 = Sys.time())
t2 - t1








