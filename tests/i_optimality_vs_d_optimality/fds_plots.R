library(tidyverse)
# library(opdesmixr)
# devtools::load_all(".")
Rcpp::sourceCpp("src/utils.cpp")
source("R/mnl_model_functions.R")
source("R/general_functions.R")

theme_set(theme_bw())

q1 = 3
J1 = 5
S1 = 3

beta1_0 = rep(0, (q1*q1*q1 + 5*q1)/6)

D_opt_0 = mixture_coord_ex(
  n_random_starts = 100,
  q = q1,
  J = J1,
  S = S1,
  beta = beta1_0,
  model = "MNL",
  n_cox_points = 50,
  max_it = 5,
  verbose = 1,
  plot_designs = F,
  opt_crit = 0,
  seed = 1,
  n_cores = 4
)

I_opt_0 = mixture_coord_ex(
  n_random_starts = 100,
  q = q1,
  J = J1,
  S = S1,
  beta = beta1_0,
  model = "MNL",
  n_cox_points = 50,
  max_it = 5,
  verbose = 0,
  plot_designs = F,
  opt_crit = 1,
  seed = 1,
  n_cores = 4
)


designs = list(
  D_opt = D_opt_0$X, I_opt =  I_opt_0$X
)

# Variance of fraction of design space
pred_vars = lapply(seq_along(designs), function(k){

  X = designs[[k]]
  inf_mat = getInformationMatrixMNL(X, beta = beta1_0)

  pred_var = lapply(1:1000, function(i){
    des_i = create_random_initial_MNL_design(q1, J1, S1, seed = i)


    vars_1 = rep(NA_real_, J1)
    for(j in 1:J1){
      f_x = getXsMNL(des_i, 1)[j,]
      vars_1[j] = t(f_x) %*% solve(inf_mat, f_x)
    }

    return(vars_1)
  }) %>%
    unlist()



  out = tibble(Design = names(designs)[k],
               pred_var = sort(pred_var)) %>%
    mutate(fraction = 1:nrow(.)/nrow(.))

  return(out)
}) %>%
  bind_rows()




medians = pred_vars %>%
  group_by(Design) %>%
  summarize(
    med = median(pred_var),
    mean = mean(pred_var))


pred_vars %>%
  ggplot() +
  geom_vline(xintercept = 0.5, linetype = "dashed", size = 0.5) +
  geom_hline(yintercept = medians$med, linetype = "dashed", size = 0.5) +
  geom_line(aes(fraction, pred_var, color = Design), size = 0.8) +
  xlab("Fraction of design space") +
  ylab("Prediction variance")


