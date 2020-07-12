# I-optimality

library(dplyr)
library(purrr)
library(opdesmixr)





plot_D_vs_I_opt = function(D_opt_design, I_opt_design){

  dim_X_D = dim(D_opt_design)
  dim_X_I = dim(I_opt_design)

  if(any(dim_X_D != dim_X_I)) stop("Dimension missmatch")

  q = dim_X_D[1]
  S = dim_X_D[3]

  if(q != 3) stop("Design must be of 3 ingredients.")

  # Convert 3 dimensional arrays into matrices by vertically binding them
  D_eff_mat = t(D_opt_design[,,1])
  for(s in 2:S){
    D_eff_mat = rbind(D_eff_mat, t(D_opt_design[,,s]))
  }

  I_eff_mat = t(I_opt_design[,,1])
  for(s in 2:S){
    I_eff_mat = rbind(I_eff_mat, t(I_opt_design[,,s]))
  }

  # Plot matrices
  ggtern::grid.arrange(
    D_eff_mat %>%
      dplyr::as_tibble() %>%
      purrr::set_names(c("c1", "c2", "c3")) %>%
      ggtern(aes(c1, c2, c3)) +
      geom_point(shape = "x", size = 4) +
      theme_minimal() +
      theme_nomask() +
      ggtitle(
        label = "D-efficient")
    ,
    I_eff_mat %>%
      dplyr::as_tibble() %>%
      purrr::set_names(c("c1", "c2", "c3")) %>%
      ggtern(aes(c1, c2, c3)) +
      geom_point(shape = "x", size = 4) +
      theme_minimal() +
      theme_nomask() +
      ggtitle(label = "I-efficient")
    ,
    ncol = 2
  )

}



##########################################################
##########################################################
## MNL model
##########################################################
##########################################################



q_1 = 3
J_1 = 6
S_1 = 2

beta_1_neutral = rep(0, (q_1*q_1*q_1 + 5*q_1)/6)

D_opt_neutral_1 = mixture_coord_ex(
  n_random_starts = 100,
  q = q_1,
  J = J_1,
  S = S_1,
  beta = beta_1_neutral,
  model = "MNL",
  n_cox_points = 50,
  max_it = 5,
  verbose = 1,
  plot_designs = T,
  opt_crit = 0,
  seed = 1,
  n_cores = 4
)

I_opt_neutral_1 = mixture_coord_ex(
  n_random_starts = 100,
  q = q_1,
  J = J_1,
  S = S_1,
  beta = beta_1_neutral,
  model = "MNL",
  n_cox_points = 50,
  max_it = 5,
  verbose = 0,
  plot_designs = T,
  opt_crit = 1,
  seed = 1,
  n_cores = 4
)



beta_1_random = create_random_beta(q_1)$beta

D_opt_random_1 = mixture_coord_ex(
  n_random_starts = 100,
  q = q_1,
  J = J_1,
  S = S_1,
  beta = beta_1_random,
  model = "MNL",
  n_cox_points = 50,
  max_it = 5,
  verbose = 0,
  plot_designs = T,
  opt_crit = 0,
  seed = 1,
  n_cores = 4
)

I_opt_random_1 = mixture_coord_ex(
  n_random_starts = 100,
  q = q_1,
  J = J_1,
  S = S_1,
  beta = beta_1_random,
  model = "MNL",
  n_cox_points = 50,
  max_it = 5,
  verbose = 0,
  plot_designs = T,
  opt_crit = 1,
  seed = 1,
  n_cores = 4
)


plot_D_vs_I_opt(D_opt_neutral_1$X, I_opt_neutral_1$X)
plot_D_vs_I_opt(D_opt_random_1$X, I_opt_random_1$X)














q_2 = 3
J_2 = 10
S_2 = 3

beta_2_neutral = rep(0, (q_2*q_2*q_2 + 5*q_2)/6)

D_opt_neutral_2 = mixture_coord_ex(
  n_random_starts = 100,
  q = q_2,
  J = J_2,
  S = S_2,
  beta = beta_2_neutral,
  model = "MNL",
  n_cox_points = 50,
  max_it = 5,
  verbose = 1,
  plot_designs = T,
  opt_crit = 0,
  seed = 1,
  n_cores = 4
)

I_opt_neutral_2 = mixture_coord_ex(
  n_random_starts = 100,
  q = q_2,
  J = J_2,
  S = S_2,
  beta = beta_2_neutral,
  model = "MNL",
  n_cox_points = 50,
  max_it = 5,
  verbose = 0,
  plot_designs = T,
  opt_crit = 1,
  seed = 1,
  n_cores = 4
)



beta_2_random = create_random_beta(q_2)$beta

D_opt_random_2 = mixture_coord_ex(
  n_random_starts = 100,
  q = q_2,
  J = J_2,
  S = S_2,
  beta = beta_2_random,
  model = "MNL",
  n_cox_points = 50,
  max_it = 5,
  verbose = 0,
  plot_designs = T,
  opt_crit = 0,
  seed = 1,
  n_cores = 4
)

I_opt_random_2 = mixture_coord_ex(
  n_random_starts = 100,
  q = q_2,
  J = J_2,
  S = S_2,
  beta = beta_2_random,
  model = "MNL",
  n_cox_points = 50,
  max_it = 5,
  verbose = 0,
  plot_designs = T,
  opt_crit = 1,
  seed = 1,
  n_cores = 4
)



plot_D_vs_I_opt(D_opt_neutral_2$X, I_opt_neutral_2$X)
plot_D_vs_I_opt(D_opt_random_2$X, I_opt_random_2$X)



