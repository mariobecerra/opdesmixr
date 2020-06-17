library(tidyverse)
library(here)
# library(opdesmixr)
# devtools::load_all(".")


theme_set(theme_bw())


out_dir_1 = "tests/i_optimality_vs_d_optimality/out/"
out_dir = paste0(out_dir_1, "rds_files/")
dir.create(here(out_dir_1))
dir.create(here(out_dir))

file_name = here(out_dir, "fds_plots.rds")

if(file.exists(file_name)){
  cat("File already exists. Reading RDS.\n")

  pred_vars = readRDS(file_name)

} else{

  Rcpp::sourceCpp("src/utils.cpp")
  source("R/mnl_model_functions.R")
  source("R/general_functions.R")

  test_ix = expand_grid(q = 3:4, J = 2:3, S = c(10, 20))

  test_ix = tibble(
    q = c(rep(3, 4), rep(4, 4)),
    J = c(rep(2, 2), rep(3, 2), rep(3, 2), rep(4, 2)),
    S = rep(c(10, 20), 4))

  pred_vars = lapply(1:nrow(test_ix), function(i){


    q1 = test_ix$q[i]
    J1 = test_ix$J[i]
    S1 = test_ix$S[i]


    group = paste0(
      "q = ", q1, ", ",
      "J = ", J1, ", ",
      "S = ", S1
    )

    cat("Doing group", group, "(", i, "out of", nrow(test_ix), ")\n")

    beta = create_random_beta(q1, seed = i)$beta

    D_opt_0 = mixture_coord_ex(
      n_random_starts = 50,
      q = q1,
      J = J1,
      S = S1,
      beta = beta,
      model = "MNL",
      n_cox_points = 30,
      max_it = 5,
      verbose = 0,
      plot_designs = F,
      opt_crit = 0,
      seed = 1,
      n_cores = 1
    )

    I_opt_0 = mixture_coord_ex(
      n_random_starts = 50,
      q = q1,
      J = J1,
      S = S1,
      beta = beta,
      model = "MNL",
      n_cox_points = 30,
      max_it = 5,
      verbose = 0,
      plot_designs = F,
      opt_crit = 1,
      seed = 1,
      n_cores = 1
    )


    designs = list(
      D_opt = D_opt_0$X, I_opt =  I_opt_0$X
    )

    # Variance of fraction of design space
    pred_vars = lapply(seq_along(designs), function(k){

      X = designs[[k]]
      inf_mat = getInformationMatrixMNL(X, beta = beta)

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
      bind_rows() %>%
      mutate(group = group)



    return(pred_vars)

  }) %>%
    bind_rows()

  saveRDS(pred_vars, file_name)

}





medians = pred_vars %>%
  group_by(group, Design) %>%
  summarize(
    med = median(pred_var),
    mean = mean(pred_var))



pred_vars %>%
  left_join(
    medians
  ) %>%
  ggplot() +
  geom_vline(xintercept = 0.5, linetype = "dashed", size = 0.5) +
  geom_hline(aes(yintercept = med), linetype = "dashed", size = 0.5) +
  geom_line(aes(fraction, pred_var, color = Design), size = 0.8) +
  xlab("Fraction of design space") +
  ylab("Prediction variance") +
  facet_wrap(~ group, scales = "free_y", nrow = 2)


