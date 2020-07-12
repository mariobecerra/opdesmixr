library(tidyverse)
library(here)

devtools::load_all(".")

ggplot2::theme_set(ggplot2::theme_bw())

# Rcpp::sourceCpp("src/utils.cpp")
# source("R/mnl_model_functions.R")
# source("R/general_functions.R")


out_dir_1 = "tests/i_optimality_vs_d_optimality/out/"
out_dir_2 = paste0(out_dir_1, "rds_files/")
out_dir = paste0(out_dir_2, "many_random_designs/")

file_names = list.files(out_dir)


efficiencies = map_df(seq_along(file_names), function(i){

  filename_i = file_names[i]

  q = str_extract(pattern = "[0-9]+", str_extract(pattern = "_.*?_", string = filename_i))
  J = str_extract(pattern = "[0-9]+", str_extract(pattern = "J.*?_", string = filename_i))
  S = str_extract(pattern = "[0-9]+", str_extract(pattern = "S.*?_", string = filename_i))

  designs = readRDS(here(out_dir, filename_i))
  D_opt = designs$D_opt
  I_opt = designs$I_opt
  beta = designs$beta

  d_eff_of_d_optimal_design = try(get_opt_crit_value_MNL(D_opt$X, beta = beta, opt_crit = 0), silent = T)
  d_eff_of_i_optimal_design = try(get_opt_crit_value_MNL(I_opt$X, beta = beta, opt_crit = 0), silent = T)
  i_eff_of_d_optimal_design = try(get_opt_crit_value_MNL(D_opt$X, beta = beta, opt_crit = 1), silent = T)
  i_eff_of_i_optimal_design = try(get_opt_crit_value_MNL(I_opt$X, beta = beta, opt_crit = 1), silent = T)


  if(
    class(d_eff_of_d_optimal_design) == 'try-error' |
    class(d_eff_of_i_optimal_design) == 'try-error' |
    class(i_eff_of_d_optimal_design) == 'try-error' |
    class(i_eff_of_i_optimal_design) == 'try-error'
  ){
    d_eff_of_d_optimal_design = NA_real_
    d_eff_of_i_optimal_design = NA_real_
    i_eff_of_d_optimal_design = NA_real_
    i_eff_of_i_optimal_design = NA_real_
  }

  out = tibble(
    d_eff_of_d_optimal_design = d_eff_of_d_optimal_design,
    d_eff_of_i_optimal_design = d_eff_of_i_optimal_design,
    i_eff_of_d_optimal_design = i_eff_of_d_optimal_design,
    i_eff_of_i_optimal_design = i_eff_of_i_optimal_design
  ) %>%
    mutate(q = as.numeric(q), J = J, S = S)

  return(out)

}) %>%
  mutate(n_params = (q^3 + 5*q)/6-1) %>%
  mutate(
    ratio_d_eff = (exp(d_eff_of_i_optimal_design - d_eff_of_d_optimal_design))^(1/n_params),
    ratio_i_eff = (exp(i_eff_of_d_optimal_design - i_eff_of_i_optimal_design))
  )



efficiencies %>%
  group_by(q, J, S) %>%
  summarize(
    n = n(),
    sum_d_eff_diff = sum(d_eff_of_d_optimal_design < d_eff_of_i_optimal_design),
    sum_i_eff_diff = sum(i_eff_of_i_optimal_design < i_eff_of_d_optimal_design)) %>%
  mutate(
    group = paste0("q =", q, ", J =", J, ", S =", S),
    perc_d_eff_diff = sum_d_eff_diff/n,
    perc_i_eff_diff = sum_i_eff_diff/n
  ) %>%
  pivot_longer(cols = starts_with("perc")) %>%
  ggplot() +
  geom_bar(aes(name, value), stat = "identity") +
  facet_wrap(~group) +
  xlab("") +
  ylab("Percentage") +
  geom_hline(yintercept = 1)



# I like the boxplot better
efficiencies %>%
  mutate(group = paste0("q = ", q, ", J = ", J, ", S = ", S)) %>%
  ggplot() +
  geom_histogram(aes(ratio_d_eff)) +
  facet_wrap(~group, scales = "free")


efficiencies %>%
  filter(complete.cases(.)) %>%
  mutate(group = paste0("q = ", q, ", J = ", J, ", S = ", S)) %>%
  ggplot() +
  geom_boxplot(aes(group, ratio_d_eff)) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  ggplot2::coord_flip() +
  xlab("") +
  ylab("Ratio of D-efficiencies")



efficiencies %>%
  mutate(group = paste0("q = ", q, ", J = ", J, ", S = ", S)) %>%
  ggplot() +
  geom_boxplot(aes(group, ratio_i_eff)) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  coord_flip() +
  xlab("") +
  ylab("Ratio of I-efficiencies")























