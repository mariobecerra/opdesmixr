library(tidyverse)
library(here)

theme_set(theme_bw())

# Rcpp::sourceCpp("src/utils.cpp")
# source("R/mnl_model_functions.R")
# source("R/general_functions.R")

out_dir_1 = "tests/i_opt_d_opt_plots/"
out_dir = paste0(out_dir_1, "rds_files/")

file_names = list.files(here(out_dir))

efficiencies = map_df(seq_along(file_names), function(i){
  eff = readRDS(here(out_dir, file_names[i]))
  if(any(is.na(eff))) eff = eff %>% filter(rep(F, nrow(eff)))
  return(eff)
})


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



