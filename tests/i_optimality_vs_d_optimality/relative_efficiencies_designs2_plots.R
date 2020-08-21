library(tidyverse)
library(here)

# Had to separate this fucking file because ggtern fucks up ggplot2 and IDK why yet
# Something's wrong in D-optimality

ggplot2::theme_set(ggplot2::theme_bw())

out_dir_1 = "tests/i_optimality_vs_d_optimality/out/"
out_dir_2 = paste0(out_dir_1, "rds_files/")


efficiencies = readRDS(here(out_dir_2, "efficiencies_all.rds"))


efficiencies %>%
  filter(complete.cases(.)) %>%
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
  filter(complete.cases(.)) %>%
  mutate(group = paste0("q = ", q, ", J = ", J, ", S = ", S)) %>%
  ggplot() +
  geom_histogram(aes(ratio_d_eff)) +
  facet_wrap(~group, scales = "free")


efficiencies %>%
  filter(complete.cases(.)) %>%
  group_by(q, J, S) %>%
  mutate(n_runs = n()) %>%
  ungroup() %>%
  filter(n_runs > 10) %>%
  mutate(group = paste0("q = ", q, ", J = ", J, ", S = ", S)) %>%
  ggplot() +
  geom_boxplot(aes(group, ratio_d_eff)) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  ggplot2::coord_flip() +
  xlab("") +
  ylab("Ratio of D-efficiencies")



efficiencies %>%
  filter(complete.cases(.)) %>%
  mutate(group = paste0("q = ", q, ", J = ", J, ", S = ", S)) %>%
  ggplot() +
  geom_boxplot(aes(group, ratio_i_eff)) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  coord_flip() +
  xlab("") +
  ylab("Ratio of I-efficiencies")























