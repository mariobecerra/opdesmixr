# library(opdesmixr)
library(tidyverse)
library(here)

Rcpp::sourceCpp("src/utils.cpp")
source("R/mnl_model_functions.R")
source("R/general_functions.R")

out_dir_1 = "tests/i_optimality_vs_d_optimality/out/"
out_dir = paste0(out_dir_1, "rds_files/")
dir.create(here(out_dir_1))
dir.create(here(out_dir))

n_rand_starts = 100
n_cox_points = 40
max_it = 3
n_cores = parallel::detectCores()
n_designs = 50

for(S in 2:7){
  for(J in 2:5){
    for(q in 5:3){
      # for(S in 5:7){
      #   for(J in 3:4){
      #     for(q in 4:3){

      time_now = substring(as.character(Sys.time()), 12, 1000)
      cat("\n\n\n\n", "q = ", q, ", J = ", J, ", S = ", S, " (", time_now, ")", "\n", sep = "")

      file_name = here(out_dir, paste0("efficiencies_q", q, "_J", J, "_S", S, ".rds"))
      if(file.exists(file_name)){
        cat(file_name, "already exists.")
      } else{
        efficiencies = lapply(1:n_designs, function(i){

          # Cada iteración de estas tarda como 18 segundos con 4 cores
          # Si n_designs = 50, y tengo 2 qs 2 Ss y 2 Js, entonces tardará 50*2*2*2 = 400 seg ~ 7 mins
          time_now = substring(as.character(Sys.time()), 12, 1000)
          cat("\nIter ", i, " (", time_now, ")\n", sep = "")

          beta = create_random_beta(q)$beta

          D_opt = mixture_coord_ex(
            n_random_starts = n_rand_starts,
            q = q,
            J = J,
            S = S,
            beta = beta,
            model = "MNL",
            n_cox_points = n_cox_points,
            max_it = max_it,
            verbose = 0,
            plot_designs = F,
            opt_crit = 0,
            seed = i,
            n_cores = n_cores
          )

          I_opt = mixture_coord_ex(
            n_random_starts = n_rand_starts,
            q = q,
            J = J,
            S = S,
            beta = beta,
            model = "MNL",
            n_cox_points = n_cox_points,
            max_it = max_it,
            verbose = 0,
            plot_designs = F,
            opt_crit = 1,
            seed = i,
            n_cores = n_cores
          )

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
            mutate(q = q, J = J, S = S)

          return(out)

        }) %>%
          bind_rows()

        cat("\tSaving RDS...")
        saveRDS(efficiencies, file_name)
        cat("Done.\n")


      } # end else


    }
  }
}









# ggplot2::theme_set(ggplot2::theme_bw())
#
# q1 = 3
# J1 = 5
# S1 = 4
#
# n_rand_starts_1 = 100
# n_cox_points_1 = 50
# max_it_1 = 4
# n_cores_1 = 4
#
# time_init = Sys.time()
# efficiencies_3_5_4 = lapply(1:5, function(i){
#
#   # Cada iteración de estas tarda como 18 segundos
#   cat("\n\nIter:", i, "\n")
#
#   beta1_1 = create_random_beta(q1)$beta
#
#   D_opt_1 = mixture_coord_ex(
#     n_random_starts = n_rand_starts_1,
#     q = q1,
#     J = J1,
#     S = S1,
#     beta = beta1_1,
#     model = "MNL",
#     n_cox_points = n_cox_points_1,
#     max_it = max_it_1,
#     verbose = 0,
#     plot_designs = F,
#     opt_crit = 0,
#     seed = i,
#     n_cores = n_cores_1
#   )
#
#   I_opt_1 = mixture_coord_ex(
#     n_random_starts = n_rand_starts_1,
#     q = q1,
#     J = J1,
#     S = S1,
#     beta = beta1_1,
#     model = "MNL",
#     n_cox_points = n_cox_points_1,
#     max_it = max_it_1,
#     verbose = 0,
#     plot_designs = F,
#     opt_crit = 1,
#     seed = i,
#     n_cores = n_cores_1
#   )
#
#   out = tibble(
#     d_eff_of_d_optimal_design = get_opt_crit_value_MNL(D_opt_1$X, beta = beta1_1, opt_crit = 0),
#     d_eff_of_i_optimal_design = get_opt_crit_value_MNL(I_opt_1$X, beta = beta1_1, opt_crit = 0),
#     i_eff_of_d_optimal_design = get_opt_crit_value_MNL(D_opt_1$X, beta = beta1_1, opt_crit = 1),
#     i_eff_of_i_optimal_design = get_opt_crit_value_MNL(I_opt_1$X, beta = beta1_1, opt_crit = 1)
#   )
#
#   return(out)
#
# }) %>%
#   bind_rows() %>%
#   mutate(q = q1, J = J1, S = S1)
# time_end = Sys.time()
# time_end - time_init
#
#
#
# efficiencies_3_5_4 %>%
#   group_by(q, J, S) %>%
#   summarize(
#     n = n(),
#     sum_d_eff_diff = sum(d_eff_of_d_optimal_design < d_eff_of_i_optimal_design),
#     sum_i_eff_diff = sum(i_eff_of_i_optimal_design < i_eff_of_d_optimal_design)) %>%
#   mutate(
#     group = paste("q =", q, "J =", J, "S =", S),
#     porc_d_eff_diff = sum_d_eff_diff/n,
#     porc_i_eff_diff = sum_i_eff_diff/n
#   ) %>%
#   pivot_longer(cols = starts_with("porc")) %>%
#   ggplot() +
#   geom_bar(aes(name, value), stat = "identity") +
#   facet_wrap(~group) +
#   theme(
#     axis.text.x = element_blank(),
#     axis.text.y = element_blank(),
#     axis.ticks = element_blank())












# q1 = 3
# J1 = 5
# S1 = 3
#
# beta1_0 = rep(0, (q1*q1*q1 + 5*q1)/6)
#
# D_opt_0 = mixture_coord_ex(
#   n_random_starts = 100,
#   q = q1,
#   J = J1,
#   S = S1,
#   beta = beta1_0,
#   model = "MNL",
#   n_cox_points = 50,
#   max_it = 5,
#   verbose = 1,
#   plot_designs = T,
#   opt_crit = 0,
#   seed = 1,
#   n_cores = 4
# )
#
# I_opt_0 = mixture_coord_ex(
#   n_random_starts = 100,
#   q = q1,
#   J = J1,
#   S = S1,
#   beta = beta1_0,
#   model = "MNL",
#   n_cox_points = 50,
#   max_it = 5,
#   verbose = 0,
#   plot_designs = T,
#   opt_crit = 1,
#   seed = 1,
#   n_cores = 4
# )












