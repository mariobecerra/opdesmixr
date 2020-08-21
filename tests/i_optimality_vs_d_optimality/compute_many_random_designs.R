library(opdesmixr)
library(tidyverse)
library(here)

# Rcpp::sourceCpp("src/utils.cpp")
# source("R/mnl_model_functions.R")
# source("R/general_functions.R")



out_dir_1 = "tests/i_optimality_vs_d_optimality/out/"
out_dir_2 = paste0(out_dir_1, "rds_files/")
out_dir = paste0(out_dir_2, "many_random_designs/")
dir.create(here(out_dir_1))
dir.create(here(out_dir_2))
dir.create(here(out_dir))

n_rand_starts = 50# n_rand_starts = 100
n_cox_points = 50
max_it = 3
n_cores = parallel::detectCores()
n_designs = 30

for(J in 3){
  for(S in c(35)){
    # for(S in c(25, 35)){ #for(S in c(10, 6, 15)){
    for(q in 6) {
      #for(q in 5:6){ #for(q in 3:5){

      time_now = substring(as.character(Sys.time()), 12, 1000)
      cat("\n\n\n\n", "q = ", q, ", J = ", J, ", S = ", S, " (", time_now, ")", "\n", sep = "")



      for(i in 1:n_designs){

        # Cada iteración de estas tarda como 10 segundos con 4 cores
        # Si n_designs = 50, y tengo 2 qs, 1 Ss y 2 Js, entonces tardará 10*50*2*2*1 = 2000 seg ~ 33 minutos
        # Al final tardó como 25 minutos
        time_now = substring(as.character(Sys.time()), 12, 1000)
        cat("\n\tIter ", i, " (", time_now, ", q = ", q, ", J = ", J, ", S = ", S, ")\n", sep = "")


        beta = create_random_beta(q)$beta
        # beta_fname = gsub("\\.|\\-", "", paste0(beta, collapse = "_"))
        beta_fname = gsub("\\.|\\-", "", runif(1))

        file_name = here(out_dir, paste0("design_", q, "_J", J, "_S", S, "_", i, "_", beta_fname, ".rds"))



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

        out = list(
          D_opt = D_opt,
          I_opt = I_opt,
          beta = beta
        )


        cat("\tD_opt_value = ", D_opt$opt_crit_value, "\n",
            "\tI_opt_value = ", I_opt$opt_crit_value, "\n")

        cat("\tSaving RDS...")
        saveRDS(out, file_name)
        cat("Done.\n")

      } # end for (i)



    }
  }
}
