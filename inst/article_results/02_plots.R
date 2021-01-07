library(opdesmixr)
library(tidyverse)
library(here)

# get_file_path = function(filepath){
#   # This function is useful to use when developing the package and working with data in the inst/ folder.
#   # If the package hasn't been installed, then system.file() won't work, then it goes to the folder.
#   # It assumes that the working directory is the project directory.
#   path_system_file = system.file(filepath, package = "opdesmixr", mustWork = F)
#   if(path_system_file == ""){
#     out_path = here::here("inst/", filepath)
#   } else{
#     out_path = path_system_file
#   }
#   return(out_path)
# }
#
# source(get_file_path("article_results/R_scripts/functions.R"))

# Assuming this sript is run before installing.
designs_folder = here("inst/misc_output/cocktail_cornell_designs/")

out_folder = here("inst/misc_output/cocktail_cornell_plots/")
dir.create(out_folder, showWarnings = F)

low_color_gradient = "yellow"
high_color_gradient = "red"

##########################################################################################
#### Load designs for cocktail experiment
##########################################################################################

beta0 = c(1.36, 1.57, 2.47, -0.43, 0.50, 1.09)

cocktail_d_opt_filename = paste0(designs_folder, "cocktail_d_optimal.rds")
cocktail_i_opt_filename = paste0(designs_folder, "cocktail_i_optimal.rds")

cocktail_D_opt = readRDS(cocktail_d_opt_filename)
cocktail_I_opt = readRDS(cocktail_i_opt_filename)



##########################################################################################
#### Load designs for Cornell's experiment
##########################################################################################

beta_2 = c(0.86, 0.21, 0, 3.07, 2.34, 3.24, -20.59)
beta_2_prime = c(0.86, 0.21, 3.07, 2.34, 3.24, -20.59)
kappas = c(0.5, 5, 10, 30)

n_random_starts_2 = 128
max_it_bayes = 20

cornell_designs_basefilename_transf = paste0(designs_folder, "cornell_experiment_transformed_betas_rs", n_random_starts_2, "_maxit", max_it_bayes)

cornell_designs_basefilename_untransf = paste0(designs_folder, "cornell_experiment_untransformed_betas_rs", n_random_starts_2, "_maxit", max_it_bayes)




## Read designs with transformed betas and save them in a list
cornell_designs_transf = lapply(kappas, function(k){

  cat("kappa =", k, "\n")

  d_opt_filename = paste0(cornell_designs_basefilename_transf, "_kappa", k, "_Dopt.rds")
  i_opt_filename = paste0(cornell_designs_basefilename_transf, "_kappa", k, "_Iopt.rds")



  if(file.exists(d_opt_filename)){
    cat("\tD optimal file exists. Loading.\n")
    cornell_beta_2_bayesian_d_opt = readRDS(d_opt_filename)
  }else{
    stop("\tD optimal file does not exist.\n")

  }

  if(file.exists(i_opt_filename)){
    cat("\tI optimal file exists. Loading.\n")
    cornell_beta_2_bayesian_i_opt = readRDS(i_opt_filename)
  }else{
    stop("\tI optimal file does not exist.\n")
  }

  return(list(
    d_opt = cornell_beta_2_bayesian_d_opt,
    i_opt = cornell_beta_2_bayesian_i_opt,
    kappa = k
  ))

})


## Read designs with untransformed betas and save them in a list
cornell_designs_untransf = lapply(kappas, function(k){

  cat("kappa =", k, "\n")

  d_opt_filename = paste0(cornell_designs_basefilename_untransf, "_kappa", k, "_Dopt.rds")
  i_opt_filename = paste0(cornell_designs_basefilename_untransf, "_kappa", k, "_Iopt.rds")



  if(file.exists(d_opt_filename)){
    cat("\tD optimal file exists. Loading.\n")
    cornell_beta_2_bayesian_d_opt = readRDS(d_opt_filename)
  }else{
    stop("\tD optimal file does not exist.\n")

  }

  if(file.exists(i_opt_filename)){
    cat("\tI optimal file exists. Loading.\n")
    cornell_beta_2_bayesian_i_opt = readRDS(i_opt_filename)
  }else{
    stop("\tI optimal file does not exist.\n")
  }

  return(list(
    d_opt = cornell_beta_2_bayesian_d_opt,
    i_opt = cornell_beta_2_bayesian_i_opt,
    kappa = k
  ))

})

##########################################################################################
#### Plot cocktail designs
##########################################################################################

### Design plots

levels_bayesian = c(0, .5625, 1.125, 1.6875, 2.25)
utilities_cocktail_bayesian_plot = opdesmixr:::scheffe_order3_q3_utilities(beta0, 100000) %>%
  mutate(utility_fact = cut(utility, levels_bayesian, right = F)) %>%
  mutate(utility_int = as.integer(utility_fact)) %>%
  mutate(utility_avg = case_when(
    utility_int == 1 ~ (levels_bayesian[1] + levels_bayesian[2] - levels_bayesian[1])/2,
    utility_int == 2 ~ (levels_bayesian[2] + levels_bayesian[3] - levels_bayesian[2])/2,
    utility_int == 3 ~ (levels_bayesian[3] + levels_bayesian[4] - levels_bayesian[3])/2,
    utility_int == 4 ~ (levels_bayesian[4] + levels_bayesian[5] - levels_bayesian[4])/2
  ))


width_cocktail_simplex = 10
height_cocktail_simplex = 10
utility_point_size_cocktail = 0.01
utility_point_shape_cocktail = "circle"

ggtern::ggsave(
  filename = paste0(out_folder, "res_cocktail_design_db_simplex.png"),
  plot = opdesmixr:::plot_choice_set_utility(
    cocktail_D_opt,
    utilities_cocktail_bayesian_plot,
    utility_point_size = utility_point_size_cocktail,
    utility_point_shape = utility_point_shape_cocktail,
    low_color_gradient = low_color_gradient,
    high_color_gradient = high_color_gradient
    ) +
    # ggtitle("Bayesian D-optimal for cocktail experiment") +
    theme(legend.position = "none"),
  width = width_cocktail_simplex,
  height = height_cocktail_simplex,
  units = "cm"
)

ggtern::ggsave(
  filename = paste0(out_folder, "res_cocktail_design_ib_simplex.png"),
  plot = opdesmixr:::plot_choice_set_utility(
    cocktail_I_opt,
    utilities_cocktail_bayesian_plot,
    utility_point_size = utility_point_size_cocktail,
    utility_point_shape = utility_point_shape_cocktail,
    low_color_gradient = low_color_gradient,
    high_color_gradient = high_color_gradient
    ) +
    # ggtitle("Bayesian I-optimal for cocktail experiment") +
    theme(legend.position = "none"),
  width = width_cocktail_simplex,
  height = height_cocktail_simplex,
  units = "cm"
)



# ggtern::grid.arrange(
#   opdesmixr:::plot_choice_set_utility(cocktail_D_opt, utilities_cocktail_bayesian_plot) +
#     ggtitle("Bayesian D-optimal for cocktail experiment")
#   ,
#   opdesmixr:::plot_choice_set_utility(cocktail_I_opt, utilities_cocktail_bayesian_plot) +
#     ggtitle("Bayesian I-optimal for cocktail experiment")
#   ,
#   ncol = 2
# )






### FDS plots

n_points_per_alternative_cocktail = 2000


pred_vars_cocktail = mnl_get_fds_simulations(
  design_array = cocktail_I_opt$X,
  beta = cocktail_I_opt$beta,
  order = 3,
  n_points_per_alternative = n_points_per_alternative_cocktail,
  transform_beta = F,
  verbose = 1) %>%
  mutate(Design = "I-optimal") %>%
  bind_rows(
    mnl_get_fds_simulations(
      design_array = cocktail_D_opt$X,
      beta = cocktail_D_opt$beta,
      order = 3,
      n_points_per_alternative = n_points_per_alternative_cocktail,
      transform_beta = F,
      verbose = 1) %>%
      mutate(Design = "D-optimal")
  )



cocktail_fds_plots_untransf = pred_vars_cocktail %>%
  ggplot() +
  geom_vline(xintercept = 0.5, linetype = "dashed", size = 0.2) +
  geom_hline(yintercept = pred_vars_cocktail %>%
               group_by(Design) %>%
               summarize(
                 med = median(pred_var),
                 mean = mean(pred_var)) %>%
               pull(med),
             linetype = "dashed", size = 0.2) +
  geom_line(aes(fraction, pred_var, linetype = Design), size = 0.8) +
  xlab("Fraction of design space") +
  ylab("Prediction variance") +
  ggtitle("Cocktail experiment") +
  theme_bw() +
  theme(legend.position = "right")


width_cocktail_fds = 18
height_cocktail_fds = 10

ggplot2::ggsave(
  filename = paste0(out_folder, "res_cocktail_fds_db_vs_ib_plot.png"),
  plot = cocktail_fds_plots_untransf,
  width = width_cocktail_fds,
  height = height_cocktail_fds,
  units = "cm"
)





# mnl_get_fds_simulations_all = function(design_array, beta, order, n_points_per_alternative = 500, transform_beta = F, verbose = 0){
#
#   q = dim(design_array)[1]
#   J = dim(design_array)[2]
#   S = dim(design_array)[3]
#
#   pred_var = suppressWarnings(as.data.frame(matrix(rep(NA_real_, nrow(beta)*J*n_points_per_alternative), ncol = J))) %>%
#     purrr::set_names(paste0("V", 1:J)) %>%
#     mutate(beta_ix = rep(1:nrow(beta), n_points_per_alternative))
#
#   progress_old = -1
#   for(k in 1:n_points_per_alternative){
#     if(verbose > 0){
#       progress = round(100*(k-1)/n_points_per_alternative, 0)
#       if(progress - progress_old >= 1){
#         cat('\r', "Progress: ", progress, "%", sep = "")
#         flush.console()
#       }
#       progress_old = progress
#     }
#     des_k = mnl_create_random_initial_design(q, J, S, seed = k)
#     vars_1 = matrix(rep(NA_real_, J*nrow(beta)), ncol = J)
#     for(j in 1:J){
#       f_x = mnl_get_Xs(des_k, 1, order = order)[j,]
#       acc = 0
#       for(i in 1:nrow(beta)){
#         inf_mat = mnl_get_information_matrix(design_array, beta = beta[i,], order = order, transform_beta = transform_beta)
#         vars_1[i,j] = t(f_x) %*% solve(inf_mat, f_x)
#       }
#
#     }
#
#     pred_var[((k-1)*nrow(beta) + 1):(k*nrow(beta)),1:J] = vars_1
#   }
#
#   # out = tibble(pred_var = sort(unlist(pred_var))) %>%
#   # mutate(fraction = 1:nrow(.)/nrow(.))
#
#   out = pred_var %>%
#     pivot_longer(cols = -beta_ix, values_to = "pred_var") %>%
#     select(-name) %>%
#     group_by(beta_ix) %>%
#     # mutate(ix = 1:n()) %>%
#     # ungroup() %>%
#     # group_by(ix) %>%
#     mutate(fraction = (0:(n()-1))/n()) %>%
#     mutate(pred_var = quantile(pred_var, probs = fraction, type = 1)) %>%
#     ungroup()
#
#   # pred_var = unlist(pred_var)
#   #
#   # out = dplyr::tibble(
#   #   fraction = (0:length(pred_var))/length(pred_var)
#   # ) %>%
#   #   dplyr::mutate(pred_var = quantile(pred_var, probs = fraction, type = 1))
#
#
#   if(verbose > 0) cat("\nFinished\n\n")
#
#   return(out)
# }
#
#
# fds_test = mnl_get_fds_simulations_all(
#   design_array = cocktail_I_opt$X,
#   beta = cocktail_I_opt$beta,
#   order = 3,
#   n_points_per_alternative = 500,
#   transform_beta = F,
#   verbose = 1)
#
# fds_test %>%
#   ggplot() +
#   geom_line(aes(fraction, pred_var, group = beta_ix), size = 0.2) +
#   xlab("Fraction of design space") +
#   ylab("Prediction variance") +
#   ggtitle("Cocktail experiment") +
#   theme_bw() +
#   theme(legend.position = "right")
#
# fds_test %>%
#   ggplot() +
#   geom_line(aes(fraction, pred_var, group = beta_ix), size = 0.2) +
#   xlab("Fraction of design space") +
#   ylab("Prediction variance") +
#   ggtitle("Cocktail experiment") +
#   theme_bw() +
#   theme(legend.position = "right") +
#   ylim(0, 4)
#
# fds_test %>%
#   group_by(beta_ix) %>%
#   mutate(ix = 1:n()) %>%
#   ungroup() %>%
#   group_by(ix, fraction) %>%
#   summarize(pred_var = mean(pred_var)) %>%
#   ungroup() %>%
#   ggplot() +
#   geom_line(aes(fraction, pred_var), size = 0.8) +
#   xlab("Fraction of design space") +
#   ylab("Prediction variance") +
#   ggtitle("Cocktail experiment") +
#   theme_bw() +
#   theme(legend.position = "right") +
#   ylim(0, 4)
#
# aaa = mnl_get_fds_simulations(
#   design_array = cocktail_I_opt$X,
#   beta = cocktail_I_opt$beta,
#   order = 3,
#   n_points_per_alternative = 500,
#   transform_beta = F,
#   verbose = 1)
#
# aaa %>%
#   ggplot() +
#   geom_line(aes(fraction, pred_var), size = 0.8) +
#   xlab("Fraction of design space") +
#   ylab("Prediction variance") +
#   ggtitle("Cocktail experiment") +
#   theme_bw() +
#   theme(legend.position = "right") +
#   ylim(0, 4)

##########################################################################################
#### Plot designs for Cornell's experiment
##########################################################################################

# Parameters for saving plots in PNG
width_cornell_simplex = 10
height_cornell_simplex = 10
utility_point_size_cornell = 0.01
utility_point_shape_cornell = "circle"

width_cornell_fds = 20
height_cornell_fds = 12

# Number of sampled points for FDS plot
fds_n_points_per_alternative_cornell = 2000

# Create utility contours
levels_cornell_bayesian = c(0, 0.375, 0.75, 1.125, 1.34)
utilities_cornell_bayesian_plot = opdesmixr:::scheffe_order3_q3_utilities(beta_2_prime, 100000) %>%
  mutate(utility_fact = cut(utility, levels_cornell_bayesian, right = F)) %>%
  mutate(utility_int = as.integer(utility_fact)) %>%
  mutate(utility_avg = case_when(
    utility_int == 1 ~ (levels_cornell_bayesian[1] + levels_cornell_bayesian[2] - levels_cornell_bayesian[1])/2,
    utility_int == 2 ~ (levels_cornell_bayesian[2] + levels_cornell_bayesian[3] - levels_cornell_bayesian[2])/2,
    utility_int == 3 ~ (levels_cornell_bayesian[3] + levels_cornell_bayesian[4] - levels_cornell_bayesian[3])/2,
    utility_int == 4 ~ (levels_cornell_bayesian[4] + levels_cornell_bayesian[5] - levels_cornell_bayesian[4])/2
  ))





##### Transformed betas


## Save PNGs of designs

for(i in seq_along(cornell_designs_transf)){

  kappa = cornell_designs_transf[[i]]$kappa
  cat("Doing kappa =", kappa, "\n")

  d_opt_plot_transf_i = opdesmixr:::plot_choice_set_utility(
    cornell_designs_transf[[i]]$d_opt$X, utilities_cornell_bayesian_plot,
    utility_point_size = utility_point_size_cornell,
    utility_point_shape = utility_point_shape_cornell,
    legend.position = "none",
    low_color_gradient = low_color_gradient,
    high_color_gradient = high_color_gradient
  ) +
    theme(plot.margin = unit(c(-5000000, -0.5, -5000000, -0.5), "cm")) # top, right, bottom, and left margins

  i_opt_plot_transf_i = opdesmixr:::plot_choice_set_utility(
    cornell_designs_transf[[i]]$i_opt$X, utilities_cornell_bayesian_plot,
    utility_point_size = utility_point_size_cornell,
    utility_point_shape = utility_point_shape_cornell,
    legend.position = "none",
    low_color_gradient = low_color_gradient,
    high_color_gradient = high_color_gradient
  ) +
    theme(plot.margin = unit(c(-5000000, -0.5, -5000000, -0.5), "cm")) # top, right, bottom, and left margins


  plot_filename_basename = paste0("res_cornell_transf_simplex_kappa_", sprintf("%03d", kappa*10), "_")

  ggtern::ggsave(
    filename = paste0(out_folder, plot_filename_basename, "D_opt.png"),
    plot = d_opt_plot_transf_i,
    width = width_cornell_simplex,
    height = height_cornell_simplex,
    units = "cm"
  )

  ggtern::ggsave(
    filename = paste0(out_folder, plot_filename_basename, "I_opt.png"),
    plot = i_opt_plot_transf_i,
    width = width_cornell_simplex,
    height = height_cornell_simplex,
    units = "cm"
  )
}


##### FDS plots

pred_vars_cornell_transf = lapply(
  seq_along(cornell_designs_transf), function(i){

    kappa_i = cornell_designs_transf[[i]]$kappa
    cat("\n\n\nkappa:", kappa_i, format(Sys.time(),'(%H:%M:%S)'), "\n")

    pred_vars_cornell_transf_i = mnl_get_fds_simulations(
      design_array = cornell_designs_transf[[i]]$i_opt$X,
      beta = cornell_designs_transf[[i]]$i_opt$beta,
      order = 3,
      n_points_per_alternative = fds_n_points_per_alternative_cornell,
      transform_beta = T,
      verbose = 1) %>%
      mutate(Design = "I-optimal") %>%
      bind_rows(
        mnl_get_fds_simulations(
          design_array = cornell_designs_transf[[i]]$d_opt$X,
          beta = cornell_designs_transf[[i]]$d_opt$beta,
          order = 3,
          n_points_per_alternative = fds_n_points_per_alternative_cornell,
          transform_beta = T,
          verbose = 1) %>%
          mutate(Design = "D-optimal")
      )

    return(pred_vars_cornell_transf_i %>%
             mutate(kappa = kappa_i))

  }) %>%
  bind_rows()



cornell_fds_plots_transf = pred_vars_cornell_transf %>%
  left_join(
    pred_vars_cornell_transf %>%
      group_by(Design, kappa) %>%
      summarize(
        med = median(pred_var),
        mean = mean(pred_var))
  ) %>%
  mutate(kappa2 = paste0("kappa = ", kappa)) %>%
  mutate(kappa2 = fct_reorder(kappa2, kappa)) %>%
  ggplot() +
  geom_vline(xintercept = 0.5, linetype = "dashed", size = 0.2) +
  geom_hline(aes(yintercept = med), linetype = "dashed", size = 0.2) +
  geom_line(aes(fraction, pred_var, linetype = Design), size = 0.8) +
  xlab("Fraction of design space") +
  ylab("Prediction variance") +
  ggtitle("Cornell's experiment") +
  theme_bw() +
  theme(legend.position = "bottom") +
  facet_wrap(~kappa2, scales = "free_y")





ggplot2::ggsave(
  filename = paste0(out_folder, "res_cornell_transf_fds_db_vs_ib_plot.png"),
  plot = cornell_fds_plots_transf,
  width = width_cornell_fds,
  height = height_cornell_fds,
  units = "cm"
)








##### Untransformed betas


## Save PNGs of designs

for(i in seq_along(cornell_designs_untransf)){

  kappa = cornell_designs_untransf[[i]]$kappa
  cat("Doing kappa =", kappa, "\n")

  d_opt_plot_untransf_i = opdesmixr:::plot_choice_set_utility(
    cornell_designs_untransf[[i]]$d_opt$X, utilities_cornell_bayesian_plot,
    utility_point_size = utility_point_size_cornell,
    utility_point_shape = utility_point_shape_cornell,
    legend.position = "none",
    low_color_gradient = low_color_gradient,
    high_color_gradient = high_color_gradient
  ) +
    theme(plot.margin = unit(c(-5000000, -0.5, -5000000, -0.5), "cm")) # top, right, bottom, and left margins

  i_opt_plot_untransf_i = opdesmixr:::plot_choice_set_utility(
    cornell_designs_untransf[[i]]$i_opt$X, utilities_cornell_bayesian_plot,
    utility_point_size = utility_point_size_cornell,
    utility_point_shape = utility_point_shape_cornell,
    legend.position = "none",
    low_color_gradient = low_color_gradient,
    high_color_gradient = high_color_gradient
  ) +
    theme(plot.margin = unit(c(-5000000, -0.5, -5000000, -0.5), "cm")) # top, right, bottom, and left margins


  plot_filename_basename = paste0("res_cornell_untransf_simplex_kappa_", sprintf("%03d", kappa*10), "_")

  ggtern::ggsave(
    filename = paste0(out_folder, plot_filename_basename, "D_opt.png"),
    plot = d_opt_plot_untransf_i,
    width = width_cornell_simplex,
    height = height_cornell_simplex,
    units = "cm"
  )

  ggtern::ggsave(
    filename = paste0(out_folder, plot_filename_basename, "I_opt.png"),
    plot = i_opt_plot_untransf_i,
    width = width_cornell_simplex,
    height = height_cornell_simplex,
    units = "cm"
  )

}


##### FDS plots

pred_vars_cornell_untransf = lapply(
  seq_along(cornell_designs_untransf), function(i){

    kappa_i = cornell_designs_untransf[[i]]$kappa
    cat("\n\n\nkappa:", kappa_i, format(Sys.time(),'(%H:%M:%S)'), "\n")

    pred_vars_cornell_untransf_i = mnl_get_fds_simulations(
      design_array = cornell_designs_untransf[[i]]$i_opt$X,
      beta = cornell_designs_untransf[[i]]$i_opt$beta,
      order = 3,
      n_points_per_alternative = fds_n_points_per_alternative_cornell,
      transform_beta = F,
      verbose = 1) %>%
      mutate(Design = "I-optimal") %>%
      bind_rows(
        mnl_get_fds_simulations(
          design_array = cornell_designs_untransf[[i]]$d_opt$X,
          beta = cornell_designs_untransf[[i]]$d_opt$beta,
          order = 3,
          n_points_per_alternative = fds_n_points_per_alternative_cornell,
          transform_beta = F,
          verbose = 1) %>%
          mutate(Design = "D-optimal")
      )

    return(pred_vars_cornell_untransf_i %>%
             mutate(kappa = kappa_i))

  }) %>%
  bind_rows()



cornell_fds_plots_untransf = pred_vars_cornell_untransf %>%
  left_join(
    pred_vars_cornell_untransf %>%
      group_by(Design, kappa) %>%
      summarize(
        med = median(pred_var),
        mean = mean(pred_var))
  ) %>%
  mutate(kappa2 = paste0("kappa = ", kappa)) %>%
  mutate(kappa2 = fct_reorder(kappa2, kappa)) %>%
  ggplot() +
  geom_vline(xintercept = 0.5, linetype = "dashed", size = 0.2) +
  geom_hline(aes(yintercept = med), linetype = "dashed", size = 0.2) +
  geom_line(aes(fraction, pred_var, linetype = Design), size = 0.8) +
  xlab("Fraction of design space") +
  ylab("Prediction variance") +
  ggtitle("Cornell's experiment") +
  theme_bw() +
  theme(legend.position = "bottom") +
  facet_wrap(~kappa2, scales = "free_y")





ggplot2::ggsave(
  filename = paste0(out_folder, "res_cornell_untransf_fds_db_vs_ib_plot.png"),
  plot = cornell_fds_plots_untransf,
  width = width_cornell_fds,
  height = height_cornell_fds,
  units = "cm"
)



# print(cornell_fds_plots_untransf)

#
# ## Create a list with the plots of the optimal designs
# cornell_optimal_designs_plots_untransf = lapply(seq_along(cornell_designs_untransf), function(i){
#
#   out_list = list(
#     d_opt_plot = opdesmixr:::plot_choice_set_utility(
#       cornell_designs_untransf[[i]]$d_opt$X, utilities_cornell_bayesian_plot,
#       legend.position = "none"
#     ) +
#       ggtitle("", subtitle = paste0("D-optimal, kappa = ", cornell_designs_untransf[[i]]$kappa))
#     ,
#     i_opt_plot = opdesmixr:::plot_choice_set_utility(
#       cornell_designs_untransf[[i]]$i_opt$X, utilities_cornell_bayesian_plot,
#       legend.position = "none"
#     ) +
#       ggtitle("", subtitle = paste0("I-optimal, kappa = ", cornell_designs_untransf[[i]]$kappa))
#
#     ,
#     kappa = cornell_designs_untransf[[i]]$kappa
#   )
#
#   return(out_list)
#
#
# })
# # Couldn't get it to print with the list only
#
# ggtern::grid.arrange(
#   cornell_optimal_designs_plots_untransf[[1]][[1]], cornell_optimal_designs_plots_untransf[[1]][[2]],
#   cornell_optimal_designs_plots_untransf[[2]][[1]], cornell_optimal_designs_plots_untransf[[2]][[2]],
#   cornell_optimal_designs_plots_untransf[[3]][[1]], cornell_optimal_designs_plots_untransf[[3]][[2]],
#   cornell_optimal_designs_plots_untransf[[4]][[1]], cornell_optimal_designs_plots_untransf[[4]][[2]],
#   top = grid::textGrob(paste0("kappa = 0.5, 5, 10, 30"), gp = grid::gpar(fontsize = 15)),
#   ncol = 2)

