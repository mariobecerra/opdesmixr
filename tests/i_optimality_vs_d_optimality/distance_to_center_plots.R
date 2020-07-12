library(tidyverse)
library(stringr)
library(here)

theme_set(theme_bw())


compute_average_distance_to_center = function(design_array){

  dim_X = dim(design_array)
  q = dim_X[1]
  S = dim_X[3]

  # Convert 3 dimensional arrays into matrices by vertically binding them
  X_final_mat = t(design_array[,,1])
  for(s in 2:S){
    X_final_mat = rbind(X_final_mat, t(design_array[,,s]))
  }

  center = rep(1/q, q)

  return(mean(sqrt(apply((X_final_mat - center)^2, 1, sum))))
}




out_dir_1 = "tests/i_optimality_vs_d_optimality/out/"
out_dir_2 = paste0(out_dir_1, "rds_files/")
out_dir = paste0(out_dir_2, "many_random_designs/")

file_names = list.files(out_dir)


average_distances = map_df(seq_along(file_names), function(i){

  filename_i = file_names[i]

  q = str_extract(pattern = "[0-9]+", str_extract(pattern = "_.*?_", string = filename_i))
  J = str_extract(pattern = "[0-9]+", str_extract(pattern = "J.*?_", string = filename_i))
  S = str_extract(pattern = "[0-9]+", str_extract(pattern = "S.*?_", string = filename_i))

  designs = readRDS(here(out_dir, filename_i))
  D_opt = designs$D_opt
  I_opt = designs$I_opt

  # If there was an error when computing efficient design
  if(is.na(D_opt$opt_crit_value) |
     is.na(I_opt$opt_crit_value) |
     is.infinite(D_opt$opt_crit_value) |
     is.infinite(I_opt$opt_crit_value)){
    d_eff_average_distance = NA_real_
    i_eff_average_distance = NA_real_
  } else{
    d_eff_average_distance = compute_average_distance_to_center(D_opt$X)
    i_eff_average_distance = compute_average_distance_to_center(I_opt$X)
  }

  out = tibble(
    d_eff_average_distance = d_eff_average_distance,
    i_eff_average_distance = i_eff_average_distance
  ) %>%
    mutate(q = q, J = J, S = S, fname = filename_i)

  return(out)

}) %>%
  mutate(diff = d_eff_average_distance - i_eff_average_distance,
         diff2 = d_eff_average_distance/i_eff_average_distance)





average_distances %>%
  filter(!is.na(diff)) %>%
  mutate(group = paste0("q = ", q, ", J = ", J, ", S = ", S)) %>%
  ggplot() +
  geom_boxplot(aes(group, diff)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  xlab("") +
  ylab(expression(d[D]~"-"~d[I])) +
  coord_flip()



# average_distances %>%
#   filter(!is.na(diff)) %>%
#   mutate(group = paste0("q = ", q, ", J = ", J, ", S = ", S)) %>%
#   ggplot() +
#   geom_boxplot(aes(group, diff2)) +
#   geom_hline(yintercept = 1, linetype = "dotted") +
#   xlab("") +
#   ylab(expression(d[D]~"/"~d[I])) +
#   coord_flip()















