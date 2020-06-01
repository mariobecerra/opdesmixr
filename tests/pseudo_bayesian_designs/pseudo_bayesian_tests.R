library(tidyverse)
library(here)

# library(opdesmixr)
Rcpp::sourceCpp("src/utils.cpp")
source("R/mnl_model_functions.R")
source("R/general_functions.R")



out_dir_1 = "tests/pseudo_bayesian_designs/out/"
out_dir = paste0(out_dir_1, "rds_files/")
dir.create(here(out_dir_1))
dir.create(here(out_dir))


theme_set(theme_bw())


