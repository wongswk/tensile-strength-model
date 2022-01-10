## Fit model to real data

rm(list = ls())

# source the file that process the raw data 
source("data_process.R")
source("results_print_with_moe.R")

## Run with the real data
library(rstan)
rstan_options(auto_write = TRUE)             # complied model saved in same directory as stan file
options(mc.cores = parallel::detectCores())  # run parallel::detectCores() chains in parallel

data = list(
  N = N,
  J = J,
  Kmax = Kmax,
  y_obs = Ym_ksi,
  min_index = min_index,
  D_array = D_array,
  Z_array = Z_array,
  Edge_array = Edge_array,
  dmax = d_max,
  K_vec = K_vec,
  MOE_vec = MOE_vec)

stan_fit_comp = stan("real_data_moe.stan", 
                     data = data,
                     iter = 10000,
                     warmup = 5000,
                     chain = 4,
                     pars=c("eta0", "eta1", "rho", "sigma", "beta", "gamma0", "gamma1"),
                     refresh=5,
                     control = list(adapt_delta = 0.99, max_treedepth = 15))

saveRDS(stan_fit_comp, file="realdata-fit.rds")
save(Ym_ksi, J,d_max,N,MOE_vec, K_vec, D_list, Z_list, Edge_list, D_list_uniform, Z_list_uniform, Edge_list_uniform, file="realdata-input.rda")
