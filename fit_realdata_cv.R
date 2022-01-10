## Fit model to real data: perform K-fold cross-validation

rm(list = ls())

# source the file that process the raw data 
source("data_process.R")
source("results_print_with_moe.R")

## Run with the real data
library(rstan)
rstan_options(auto_write = TRUE)             # complied model saved in same directory as stan file
options(mc.cores = parallel::detectCores())  # run parallel::detectCores() chains in parallel

set.seed(123)

# K-fold CV
K <- 5
validSetSplits <- sample((1:N)%%K + 1)
stan_fit_cv <- list()

for (k in 1:K) {
  show(k)

  data = list(
    N = sum(validSetSplits!=k),
    J = J,
    Kmax = Kmax,
    y_obs = Ym_ksi[validSetSplits!=k],
    min_index = min_index[validSetSplits!=k],
    D_array = D_array[,,validSetSplits!=k,drop=FALSE],
    Z_array = Z_array[,,validSetSplits!=k,drop=FALSE],
    Edge_array = Edge_array[,,validSetSplits!=k,drop=FALSE],
    dmax = d_max,
    K_vec = K_vec[validSetSplits!=k],
    MOE_vec = MOE_vec[validSetSplits!=k])
  
  stan_fit_cv[[k]] = stan("real_data_moe.stan", 
                       data = data,
                       iter = 10000,
                       warmup = 5000,
                       chain = 4,
                       pars=c("eta0", "eta1", "rho", "sigma", "beta", "gamma0", "gamma1"),
                       refresh=5,
                       control = list(adapt_delta = 0.99, max_treedepth = 15))

  save(validSetSplits, stan_fit_cv, file="realdata-cv-fit.rda")
}



#save(Ym_ksi, J,d_max,N,MOE_vec, min_index, K_vec, D_list, Z_list, Edge_list, D_list_uniform, Z_list_uniform, Edge_list_uniform, file="realdata-input.rda")
