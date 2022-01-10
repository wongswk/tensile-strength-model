# Simulation study:  simulate data and fit Bayesian model via Stan
# - set sample size N on  line 24

# load the data simulation function
rm(list = ls())
source("stan_sim.R")
source("results_print_with_moe.R")

library(rstan)
library(ggmcmc)
library(HDInterval)

# valid test span
grip = 24 # machine grip in inches 
# setup x range
x_start = grip
x_end = 144 - grip
# setup y range
y_start =  0 
y_end =  5.5

# Partition number 
J = 24L
# Sample size 
N = 120
# maximum effective distance 
d_max = 96

# Parameters
eta0 = 3
eta1 = 1.5
rho = 0.7
sigma = 0.8
beta = 0.5
gamma0 = 0.25
gamma1 = 0.15
lambda = 0.01
seed <- 12345

set.seed(seed)

sim_data_list  = simDataSTAN(eta0 = eta0,
                        eta1 = eta1,
                        rho = rho,
                        sigma = sigma,
                        beta = beta, 
                        gamma0 = gamma0,
                        gamma1 = gamma1,
                        N = N, 
                        J = J, 
                        d_max = d_max, 
                        lambda = lambda, 
                        x_start = x_start, 
                        x_end = x_end, 
                        y_start = y_start, 
                        y_end = y_end)

MOE_vec = as.array(sim_data_list[['MOE_vec']])
K_vec = as.array(sim_data_list[['K_vec']])  # number of knots for each lumber
Kmax = max(K_vec)                           # maximum number of knots among all samples
list_Z = sim_data_list[['list_Z']]          # list of knots effects - percentage 
list_D = sim_data_list[['list_D']]          # list of distance matrices
list_E = sim_data_list[['list_E']]          # list of the edge knot indicator
Ym = sim_data_list[['Ym']]                  # 'observed' minimum.
min_index = sim_data_list[['min_index']]    # 'observed' indices
Y_matrix = sim_data_list[['Y_matrix']]      # full simulated data matrix
scaled_weighted_effects = sim_data_list[['scaled_weighted_effects']] # total knots effects weighted by the distance matrix for per cell per lumber

list_Z_uniform = function(N, Z_list){
  Z_list_uniform = list()
  for (i in 1:N) {
    Z_i = Z_list[[i]]
    Z_new = matrix(0, nrow = Kmax, ncol = 1)
    Z_new[1:length(Z_i),] = Z_i
    Z_list_uniform[[i]] = Z_new
  }
  return(Z_list_uniform)
}

list_Edge_uniform = function(N, Edge_list){
  Edge_list_uniform = list()
  for (i in 1:N) {
    E_i = Edge_list[[i]]
    E_new = matrix(-1, nrow = Kmax, ncol = 1)
    E_new[1:length(E_i),] = E_i
    Edge_list_uniform[[i]] = E_new
  }
  return(Edge_list_uniform)
}


list_D_uniform = function(N,D_list){
  D_list_uniform = list()
  for (i in 1:N) {
    D_i = D_list[[i]]
    D_new = matrix(1e4, nrow = J, ncol = Kmax)
    D_new[1:nrow(D_i), 1:ncol(D_i)] = D_i
    D_list_uniform[[i]]= D_new
  }
  return(D_list_uniform)
}


Z_list_uniform = list_Z_uniform(N = N, Z_list = list_Z)
Z_array = array(as.numeric(unlist(Z_list_uniform)), dim=c(Kmax, 1,N))
D_list_uniform = list_D_uniform(N = N, D_list = list_D)
D_array = array(as.numeric(unlist(D_list_uniform)), dim=c(J, Kmax,N))
Edge_list_uniform = list_Edge_uniform(N = N, Edge_list = list_E)
Edge_array = array(as.numeric(unlist(Edge_list_uniform)), dim=c(Kmax, 1,N))


library(rstan)
rstan_options(auto_write = TRUE)             # complied model saved in same directory as stan file
options(mc.cores = parallel::detectCores())  # run parallel::detectCores() chains in parallel

data = list(
  N = N,
  J = J,
  Kmax = Kmax,
  y_obs = Ym,
  min_index = min_index,
  D_array = D_array,
  Z_array = Z_array,
  Edge_array = Edge_array,
  dmax = d_max,
  K_vec = K_vec,
  MOE_vec = MOE_vec)

saveRDS(sim_data_list, file = paste0("sim_data_list_N", N, ".rds"))

stan_fit_comp = stan("real_data_moe.stan", 
                     data = data,
                     iter = 10000,
                     warmup = 5000,
                     chain = 4,
                     pars=c("eta0", "eta1", "rho", "sigma", "beta", "gamma0", "gamma1"),
                     refresh=5,
                     control = list(adapt_delta = 0.99, max_treedepth = 15))

#print(stan_fit_comp, digits = 5)
#result_print(stan_fit_comp)

saveRDS(stan_fit_comp, file = paste0("stan_fit_sim_N", N, ".rds"))





