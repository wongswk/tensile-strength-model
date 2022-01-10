## Helper functions for computing predictive distributions

get_Wlist = function(N, beta, d_max, list_D){
  list_W = list()
  for (i in 1:N) {
    D = list_D[[i]]
    W = exp(-beta*D)*as.numeric(D<= d_max)
    list_W[[i]] = W
  }
  return(list_W)
}

get_list_of_Wlist = function(beta_post, N_test, N_samples, d_max, list_D){
  list_of_Wlist = list()
  for(n in 1:N_samples){
    list_of_Wlist[[n]] = get_Wlist(N = N_test, beta = beta_post[n],
                                   d_max = d_max, list_D = list_D)
  }
  return(list_of_Wlist)
}

get_scaled_Zlist = function(gamma0, gamma1,list_Z, list_E){
  N = length(list_Z)
  list_scaled_Z = list()
  for (i in 1:N) {
    Z_i = list_Z[[i]]
    K_vec_i = nrow(Z_i)
    E_i = list_E[[i]]
    gamma_i = diag(gamma1, nrow = K_vec_i, ncol = K_vec_i)%*%E_i + 
      diag(gamma0, nrow = K_vec_i, ncol = K_vec_i)%*%(1-E_i)
    scaled_Z_i = gamma_i*Z_i
    list_scaled_Z[[i]] = scaled_Z_i
  }
  return(list_scaled_Z)
}  

get_list_of_scaled_Zlist = function(gamma0_post, gamma1_post, N_test, 
                                    N_samples,list_Z, list_E){
  list_of_scaled_Zlist = list()
  for (n in 1:N_samples) {
    list_of_scaled_Zlist[[n]] = get_scaled_Zlist(gamma0 = gamma0_post[n],
                                                 gamma1 = gamma1_post[n],
                                                 list_Z = list_Z,
                                                 list_E = list_E)
  }
  return(list_of_scaled_Zlist) 
}


get_weighted_effects_list = function(beta_post,
                                     gamma0_post,
                                     gamma1_post,
                                     N_test, 
                                     N_samples,
                                     d_max, 
                                     list_D, 
                                     list_E, 
                                     list_Z){
  list_of_Wlist = get_list_of_Wlist(beta_post = beta_post, 
                                    N_test = N_test,
                                    N_samples = N_samples,
                                    d_max = d_max,
                                    list_D = list_D)
  
  list_of_scaled_Zlist = get_list_of_scaled_Zlist(gamma0_post = gamma0_post, 
                                                  gamma1_post = gamma1_post, 
                                                  N_test = N_test, 
                                                  N_samples = N_samples,
                                                  list_Z = list_Z, 
                                                  list_E = list_E)
  
  list_of_weighted_effects_list = list() # list of weighted effects list 
  for(i in 1:N_test){
    weighted_effects_list = list()
    for(n in 1:N_samples){
      scaled_Zlist = list_of_scaled_Zlist[[n]]
      W_list = list_of_Wlist[[n]]
      weighted_effects_list[[n]] = W_list[[i]] %*% scaled_Zlist[[i]]
    }
    list_of_weighted_effects_list[[i]] = weighted_effects_list
  }
  return(list_of_weighted_effects_list)
}


# Function to simualte Y based on new data 

get_pred_quant = function(J,
                          d_max,
                          N_test,
                          MOE_vec_test,
                          K_vec_test,
                          list_D_test,
                          list_Z_test,
                          list_E_test,
                          stan_fit_obj){
  
  post_spls = as.data.frame(stan_fit_obj)
  eta0_post = post_spls$eta0
  eta1_post = post_spls$eta1
  rho_post = post_spls$rho
  sigma_post = post_spls$sigma
  beta_post = post_spls$beta
  gamma0_post = post_spls$gamma0
  gamma1_post = post_spls$gamma1
  N_samples = dim(post_spls)[1]
  
  list_of_weighted_effects_list = get_weighted_effects_list(beta_post = beta_post,
                                                            gamma0_post = gamma0_post,
                                                            gamma1_post = gamma1_post,
                                                            N_test = N_test,
                                                            N_samples = N_samples,
                                                            d_max = d_max,
                                                            list_D = list_D_test, 
                                                            list_E = list_E_test, 
                                                            list_Z = list_Z_test)
  list_Y_all = list()
  list_Y_obs = list()
  for (i in 1:N_test) {
    show(i)
    K_i = K_vec_test[i]
    MOE_i = MOE_vec_test[i]
    Y_all = matrix(NA, nrow = N_samples, ncol = J)
    Y_obs = matrix(NA, nrow = N_samples, ncol = 1)
    for(n in 1:N_samples){
      #Y_all[n,1] <- -1
      #while(Y_all[n,1] < 0) {
        Y_all[n, 1] = rnorm(mean = eta0_post[n] + eta1_post[n]*MOE_i - list_of_weighted_effects_list[[i]][[n]][1],
                          sd = sigma_post[n]/sqrt(1-rho_post[n]^2), n = 1)
      #}
      for (j in 2:J){
        #Y_all[n,j] <- -1
        #while(Y_all[n,j] < 0) {
          Y_all[n, j] = rnorm(mean = (1-rho_post[n])*(eta0_post[n] + eta1_post[n]*MOE_i)+
                              rho_post[n]*Y_all[n,j-1] - 
                              list_of_weighted_effects_list[[i]][[n]][j] + 
                              rho_post[n]*list_of_weighted_effects_list[[i]][[n]][j-1],
                            sd = sigma_post[n], n = 1)
        #}
      }
      Y_obs[n, 1] = max(min(Y_all[n,]),0)
    }
    list_Y_all[[i]] = Y_all
    list_Y_obs[[i]] = Y_obs
  }
  result = list(list_Y_all = list_Y_all,
                list_Y_obs = list_Y_obs)
  return(result)
}
