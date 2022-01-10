## Summarize cross-validation results

library(rstan)
source("prediction-functions.R")
load("realdata-cv-fit.rda")
load("realdata-input.rda")

source("results_print_with_moe.R")

K <- max(validSetSplits)
postpreds <- list()
for (k in 1:K) {
  show(k)
  result_print(stan_fit_cv[[k]])

  temp = get_pred_quant(J = J, 
                      d_max = d_max, 
                      N_test = sum(validSetSplits==k),
                      MOE_vec_test = MOE_vec[validSetSplits==k],
                      K_vec_test = K_vec[validSetSplits==k],
                      list_D_test = D_list_uniform[validSetSplits==k],
                      list_Z_test = Z_list_uniform[validSetSplits==k],
                      list_E_test = Edge_list_uniform[validSetSplits==k],
                      stan_fit_obj = stan_fit_cv[[k]])
  
  list_Y_all = temp$list_Y_all
  list_Y_obs = temp$list_Y_obs
  
  post_predictions = data.frame(matrix(nrow = sum(validSetSplits==k), ncol = 4))
  colnames(post_predictions) = c('true_val','post_pred','post_lwr','post_upr')
  post_predictions$true_val = Ym_ksi[validSetSplits==k]
  for (i in 1:length(list_Y_obs)) {
    Y_obs_i = as.data.frame(list_Y_obs[[i]])
    Y_obs_i_quantile = quantile(Y_obs_i$V1, probs = c(0.025, 0.975))
    post_predictions$post_pred[i] = mean(Y_obs_i$V1)
    post_predictions$post_lwr[i] = Y_obs_i_quantile[1]
    post_predictions$post_upr[i] = Y_obs_i_quantile[2]
  }
  
  show(mean((post_predictions$true_val - post_predictions$post_pred)^2))
  
  postpreds[[k]] <- post_predictions
  
}

# Metrics
combinedres <- rbind(postpreds[[1]],postpreds[[2]],postpreds[[3]],postpreds[[4]],postpreds[[5]])
mean(combinedres$post_pred)
show(mean((combinedres$true_val - combinedres$post_pred)^2))
show(mean((combinedres$post_upr - combinedres$post_lwr)))
show(mean(abs(combinedres$true_val - combinedres$post_pred)))


## Comparisons with basic regression models
regpreds <- list()
max_knot_size <- sapply(Z_list_uniform, max)
for (k in 1:K) {
  show(k)

  data_train = data.frame(
    y_obs = Ym_ksi[validSetSplits!=k],
    max_knot_size = max_knot_size[validSetSplits!=k],
    MOE_vec = MOE_vec[validSetSplits!=k])  
  
  data_pred = data.frame(y_obs = Ym_ksi[validSetSplits==k], 
                          max_knot_size = max_knot_size[validSetSplits==k],
                          MOE_vec = MOE_vec[validSetSplits==k])
  
  lm1 <- lm(y_obs~MOE_vec, data=data_train)
  lm2 <- lm(y_obs~MOE_vec+max_knot_size, data=data_train)

  data_pred$reg_pred1 <- predict(lm1, newdata = data_pred, interval="predict")
  data_pred$reg_pred2 <- predict(lm2, newdata = data_pred, interval="predict")

  regpreds[[k]] <- data_pred 
}

# MSPE
combinedreg <- rbind(regpreds[[1]],regpreds[[2]],regpreds[[3]],regpreds[[4]],regpreds[[5]])
show(c(mean((combinedreg$y_obs - combinedreg$reg_pred1[,1])^2), mean(combinedreg$reg_pred1[,3] - combinedreg$reg_pred1[,2])))
show(c(mean((combinedreg$y_obs - combinedreg$reg_pred2[,1])^2), mean(combinedreg$reg_pred2[,3] - combinedreg$reg_pred2[,2])))
