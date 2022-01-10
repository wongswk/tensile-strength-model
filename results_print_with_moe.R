## Helper functions for plotting and summarizing MCMC samples

### Plot distribution for parameters ###
result_print = function(stan_fit_comp, eta0 = NA, eta1 = NA, sigma=NA , rho=NA , gamma0 =NA, gamma1= NA, beta =NA){
  post = as.matrix(stan_fit_comp)
  hdi = apply(post, MARGIN  = 2,FUN = function(x){HDInterval::hdi(x)})
  ncols = ceiling(sum(as.numeric(!is.na(c(eta0, eta1, sigma, rho, gamma0, gamma1, beta))))/2)
  par(mfrow=c(3,3), oma=c(0,0,2,0))
  # eta0
  if(!is.na(eta0)){
    hist(post[,"eta0"], breaks = 20, col = "#808080", border = FALSE,
         main = expression(paste("Distribution of ", eta0)),
         xlab = expression(eta0))
    eta095 <- hdi[,'eta0']
    for (i in 1:2) {
      abline(v = eta095[i], lty = 2, lwd = 2, col = "#40D2FE")
    }
    abline(v = eta0, lty = 2, lwd = 2, col = "red")
  } else {
    hist(post[,"eta0"], breaks = 20, col = "#808080", border = FALSE,
         main = expression(paste("Distribution of ", eta0)),
         xlab = expression(eta0))
    eta095 <- hdi[,'eta0']
    for (i in 1:2) {
      abline(v = eta095[i], lty = 2, lwd = 2, col = "#40D2FE")
    }
  }  
  
  # eta1
  if(!is.na(eta1)){
    hist(post[,"eta1"], breaks = 20, col = "#808080", border = FALSE,
         main = expression(paste("Distribution of ", eta1)),
         xlab = expression(eta1))
    eta195 <- hdi[,'eta1']
    for (i in 1:2) {
      abline(v = eta195[i], lty = 2, lwd = 2, col = "#40D2FE")
    }
    abline(v = eta1, lty = 2, lwd = 2, col = "red")
  } else {
    hist(post[,"eta1"], breaks = 20, col = "#808080", border = FALSE,
         main = expression(paste("Distribution of ", eta1)),
         xlab = expression(eta1))
    eta195 <- hdi[,'eta1']
    for (i in 1:2) {
      abline(v = eta195[i], lty = 2, lwd = 2, col = "#40D2FE")
    }
  }  
  
  if(!is.na(rho)){
    # distribution for rho
    hist(post[,"rho"], breaks = 20, col = "#808080", border = FALSE,
         main = expression(paste("Distribution of ", rho)),
         xlab = expression(rho))
    rho95 <- hdi[,'rho']
    for (i in 1:2) {
      abline(v = rho95[i], lty = 2, lwd = 2, col = "#40D2FE")
    }
    abline(v = rho, lty = 2, lwd = 2, col = "red")
  } else {
    # distribution for rho
    hist(post[,"rho"], breaks = 20, col = "#808080", border = FALSE,
         main = expression(paste("Distribution of ", rho)),
         xlab = expression(rho))
    rho95 <- hdi[,'rho']
    for (i in 1:2) {
      abline(v = rho95[i], lty = 2, lwd = 2, col = "#40D2FE")
    }
  }
  
  if(!is.na(sigma)){
    # distribution for sigma
    hist(post[,"sigma"], breaks = 20, col = "#808080", border = FALSE,
         main = expression(paste("Distribution of ", sigma)),
         xlab = expression(sigma))
    sigma95 <- hdi[,'sigma']
    for (i in 1:2) {
      abline(v = sigma95[i], lty = 2, lwd = 2, col = "#40D2FE")
    }
    abline(v = sigma, lty = 2, lwd = 2, col = "red")
  } else {
    # distribution for sigma
    hist(post[,"sigma"], breaks = 20, col = "#808080", border = FALSE,
         main = expression(paste("Distribution of ", sigma)),
         xlab = expression(sigma))
    sigma95 <- hdi[,'sigma']
    for (i in 1:2) {
      abline(v = sigma95[i], lty = 2, lwd = 2, col = "#40D2FE")
    }
  }
  
  if(!is.na(beta)){
    # distribution for beta
    hist(post[,"beta"], breaks = 20, col = "#808080", border = FALSE,
         main = expression(paste("Distribution of ", beta)),
         xlab = expression(beta))
    beta95 <- hdi[,'beta']
    for (i in 1:2) {
      abline(v = beta95[i], lty = 2, lwd = 2, col = "#40D2FE")
    }
    abline(v = beta, lty = 2, lwd = 2, col = "red")
  } else{
    # distribution for beta
    hist(post[,"beta"], breaks = 20, col = "#808080", border = FALSE,
         main = expression(paste("Distribution of ", beta)),
         xlab = expression(beta))
    beta95 <- hdi[,'beta']
    for (i in 1:2) {
      abline(v = beta95[i], lty = 2, lwd = 2, col = "#40D2FE")
    }
  }
  
  if(!is.na(gamma0)){
    # distribution for gamma0
    hist(post[,"gamma0"], breaks = 20, col = "#808080", border = FALSE,
         main = expression(paste("Distribution of ", gamma0)),
         xlab = expression(gamma0))
    gamma095 <- hdi[,'gamma0']
    for (i in 1:2) {
      abline(v = gamma095[i], lty = 2, lwd = 2, col = "#40D2FE")
    }
    abline(v = gamma0, lty = 2, lwd = 2, col = "red")
  } else {
    # distribution for gamma0
    hist(post[,"gamma0"], breaks = 20, col = "#808080", border = FALSE,
         main = expression(paste("Distribution of ", gamma0)),
         xlab = expression(gamma0))
    gamma095 <- hdi[,'gamma0']
    for (i in 1:2) {
      abline(v = gamma095[i], lty = 2, lwd = 2, col = "#40D2FE")
    }
  }
  # dev.off()
  
  if(!is.na(gamma1)){
    # distribution for gamma1
    hist(post[,"gamma1"], breaks = 20, col = "#808080", border = FALSE,
         main = expression(paste("Distribution of ", gamma1)),
         xlab = expression(gamma1))
    gamma195 <- hdi[,'gamma1']
    for (i in 1:2) {
      abline(v = gamma195[i], lty = 2, lwd = 2, col = "#40D2FE")
    }
    abline(v = gamma1, lty = 2, lwd = 2, col = "red")
  } else {
    # distribution for gamma1
    hist(post[,"gamma1"], breaks = 20, col = "#808080", border = FALSE,
         main = expression(paste("Distribution of ", gamma1)),
         xlab = expression(gamma1))
    gamma195 <- hdi[,'gamma1']
    for (i in 1:2) {
      abline(v = gamma195[i], lty = 2, lwd = 2, col = "#40D2FE")
    }
  }
}