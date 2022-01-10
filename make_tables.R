## Summarize simulation study parameter estimates (Table 1)
library(rstan)
N <- c(120,360,720)

simlist <- list()
parslist <- list()
restable <- list()
for (i in 1:length(N)) {
  simlist[[i]] <- readRDS(file = paste0("stan_fit_sim_N", N[i], ".rds"))
  
  #print(simlist[[i]] ,  digits = 2)
  parslist[[i]] <- extract(simlist[[i]])
  restable[[i]] <- t(round(sapply(parslist[[i]], function(x) quantile(x,c(0.5,0.025,0.975))),2))  
}

truepar <- c(3,1.5,0.7,0.8,0.5, 0.25,0.15, NA)

library(xtable)
print(xtable(cbind(truepar, restable[[1]], restable[[2]], restable[[3]])))


## Make MOE, UTS, cell index histograms (Figure 3)
load("realdata-input.rda")
pdf("hist-MOE-UTS-index.pdf", width=8, height=3.5)
#par(mfrow=c(1,3), mar=c(5,4,1,0))
layout(matrix(c(1,2,3,3), nrow = 1, ncol = 4, byrow = TRUE))
par(mar=c(5,4,1,0))
hist(MOE_vec, main="", xlab=expression(bold("MOE (psi "%*%~10^6~")")))
hist(Ym_ksi, main="", xlab=expression(bold("UTS (psi "%*%~10^3~")")))
hist(min_index, breaks=seq(0.5,24.5,by=1), xlim=c(0,25), main="", xlab=expression(bold("Cell index")))
dev.off()

## Summarize real data parameter estimates (Table 2)
stan_fit_real <- readRDS("realdata-fit.rds")
realparslist <- extract(stan_fit_real)
realrestable <- t(round(sapply(realparslist, function(x) c(mean(x),quantile(x,c(0.5,0.025,0.975)))),4))  
print(xtable(realrestable))


     