## Households MAIN model code
setwd('/Users/AHenderson/OneDrive - London School of Hygiene and Tropical Medicine/PhD/SISMID/SISMID_MCMC2')

##run functions
source("Code/5households.R")

cases <- c(0,1,2,3)
households <- c(29,9,2,2)
mat <- matrix(NA, nrow=4, 2)
mat[,1] <- cases
mat[,2] <- households
mat <- as.data.frame(mat)
colnames(mat) <- c("No.of.Cases", "No.of.Households")
mat

## COMPUTE log - likelihood
# assume p = 0.2 and alpha = 0.33
# m initial susceptibles
# a initial infectives
# j = number infected at the end in the household...
# alpha infection rate
# the length of the fixed infectious period
p <- 0.2; alpha <- 0.33
## work out q(j,k) where k=3...
for(k in 0:3){
pTk <- compute.marginal.prob(k=k, n=3, c=1, alpha=0.33, p=0.2)
print(paste0("Pr(T=",k,") = ",pTk))
}

# so to sum into a log-likelihood:
loglik <- function(alpha, p){
loglik <- 0
for(k in 0:3){
  loglik <- loglik + mat[k+1,2]*log(compute.marginal.prob(k=k, n=3, c=1, alpha, p))
}
loglik
}
loglik(1, 0.8)

## or sexy smooth coding:
#loglik2 <- function(alpha, p){
# sum(mat[,1] * log(sapply(mat[,2], compute.marginal.prob, n=3, c=1, alpha, p)))
#}
#loglik2(alpha=1, p=0.8)
#

###############################################################
###############################################################
### the MCMC code
mcmc.households <- function(data, iter){
  #init
  p.cur <- init['p']
  alpha.cur <- init['alpha']

  # init priors
  p.prior.cur <- dbeta(p.cur, shape1=p.prior[1],shape2=p.prior[2])
  alpha.prior.cur <- dgamma(alpha.cur, shape=alpha.prior[1], rate=alpha.prior[2])
  
  # init lik
  lik.cur.p <- loglik(alpha=alpha.cur, p=p.cur)
  lik.cur.alpha <- loglik(alpha=alpha.cur, p=p.cur)
  
  p <- rep(NA, iter)
  p[1] <- p.cur
  alpha <- rep(NA, iter)
  alpha[1] <- alpha.cur
  accepttab <- rep(NA, iter)
  
  for(i in 2:iter){
  #proposalm P
  p.prop <- runif(1,0,1)
  
  #likelihoods
  lik.prop.p <- loglik(alpha=alpha.cur, p=p.prop)
  
  #priors
  p.prior.prop <- dbeta(p.prop, shape1=p.prior[1],shape2=p.prior[2])
  
  #MH number
  MH.alg <- (p.prior.prop/p.prior.cur)*   ## prior
              exp(lik.prop.p - lik.cur.p)*    ## lik
              (dunif(p.prop,0,1)/dunif(p.cur,0,1)) ##q ratio
  #selection decision
  if(runif(1) < min(1, MH.alg)){
    lik.cur.p <- lik.prop.p
    p.prior.cur <- p.prior.prop
    p[i] <- p.prop 
    p.cur <- p.prop
    accepttab[i] <- 1
  } else{
    p[i] <- p.cur
    accepttab[i] <- 0
  }
  
  # repeat for ALPHA
  alpha.prop <- rnorm(1,alpha.cur,sd=delta)
  if(alpha.prop>0){
  
  #likelihoods
  lik.prop.alpha <- loglik(alpha=alpha.prop, p=p.cur)
  
  #priors
  alpha.prior.prop <- dgamma(alpha.prop, shape=alpha.prior[1],rate=alpha.prior[2])
  
  #MH number
  MH.alg <- (alpha.prior.prop/alpha.prior.cur)*         ## prior
                exp(lik.prop.alpha - lik.cur.alpha)*    ## lik
                (dnorm(alpha.prop,alpha.cur,delta)/dnorm(alpha.cur,alpha.cur,delta)) ##q ratio
  #selection decision
  if(runif(1) < min(1, MH.alg)){
    lik.cur.alpha <- lik.prop.alpha
    alpha.prior.cur <- alpha.prior.prop
    alpha[i] <- alpha.prop 
    alpha.cur <- alpha.prop
    accepttab[i] <- 1
  } else{
    alpha[i] <- alpha.cur
    accepttab[i] <- 0
  } 
  }else{
    alpha[i] <- alpha.cur
    accepttab[i] <- 0
  }
  }
  return(list(p=p,alpha=alpha, accept=accepttab))
  
}


p.prior <- c(1, 1)
alpha.prior <- c(1, 10^-3)
delta <- 1 # tuning parameter
init <- c(p=0.8, alpha=1)
mcmc.run <- mcmc.households(data=mat, iter=10000)

par(mfrow=c(2,2))
plot(mcmc.run$p, type='l')
hist(mcmc.run$p)
plot(mcmc.run$alpha, type='l')
hist(mcmc.run$alpha)

par(mfrow=c(1,1))
plot(mcmc.run$p, mcmc.run$alpha, col=2)

###############################################################
###############################################################
### Now do the mixing by block sampling from Mutlivariate normal 

source('Code/MCMC_mvtnorm.R')
p.prior <- c(1, 1)
alpha.prior <- c(1, 10^-3)
delta <- 1 # tuning parameter
init <- c(p=0.8, alpha=1)
mcmc.households.mvnorm <- mcmc.households.mvnorm(data=mat, iter=10000)

par(mfrow=c(2,2))
plot(mcmc.households.mvnorm$theta[,1], type='l', main=expression('Trace: p'))
hist(mcmc.households.mvnorm$theta[,1], main=expression('Histogram: p'))
abline(v=mean(mcmc.households.mvnorm$theta[,1]), col=2)
plot(mcmc.households.mvnorm$theta[,2], type='l',  main=paste0('Trace: ',expression(alpha)))
hist(mcmc.households.mvnorm$theta[,2],  main=expression('Hist: alpha'))
abline(v=mean(mcmc.households.mvnorm$theta[,2]), col=2)

## correlation plot
par(mfrow=c(1,1))
plot(mcmc.households.mvnorm$theta[,1], mcmc.households.mvnorm$theta[,2], col=2)

## plot the 1-by-1 sampling and block sampling traces and correlation plots
par(mfrow=(c(3,2)))
## p
plot(mcmc.run$p, type='l', main=expression('Trace 1-by-1 sampling: p'),ylab='p')
plot(mcmc.households.mvnorm$theta[,1], type='l', main=expression('Trace Block sampling: p'), ylab='p')
## alpha
plot(mcmc.run$alpha, type='l', main=expression('Trace 1-by-1 sampling: alpha'),ylab='alpha')
plot(mcmc.households.mvnorm$theta[,2], type='l',  main=paste0('Trace block: ',expression(alpha)), ylab='alpha')
## corr
plot(mcmc.run$p, mcmc.run$alpha, col=2,pch=4, cex=0.4, xlab='p', ylab='alpha')
plot(mcmc.households.mvnorm$theta[,1], mcmc.households.mvnorm$theta[,2], col=2, pch=4, cex=0.6, xlab='p', ylab='alpha')

dev.copy(pdf,('Output/samplingmethods_comparison.pdf'), width=6, height=6)
dev.off()
###############################################################
###############################################################
### multiple chains for the hell of it 
hist(exp(-mcmc.run$alpha))
summary(exp(-mcmc.run$alpha)) ## p = 0.8113

p.init <- c(0.8, 0.25, 0.5)
alpha.init <- c(0.25, 0.33, 1)
MCMC.runs <- 10000
mcmc.posterior.p <- array(NA, dim=c(MCMC.runs, 3))
mcmc.posterior.alpha <- array(NA, dim=c(MCMC.runs, 3))

for(iiM in 1:3){
  init <- c(p=p.init[iiM], alpha=alpha.init[iiM])
  mcmc.run <- mcmc.households(data=mat, iter=MCMC.runs)
  mcmc.posterior.p[,iiM] <- mcmc.run$p
  mcmc.posterior.alpha[,iiM] <- mcmc.run$alpha
}

par(mfrow=c(2,2))
plot(mcmc.posterior.p[,1], type='l')
lines(mcmc.posterior.p[,2], col=2)
lines(mcmc.posterior.p[,3], col=3)
hist(mcmc.posterior.p)
plot(mcmc.posterior.alpha[,1], type='l')
lines(mcmc.posterior.alpha[,2], col=2)
lines(mcmc.posterior.alpha[,3], col=3)
hist(mcmc.posterior.alpha)