## to easily calc the density of the gamma dist...
log.gamma.density <- function(theta, lambda, nu) {
  
  out <- (lambda-1)*log(theta) - theta*nu
  return(out)
}
##
##
## MCMC loop - alpha from MH, beta from gibbs
##
mcmc <- function(y, iter, sigma.rw, lambda.alpha, lambda.beta, nu.alpha, nu.beta) {
  
  # find n from the data (ie the number of observations)
  n <- length(y)
  
  # compute the sum of the observations and the sum of the log of the observations
  sum.y <- sum(y)
  sum.log.y <- sum(log(y))
  
  # first create a table which will store the samples; first column will be alpha, second will be beta
  out <- matrix(NA, nrow=iter, ncol=2)
  
  # initial values
  alpha.cur <- 1
  beta.cur <- 1
  out[1,] <- c(alpha.cur, beta.cur)
  
  # mcmc loop starts here
  for (i in 2:iter) {
    
    ###############
    # update alpha (assume beta is fixed)
    ###############
    
    # propose a new value for alpha
    alpha.can <- rnorm(1, alpha.cur, sigma.rw)
    
    # if it is negative reject straight away else compute the M-H ratio
    if (alpha.can > 0) {
      
      # evaluate the loglikelihood at the current values of alpha.
      loglik.cur <- n*alpha.cur*log(beta.cur) - n*lgamma(alpha.cur) + (alpha.cur-1)*sum.log.y 
      
      # compute the log-likelihood at the candidate value of alpha
      loglik.can <- n*alpha.can*log(beta.cur) - n*lgamma(alpha.can) + (alpha.can-1)*sum.log.y 
      
      
      # log prior densities
      log.prior.alpha.cur <- log.gamma.density(alpha.cur, lambda.alpha, nu.alpha)
      log.prior.alpha.can <- log.gamma.density(alpha.can, lambda.alpha, nu.alpha)
      
      logpi.cur <- loglik.cur + log.prior.alpha.cur
      logpi.can <- loglik.can + log.prior.alpha.can
      
      # M-H ratio
      
      # draw from a U(0,1)
      u <- runif(1)
      
      if (log(u) < logpi.can - logpi.cur) {
        alpha.cur <- alpha.can
      }
    }
    
    ###############
    # update beta
    ###############
    
    beta.cur <- rgamma(1, alpha.cur*n + lambda.beta, rate = sum.y + nu.beta)  
    
    # store the samples
    out[i,] <- c(alpha.cur, beta.cur)
    
  }
  
  return(out)
}

##'
##'Generate some data
set.seed(05081991)
data <- rgamma(500, 4, 2)
## run mcmc
res <- mcmc(y = data, iter = 10000, sigma.rw = 1, lambda.alpha = 1, lambda.beta = 1, nu.alpha = 1e-4, nu.beta = 1e-4)
## plot traces
par(mfrow=c(1,2))
plot(res[,1],type='l', main=expression(paste(alpha)))
plot(res[,2],type='l', main=expression(paste(beta)))

## burnin 1000 and plot hist
burnin <- 1000
resburn <- res[-c(1:burnin),]
plot(resburn[,1],type='l', main=expression(paste(alpha)))
plot(resburn[,2],type='l', main=expression(paste(beta)))
hist(resburn[,1], prob=TRUE, col=4, xlab=expression(paste(alpha)), main=expression(paste("Posterior distribution of ", alpha)))
hist(resburn[,2], prob=TRUE, col =4, xlab=expression(paste(beta)), main=expression(paste("Posterior distribution of ", beta)))

## joint dist of parameters
par(mfrow=c(1,1))
plot(res, xlab=paste(expression(alpha)), ylab=paste(expression(beta)))
