
# this is a function to compute the loglikelihood for any pair of alpha, p for the given data in the hand out;
loglik <- function(alpha, p) {
  
  29*log(compute.marginal.prob(0,3,1,alpha, p)) +  9*log(compute.marginal.prob(1,3,1,alpha, p)) +  2*log(compute.marginal.prob(2,3,1,alpha, p)) +  2*log(compute.marginal.prob(3,3,1,alpha, p)) 
  
}

# this is an MCMC algorithm for Asian-Data household model
# iter: number of MCMC iteration
# sigma: the standard deviation for the random walk Metropolis algorithm

mcmc.hh <- function(iter, sigma) {
  
  
  # values for the prior hyper-parameters
  lambda.p <- 1.0
  nu.p <- 1.0
  
  lambda.alpha <- 1.0  
  nu.alpha <- 10^(-3)
  
  
  # initial values for the Markov Chain
  alpha.cur <- 0.5;
  p.cur <- 0.8;
  
  # create a matrix to store the output.
  res <- matrix(NA, nrow = iter, ncol=2);
  res[1,] <- c(alpha.cur, p.cur);
  
  
  for (i in 2:iter) {
    
    
    ##############
    # update alpha
    ##############
    
    # propose new value
    alpha.can <- rnorm(1, alpha.cur, sigma)
    
    # alpha is strictly positive, therefore any proposed value which is negative is automatically rejected.
    if (alpha.can > 0.0) {
      
      # otherwise, if it is positive the accept it with some probability:
      
      # compute the log-densities
      log.pi.cur <- ############# what should we put in here?
        log.pi.can <- ############# what should we put in here?
        
        # log-qratio is zero since we are doing a random Walk Metropolis
        log.q.ratio <- 0;
      
      u <- runif(1)
      if (log(u) < log.pi.can - log.pi.cur) {
        alpha.cur <- alpha.can
      }
    }
    
    
    ############
    # propose p
    ############
    p.can <- runif(1, 0, 1)
    
    
    # compute the log-densities
    log.pi.cur <- ################# what should we put in here?
      log.pi.can <- ################# what should we put in here?
      
      # log-qratio is zero
      log.q.ratio <- 0
    u <- runif(1)
    if (log(u) < log.pi.can - log.pi.cur) {
      p.cur <- p.can
    }
    
    
    # store the output
    res[i,] <- c(alpha.cur, p.cur)
    
  }
  
  res                                             
  
}