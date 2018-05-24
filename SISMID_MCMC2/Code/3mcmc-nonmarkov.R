
mcmcSIR.NonMarkov <- function(data, iter) {
  
  #Data frame "data" should contain two columns; first is the label of the individuals and the second contains the removal times.
  # Get the final size and the size of population
  fs <- sum(data$removal != Inf)
  N <- nrow(data)
  
  # Get a vector for the removal times
  removal.times <- data$removal[1:fs]
  
  # First of all, we need to simulate the initial values for the infection
  # times to iniate our MCMC sampler.
  
  # an easy way to place all the infection times just before the first removal time.
  infection.times <- seq(0, min(removal.times), len=fs)
  
  # We "augment" the data frame "data" with the infection times.
  data$infection <- c(infection.times, rep(Inf, N-fs))
  data <- data[,c(1,3,2)] # re order
  
  # Assume  for the parameters beta and gamma
  lambda.beta <- 1.0
  nu.beta <- 10^(-3)
  
  lambda.gamma <- 1.0  
  nu.gamma <- 10^(-3)
  
  # compute the current values of the double sum and the log(product).
  double.sum.cur <- compute.total.pressure(data, fs);  ## double sum i.e. integral 
  log.prod.cur <- compute.log.prod(data,infection.times, removal.times, fs); ## log(product(Itimes))
  sum.R.minus.I.cur <- sum(removal.times-infection.times) ## the other integral bit
  sum.log.R.minus.I.cur <- sum(log(removal.times-infection.times))
  
  # initial values for beta and gamma - draw from the full conditionals.
  # gen alpha
  alpha.cur <- 2
  beta.cur <- rgamma(1, fs - 1 + lambda.beta, 1.0)/(nu.beta + double.sum.cur/N)
  gamma.cur <- rgamma(1, alpha.cur * fs + lambda.gamma, 1.0)/(nu.gamma + sum.R.minus.I.cur);
  
  #gammadensity.cur <- sum(log(dgamma(removal.times-infection.times, alpha.cur, gamma.cur)))
  #gammadensity.cur <- (alpha.cur-1)*sum(log(removal.times-infection.times)) - gamma.cur*sum.R.minus.I.cur
  #print(gammadensity.cur)
  # create a matrix to store the values for beta, gamma and I.
  res <- matrix(NA, nrow=iter, ncol=4)
  res[1,] <- c(beta.cur, gamma.cur, sum(infection.times), accept=0)
  
  ##########################
  # MCMC loop starts here  # 
  ##########################
  
  for (i in 2:iter) {
    ########################
    # update infection times
    ########################
    
    # first choose an individual whose infection time will be updated.
    choose.ind <- sample(1:fs, 1)
    
    # record the current value and propose a candidate value
    inf.time.cur <- infection.times[choose.ind]
    (inf.time.can <- removal.times[choose.ind] - rexp(1, gamma.cur)); 
    
    # compute the (new) value of the double sum, log(prod) and sum.R.minus.I
    
    # first update the data frame data which is used for the computation of the quantities of interest
    data$infection[choose.ind] <- inf.time.can;
    infection.times[choose.ind] <- inf.time.can;
    
    # compute the log product first. The proposed move should be consistent with the observed data,
    # i.e. results into the same final size. Otherwise, the move is rejected straight away.
    
    log.prod.can <- compute.log.prod(data, infection.times, removal.times, fs);
    #print(log.prod.can)
    if (log.prod.can != -Inf) { 
      
      # compute the other two quantities
      double.sum.can <- compute.total.pressure(data, fs);
      sum.R.minus.I.can <- sum(removal.times-infection.times);
      sum.log.R.minus.I.can <- sum(log(removal.times-infection.times))
      
      # compute the log(q_ratio);
      # prob of choosing i-time is not symmetric so compute ratio
      #log.q.ratio <- log(dexp(removal.times[choose.ind] - inf.time.cur, gamma.cur)) - log(dexp(removal.times[choose.ind] - inf.time.can, gamma.cur))
      log.q.ratio <- log(dexp(removal.times[choose.ind] - inf.time.cur, gamma.cur)) - log(dexp(removal.times[choose.ind] - inf.time.can, gamma.cur))
      
      #gammadensity.can <- sum(log(dgamma(removal.times-infection.times, alpha.cur, gamma.cur)))
      #gammadensity.can <- (alpha.cur-1)*sum(log(removal.times-infection.times)) - gamma.cur*sum.R.minus.I.can
      
      #log.pi.can <- log.prod.can - (beta.cur/N)*double.sum.can + gammadensity.can;
      #log.pi.cur <- log.prod.cur - (beta.cur/N)*double.sum.cur + gammadensity.cur;
      
      log.pi.can <- log.prod.can - (beta.cur/N)*double.sum.can + (alpha.cur-1)*sum.log.R.minus.I.can - gamma.cur*sum.R.minus.I.can
      log.pi.cur <- log.prod.cur - (beta.cur/N)*double.sum.cur + (alpha.cur-1)*sum.log.R.minus.I.cur - gamma.cur*sum.R.minus.I.cur
      
      #print(c(log.pi.can - log.pi.cur + log.q.ratio, log.pi.can, log.pi.cur, log.q.ratio))
      # accept/reject move
      u <- runif(1);
      
      if (log(u) < log.pi.can - log.pi.cur + log.q.ratio) {
        
        log.prod.cur <- log.prod.can;
        sum.log.R.minus.I.cur <- sum.log.R.minus.I.can;
        double.sum.cur <- double.sum.can;
        sum.R.minus.I.cur <- sum.R.minus.I.can;
        accepttab <- 1;
      }
      else {
        data$infection[choose.ind] <- inf.time.cur;
        infection.times[choose.ind] <- inf.time.cur;
        accepttab <- 0;
      }      
    }
    else {
      data$infection[choose.ind] <- inf.time.cur;
      infection.times[choose.ind] <- inf.time.cur;
      accepttab <- NA;
    }
    
    #################    
    # update beta
    #################
    beta.cur <- rgamma(1, fs - 1 + lambda.beta, 1.0)/(nu.beta + double.sum.cur/N)
    
    #################    
    # update gamma
    #################    
    gamma.cur <- rgamma(1, (alpha.cur*fs) + lambda.gamma, 1.0)/(nu.gamma + sum.R.minus.I.cur)
    
    
    # store the current values of beta, gamma and the sum of the infection times.
    res[i,] <- c(beta.cur, gamma.cur, sum(infection.times), accept = accepttab)
    
  }
  
  res
  
}
