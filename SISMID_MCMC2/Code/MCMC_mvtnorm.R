mcmc.households.mvnorm <- function(data, iter){
  #init
  p.cur <- init['p']
  alpha.cur <- init['alpha']
  
  # init priors
  p.prior.cur <- dbeta(p.cur, shape1=p.prior[1],shape2=p.prior[2])
  alpha.prior.cur <- dgamma(alpha.cur, shape=alpha.prior[1], rate=alpha.prior[2])
  
  #i=1 init values, prior ddensity and likelihood
  p <- rep(NA, iter)
  p[1] <- p.cur
  alpha <- rep(NA, iter)
  alpha[1] <- alpha.cur
  accepttab <- rep(NA, iter)
  accepttab[1] <- 0
  
  theta <- cbind(p, alpha)
  covmat <- diag(c(0.1,1))
  
  theta.cur.prior <- dmvnorm(theta[1,], mean=theta[1,], sigma=covmat)
  lik.cur <- loglik(alpha=theta[1,2], p=theta[1,1])
  
  
  for(i in 2:iter){
    #proposal theta
    theta.prop <- rmvnorm(1, mean=theta[i-1,], sigma=covmat)
    
    if(theta.prop[1]>0 & theta.prop[1]<=1 & theta.prop[2] > 0){

      #likelihoods
      lik.prop <- loglik(alpha=theta.prop[2], p=theta.prop[1])
      
      #priors
      theta.prop.prior <- dmvnorm(theta.prop, mean=theta[i-1,], sigma=covmat)
      
      #MH number
      MH.alg <- (theta.prop.prior/theta.cur.prior)*exp(lik.prop - lik.cur)
      
      #selection decision
      if(runif(1) < min(1, MH.alg)){
        lik.cur <- lik.prop
        theta.cur.prior <- theta.prop.prior
        theta[i,] <- theta.prop 
        accepttab[i] <- 1
      } else{
        theta[i,] <- theta[i-1,] 
        accepttab[i] <- 0
      } 
    }else{
      theta[i,] <- theta[i-1,]  
      accepttab[i] <- 0
    }
  }
  return(list(theta=theta, accept=accepttab))
  
}
