############
## MCMC 1 - LAB2 Gibbs Sampler
setwd("/Users/AHenderson/OneDrive - London School of Hygiene and Tropical Medicine/PhD/SISMID/SISMID_MCMCI")

chainGibbs = function(n1, n11, N3, mcmc.size,alpha,beta){
  # Reserve space
  q    = rep(0,mcmc.size) 
  n111 = rep(0,mcmc.size)
  
  # Initialize the model unknowns
  q[1]    = 0.5;                         # just an initial guess
  n111[1] = round(275*2*q[1]/(2*q[1]+1)) # the (rounded) expected value of n_111, 
  # given q =0.5 
  
  # The observations (cf. the lecture)
  n1  = n1   # frequency of chain 1
  n11 = n11   # frequency of chain 1->1
  N3  = N3  # frequency of chains with outbreak size 3 
  # (the total frequency of chains 1->2 and 1->1->1)

  
  # Draw MCMC samples
  for (i in 2:mcmc.size){
    
    # Draw an iterate of q from its full conditional distribution
    q[i] = rbeta(1, 2*n1 + 2*n11 + n111[i-1] + alpha, n11 + 2*N3 + beta)
    
    # Draw an iterate of n111 from its full conditional distribution
    n111[i] = rbinom(1, N3, (2*q[i])/(2*q[i] + 1))
    
  }
  
  # The output: the MCMC samples
  chainGibbs = list(q=q,n111=n111)
  
}
test <-  chainGibbs(34, 25, 275, 5000, 1, 1)

#pdf(file="/Users/betz/Documents/TexWork/MCMC/betz/Bayesintro/chaingibbs1.pdf", height=4.5, width=8.9)
# Set up to have 2 plots in one figure
par(mfrow=c(1,2),oma=c(0,0,0,0))

# Assume a burn-in of 500 iterations, so just plot those greater than 500. 
burnin=500
mcmc.runs=5000
hist(test$q[burnin+1:mcmc.runs],main="",xlab="q")
hist(test$n111[burnin+1:mcmc.runs],main="",xlab="n111")

# Summary of the results 
summary(test$q[burnin+1:mcmc.runs])
summary(test$n111[burnin+1:mcmc.runs])

dev.copy(pdf, "Output/1_gibbs1.pdf", width=6, height=6)
dev.off()


####################################
## play around with value of N3
par(mfrow=c(3,2))
test2 <-  chainGibbs(34, 25, 27500, 5000, 1, 1)
test3 <-  chainGibbs(34, 25, 50, 5000, 1, 1)
print("N3=275")
summary(test$q[burnin+1:mcmc.runs])
summary(test$n111[burnin+1:mcmc.runs])

print("N3=2750")
summary(test2$q[burnin+1:mcmc.runs])
summary(test2$n111[burnin+1:mcmc.runs])

print("N3=20")
summary(test3$q[burnin+1:mcmc.runs])
summary(test3$n111[burnin+1:mcmc.runs])

hist(test$q[burnin+1:mcmc.runs],xlab="q", main='N3=275')
hist(test$n111[burnin+1:mcmc.runs],main="",xlab="n111")

hist(test2$q[burnin+1:mcmc.runs],xlab="q", main='N3=2750')
hist(test2$n111[burnin+1:mcmc.runs],main="",xlab="n111")

hist(test3$q[burnin+1:mcmc.runs],xlab="q", main='N3=50')
hist(test3$n111[burnin+1:mcmc.runs],main="",xlab="n111")

dev.copy(pdf, "Output/1_gibbs2.pdf", width=6, height=6)
dev.off()

####################################
## Play around with alpha and beta
test2 <-  chainGibbs(34, 25, 275, 5000, 0.1, 0.1)
test3 <-  chainGibbs(34, 25, 275, 5000, 1, 0.1)

hist(test$q[burnin+1:mcmc.runs],xlab="q", main='alpha=1, beta=1')
hist(test$n111[burnin+1:mcmc.runs],main="",xlab="n111")

hist(test2$q[burnin+1:mcmc.runs],xlab="q", main='alpha=0.1, beta=0.1')
hist(test2$n111[burnin+1:mcmc.runs],main="",xlab="n111")

hist(test3$q[burnin+1:mcmc.runs],xlab="q", main='alpha=1, beta=0.1')
hist(test3$n111[burnin+1:mcmc.runs],main="",xlab="n111")

dev.copy(pdf, "Output/1_gibbs3.pdf", width=6, height=6)
dev.off()
