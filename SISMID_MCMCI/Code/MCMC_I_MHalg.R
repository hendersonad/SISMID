## This script illustrates the Metropolis-Hastings algorithm for
## approximating the standard normal distribution
setwd("/Users/AHenderson/OneDrive - London School of Hygiene and Tropical Medicine/PhD/SISMID/SISMID_MCMCI")
set.seed(05081991)
## Your task: Add code the following function. The input of  
## the function is the current value of the random walk and 
## the tuning parameter (\delta in the lab notes). The output
## should be a vector with the first component being the next 
## value of the random walk and the second component being 0 if 
## the proposed value was rejected and 1 if accepted.
cur_value=3
tuning_par=5
unif_rw_next = function(cur_value, tuning_par){
  rw.1 <- runif(1, -tuning_par, tuning_par)
  
  proposal_value = cur_value + rw.1
  #mh.value = PIy.q(x|y) / PIx.q(y|x) but symmetric so the Qs cancel and only care about PI
  # PI is the density at point x or point y (ignoring normalising constant)
  mh.value = exp((-proposal_value^2)/2)/
              exp((-cur_value^2)/2)
  log.mh = (cur_value^2 - proposal_value^2)/2  # usually use log scale for numerical stability 
  
  return_value = c(cur_value, 0)
  if(log(runif(1)) < min(log.mh, 1)){
    return_value[1] <- proposal_value
    return_value[2] <- 1
  } else {
    return_value[1] <- cur_value
    return_value[2] <- 0
  }
  return(return_value)      
}
# test it:
unif_rw_next(10, 2)

mcmc_size = 10000
start.value = 3.0


mcmc_out = matrix(0, nrow=(mcmc_size), ncol=2)
colnames(mcmc_out) = c("state", "acc.status")

cur_value = c(start.value,1)

## Your task: add code to the below for loop to 
## fill in the mcmc.out matrix defined above with 
## the first column recording the state of the random 
## walk and the second column recording the acceptance status
## of each Metropolis-Hastings move. Don't forget to play 
## with the tuning parameter \delta.

for (i in 2:mcmc_size){
  delta <- 5
  
  mcmc_out[i,] <- unif_rw_next(mcmc_out[i-1,1], delta)
  
}

## acceptance probability
mean(mcmc_out[,"acc.status"])

## mean of the target distribution
mean(mcmc_out[,"state"])

## trace plot
plot(1:mcmc_size, mcmc_out[,"state"], type="l", xlab="Iteration", ylab="MCMC State")

## histogram
hist(mcmc_out[,"state"], xlab="MCMC State", main="Target Distribution Histogram")

dev.copy(pdf, "Output/3_MHalgorithm.pdf", width=6, height=6)
dev.off()
