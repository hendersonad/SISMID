## This script illustrates the Metropolis-Hastings algorithm for
## approximating the posterior distribution of the time of infection
## in a simple SIS model
setwd("/Users/AHenderson/OneDrive - London School of Hygiene and Tropical Medicine/PhD/SISMID/SISMID_MCMCI")
set.seed(05081991)

sis_log_like = function(inf_time, inf_rate, clear_rate, total_time){
  return(log(inf_rate) - inf_rate*inf_time - clear_rate*(total_time-inf_time))
}

# finish this function
sis_proposal = function(cur_inf_time, total_time, win_half_len){
  rw.1 <- runif(1, cur_inf_time-win_half_len, cur_inf_time+win_half_len)
  
  # reflective sampling
  if(rw.1 < total_time &  rw.1 >0){
    proposal_inf_time = rw.1
  }else if (rw.1 > total_time){
    proposal_inf_time = (2*total_time) - rw.1
  }else if (rw.1 < 0){
    proposal_inf_time = -rw.1
  }

  return(proposal_inf_time)
}
#Test
inf_rate=0.1; clear_rate=0.2; total_time=1.0; win_half_len=0.2;
sis_proposal(0.1, 0.4, 0.05)


inf_time_mcmc = function(start_inf_time, inf_rate, clear_rate, total_time, win_half_len, chain_len){
  
  result_mat = matrix(0, chain_len, 3)
  colnames(result_mat) = c("inf_time", "log_like", "acc_ind")
  
  cur_inf_time= start_inf_time
  result_mat[1,1] = start_inf_time
  result_mat[1,2] = sis_log_like(start_inf_time, inf_rate, clear_rate, total_time)
  
  for (i in 2:chain_len){
    ## 1. Generate a new value of the infection time using 
    ##  the function sis.proposal()
    ## 2. Decide whether to accept or reject the proposed value by computing
    ## the Metropolis-Hastings ratio
    ## 3. Save the current or proposed value in result.mat[i,1]
    ##    Save the complete-data log-likelihood evaluated either at the current
    ##    or proposed value of the infection time in result.mat[i,2]
    ##    Save the indicator of the acceptance in result.mat[i,3]
    proposal_inf_time <- sis_proposal(result_mat[i-1,1], total_time, win_half_len)
    
    log.mh.alg = sis_log_like(proposal_inf_time, inf_rate, clear_rate, total_time) - 
                  sis_log_like(result_mat[i-1,1], inf_rate, clear_rate, total_time)
                  
    if(log(runif(1)) < min(1,log.mh.alg)){
      result_mat[i,1] <- proposal_inf_time
      result_mat[i,2] <- sis_log_like(proposal_inf_time, inf_rate, clear_rate, total_time)
      result_mat[i,3] <- 1
    }else{
      result_mat[i,1] <- result_mat[i-1,1]
      result_mat[i,2] <- result_mat[i-1,2]
      result_mat[i,3] <- 0
    }
  }
  
  return(result_mat)
}


#set chain length, starting points and tuning parameter
par(mfrow=c(2,2))
inf_rate_range = c(0.1,  2, 0.5, 5)
clear_rate_range = c(2, 0.1, 0.5, 5)
for(i in 1:4){
set.seed(05081991)
chain_len = 10000
start_inf_time = 0.25;
inf_rate = inf_rate_range[i];
clear_rate = clear_rate_range[i];
total_time = 1;

delta <- 0.5

## run the above functions
sample = inf_time_mcmc(start_inf_time=0.25, inf_rate=inf_rate, clear_rate=clear_rate, 
                        total_time=1.0, win_half_len=delta, chain_len=chain_len)

burnin <- 1000
summary(sample[burnin:chain_len,])

plot(burnin:chain_len, sample[burnin:chain_len,"inf_time"], type="l", xlab="Iteration", ylab="Infection Time")

hist(sample[burnin:chain_len,1], xlab='Infection time', main=paste0("Inf_rate=", inf_rate, " clear_rate=", clear_rate))

if(i %% min(i,2) == 0){
dev.copy(pdf, paste0("Output/3_MHalgorithm_continuous", i,".pdf"), width=6, height=6)
dev.off()
}
}

