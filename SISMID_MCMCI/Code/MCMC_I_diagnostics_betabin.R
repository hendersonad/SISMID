## This script illustrates some diagnostic tools available in package coda
## Author: Vladimir N. Minin
## last update: 07/16/17
setwd("/Users/AHenderson/OneDrive - London School of Hygiene and Tropical Medicine/PhD/SISMID/SISMID_MCMCI")

## first we need to load coda and mcmcse packages (you need to install them first)
library(coda)
library(mcmcse)

## now let's load our M-H and Gibbs sampler examples
source("https://raw.githubusercontent.com/vnminin/SISMID_MCMC_I/master/2017/code/chainGibbs.R")
source("/Users/AHenderson/OneDrive - London School of Hygiene and Tropical Medicine/PhD/SISMID/SISMID_MCMCI/Code/MCMC_I_BetaBin_functionsonly.R")

dev.off()

## run one M-H example chain
rat_data = read.table("https://raw.githubusercontent.com/vnminin/SISMID_MCMC_I/master/2016/code/rat_tumor.txt", header=TRUE)


rat_results = mcmc_sampler(my_data=rat_data,
                           init_alpha=1.0,
                           init_beta=1.0,
                           prior_inten_alpha=0.1,
                           prior_inten_beta=0.1,
                           tuning_alpha=0.7,
                           tuning_beta=0.7,
                           mcmc_size=110000,
                           mcmc_burnin=10000,
                           mcmc_subsample=10)
head(rat_results)


## convert the output into coda format - save data as MCMC objects
## given matrix (data created from gibbs.chain) it will make it coda-readable
coda.gibbs.chain = mcmc(cbind(alpha=rat_results[,'alpha'],beta=rat_results[,'beta']))
class(coda.gibbs.chain) ## 'mcmc'

head(coda.gibbs.chain) #col1 = q(escape prob), col2 = Number of households sampled
summary(coda.gibbs.chain)

## we can also compute Monte Carlo error for the quantiles
## we care most about the quantiles - what are the extremes that are plausible from data?
## do we just trust the quantile reported by summary()??
## Probably shouldn't so what is the MCerror on these quantiles...
## i.e. if you just rerun the code then could the quantile jump to 0.1? or is it accurate...
mcse.q.mat(coda.gibbs.chain, 0.025) 
mcse.q.mat(coda.gibbs.chain, 0.975) # MCerror tiny for both

## plot traceplots and (hopefully) posterior densities
plot(coda.gibbs.chain) # coda knows the class is MCMC so does funky useful plots - trace and density

## look at what coda can do
#help(package=coda)

## look at the menu options
#codamenu()
##' awesome little interactive menu with all the options from the coda
##' package...

## all commands are availabe outside of the menu

## plot autocorrelations plots
autocorr.plot(coda.gibbs.chain)
##' Markov chain so only want dependency on t-1... definitely not on lag2,3,4,...
##' at lag=0 AC should always be 1
##' at lag=1 AC should have some dependency because it's a chain
##' beyond that should be very small then chain is mixing well
##' 
##' If there is strong AC for many lAgs then the less IID your samples are
##' i.e. MCMC IF samples were IID then would take many fewer samples (effective sample size)

## calculate effective sample size
effectiveSize(coda.gibbs.chain)
##' we ran chain for 5000 iterations
##' IF we were able to run IID then our results correspond to 3000 iterations 


################################################################
##' CONVERGENCE
##' There is no formula to answer question: Have my chains converged? 
##' So we have other checks

## Run 50 chains with overdispersed starting values
coda.gibbs.chains = list()

for (i in 1:50){
                ## new function to sample init values for n11 between 0 and 275
  rat_results = mcmc_sampler(my_data=rat_data,
                             init_alpha=runif(1),
                             init_beta=runif(1),
                             prior_inten_alpha=0.1,
                             prior_inten_beta=0.1,
                             tuning_alpha=0.7,
                             tuning_beta=0.7,
                             mcmc_size=110000,
                             mcmc_burnin=10000,
                             mcmc_subsample=10)
  coda.gibbs.chains[[i]] = mcmc(cbind(alpha=rat_results[,'alpha'],beta=rat_results[,'beta']))

}

coda.gibbs.list = mcmc.list(coda.gibbs.chains)

plot(coda.gibbs.list)

## compute and plot Gelman-Rubin potential reduction factor
gelman.diag(coda.gibbs.list)
## potential scale reduction factor - how far could chains be from each other
## if = 1 then converged
## think this is the between chain and within chain ratio

gelman.plot(coda.gibbs.list)
## Shows G-R value if had stopped at iteration = n


dev.copy(pdf, "Output/1_mcmc_multichain.pdf", width=6, height=6)
dev.off()

##' note: convergence isn't mixing. 
##' You could have converged but be mixing very slowly 
##' Usually they work in tandem 