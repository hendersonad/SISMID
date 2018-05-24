
## first we need to load coda and mcmcse packages (you need to install them first)
library(coda)
library(mcmcse)

coda.gibbs.1by1 = mcmc(cbind(mcmc.run$alpha[],mcmc.run$p[]))
coda.gibbs.block = mcmc(cbind(mcmc.households.mvnorm$theta[,2], mcmc.households.mvnorm$theta[,1]))

head(coda.gibbs.1by1) #col1 = q(escape prob), col2 = Number of households sampled
head(coda.gibbs.block) #col1 = q(escape prob), col2 = Number of households sampled
summary(coda.gibbs.1by1)
summary(coda.gibbs.block)

## we can also compute Monte Carlo error for the quantiles
## we care most about the quantiles - what are the extremes that are plausible from data?
## do we just trust the quantile reported by summary()??
## Probably shouldn't so what is the MCerror on these quantiles...
## i.e. if you just rerun the code then could the quantile jump to 0.1? or is it accurate...
mcse.q.mat(coda.gibbs.1by1, 0.025) 
mcse.q.mat(coda.gibbs.1by1, 0.975) 
mcse.q.mat(coda.gibbs.block, 0.025) 
mcse.q.mat(coda.gibbs.block, 0.975) 

## plot traceplots and (hopefully) posterior densities
plot(coda.gibbs.1by1) # coda knows the class is MCMC so does funky useful plots - trace and density
plot(coda.gibbs.block) # coda knows the class is MCMC so does funky useful plots - trace and density

## plot autocorrelations plots
autocorr.plot(coda.gibbs.1by1)
autocorr.plot(coda.gibbs.block)

effectiveSize(coda.gibbs.1by1) ##817 and 511
effectiveSize(coda.gibbs.block) ## 95 and 125!!!!