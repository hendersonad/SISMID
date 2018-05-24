setwd('/Users/AHenderson/OneDrive - London School of Hygiene and Tropical Medicine/PhD/SISMID/SISMID_MCMC2')

source("Code/2coding.R")
source("Code/3mcmc-nonmarkov.R")
data <- read.table("Data/data.txt", header=TRUE)

set.seed(05081991)
mcmc.run <- mcmcSIR.NonMarkov(data, 1e4)
head(mcmc.run)
colnames(mcmc.run) <- c("Beta", "Gamma", "Sum inf times", "Acceptance rate")

par(mfrow=(c(4,1)))
plot(mcmc.run[,1], type='l')
plot(mcmc.run[,2], type='l')
plot(mcmc.run[,3], type='l')
hist(mcmc.run[,4])


summary(mcmc.run[,1])
summary(mcmc.run[,2])
summary(mcmc.run[,4])



