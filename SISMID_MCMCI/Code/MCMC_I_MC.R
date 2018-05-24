## This script illustrates the Importance-Sampling algorithm and Ehrenfest model of diffusion
setwd("/Users/AHenderson/OneDrive - London School of Hygiene and Tropical Medicine/PhD/SISMID/SISMID_MCMCI")

## IMPORTANCE-SAMPLING
## define a threshold value and number of Monte Carlo samples
my_const = 4.5
sim_size = 10000

## true probability of interest
(true_prob = pnorm(my_const,lower.tail=FALSE))


## naive Monte Carlo estimate

## Your task: create a naive and an importance sampling 
## estimate of the normal tail probability. 
## To generate realizations from the shifted exponential 
## use `rexp()` to generate regular exponentials and 
## then add my.const to them. Also, remember that you 
## don't have to code the formula for the normal 
## density, because it is available via `dnorm()'. 
##If you finish early, get Monte Carlo errors for 
## naive and important sampling schemes. 

# naive - random draw from standard normal 
set.seed(05081991)
# Count how many of 10,000 draws are greater than 4.5
naive <- ifelse(rnorm(sim_size, mean=0, sd = 1)>my_const , 1, 0)
(naive.p <- mean(naive))

# importance sampling - shift the dist and draw from exponential
set.seed(05081991)
g.y <- rexp(sim_size, 1)+my_const # g(Y)
# For weighting we need density of draws from std. norm.
theta.y <- dnorm(g.y, 0, 1)
# then re-estimate from the new dist and weight by (density/exp-(y-c))
imp.sampl <- ifelse(g.y>my_const , 1, 0)*theta.y/exp(my_const-g.y)
(imp.sampl.p <- mean(imp.sampl))

par(mfrow=c(2,1))
hist(rnorm(sim_size, mean=0, sd = 1), xlim=c(-4, 14))
hist(g.y,xlim=c(-4, 14))
dev.copy(pdf, "Output/2_impsampl.pdf", width=6, height=6)
dev.off()
# compare means
print(true_prob)
print(naive.p)
print(imp.sampl.p) # so importance sampling is working better

# compare variance of the two estimators
(var.naive <- var(naive)/sim_size)
(var.imp.sampl <- var(imp.sampl)/sim_size) # something sort of logical

## Computing Monte Carlo confidence intervals
## Naive (not very useful, usually variance estimate is degenerate):
c(naive.p-1.96*sqrt(var.naive), naive.p+1.96*sqrt(var.naive)) # 0s so useless

## Importance Sampling:
c(imp.sampl.p-1.96*sqrt(var.imp.sampl), imp.sampl.p+1.96*sqrt(var.imp.sampl))


##' ##
##' EHRENFEST model of diffusion

## this function randomly draws a new state of the Ehrenfest model
## "Transition Kernel"
next_state = function(cur.state, num.mol){
  return.value=NULL
  i <- cur.state
  N <- num.mol
  p <- runif(1,0,1)
  
  if(p < i/N){return.value = cur.state - 1} else{
              return.value = cur.state + 1
              }
  return(return.value)
}

newval <- NULL
newval2 <- NULL
for(i in 1:1000){
newval[i] <- next_state(99,100)
newval2[i] <- next_state(100,100)
}
mean(newval) ## 98.02 - almost always goes down
mean(newval2) ## 99 exactly - so ALWAYS goes down 

## set the number of molecules and the number of iterations
my_num_mol = 100
sim_size = 1000

## initialize the chain by drawing the initial state uniformly at random from all possible states. R function `sample()' will be handy.
my_draws = numeric(sim_size)
set.seed(05081991)
my_draws[1] = sample(0:my_num_mol, 1)

## run the Markov chain

for (i in 2:sim_size){
  ## use next.state function to evolve the Markov chain one step at a time
  my_draws[i] = next_state(my_draws[i-1],my_num_mol)
}

## plot the chain
par(mfrow=c(1,1))
plot(1:sim_size, my_draws, type="l", ylab="Ehhenfest State", xlab="Time Step")

## get the "time averages"

mean(my_draws) ## 47.072
var(my_draws) ## 58.70152

## get the "space averages"
my_num_mol/2 # 50
my_num_mol/4 # 25

# show the impact of increasing sim size
par(mfrow=c(2,1))
plot(1:sim_size, my_draws, type="l", ylab="Ehhenfest State", xlab="Time Steps = 1000")
#now increase sim size
sim_size = 10000
my_draws2 = numeric(sim_size)
set.seed(05081991)
my_draws2[1] = sample(seq(0,my_num_mol,1), 1)
for (i in 2:sim_size){
  ## use next.state function to evolve the Markov chain one step at a time
  my_draws2[i] = next_state(my_draws2[i-1],my_num_mol)
}
mean(my_draws2) ## 49.9424
var(my_draws2)  ## 30.1557
plot(1:sim_size, my_draws2, type="l", ylab="Ehhenfest State", xlab="Time Steps = 1000000")

dev.copy(pdf, "Output/2_mcsims.pdf", width=6, height=6)
dev.off()

