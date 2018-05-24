#SISMID_lik_1
# install appropriate software:
update.packages()
#update.packages() - updates all packages...
source("https://kingaa.github.io/sbied/prep/packages.R")
source("https://kingaa.github.io/sbied/prep/pompTest.R")

## Ricker model
library(ggplot2)
library(reshape2)
library(pomp)
stopifnot(packageVersion("pomp")>="1.12")
pompExample(ricker)
plot(ricker)
x <- simulate(ricker)
class(x) # created a POMP
plot(x) # now we see 3 series
#we can turn the pomp object into a Dframe
y <- as.data.frame(ricker)
head(y) #y
head(simulate(ricker,as.data.frame=TRUE)) #y N e 
#Now run several sims at once
x <- simulate(ricker,nsim=10)
class(x) # now a list of 10 pomp objects
sapply(x,class)
x <- simulate(ricker,nsim=10,as.data.frame=TRUE)
head(x)
tail(x)
str(x) #'sim' is a var 1-10 with simulation number
#Plot several sims against the data
x <- simulate(ricker,nsim=9,as.data.frame=TRUE,include.data=TRUE)
ggplot(data=x,aes(x=time,y=y,group=sim,color=(sim=="data")))+
  geom_line()+guides(color=FALSE)+
  facet_wrap(~sim,ncol=2)
#the deterministic map is the 'skeleton' of the stoch map
#we can compute a trajectory of skeleton using trajectory
y <- trajectory(ricker)
dim(y)
dimnames(y)
plot(time(ricker),y["N",1,],type="l")
#We can extract or set the parameters in the pomp object using `coef`:
coef(ricker) # r=44.70118, N=7
coef(ricker,"phi") 
coef(ricker) <- c(phi=20,c=1,r=44,sigma=0.3,N.0=10,e.0=0)
coef(ricker)
coef(ricker,c("phi","c")) <- c(10,2)
coef(ricker) # r=44, N=10
#' There are a number of other examples included with the package.
#' Do `pompExample()` to see a list of these.

## Inference algorithms in pomp
#lots of inference algorithms available.
#for now use a simple particle filter. 
#' It can be used to evaluate the likelihood at a particular set of parameters.
#' One uses the `Np` argument to specify the number of particles to use:
pompExample(ricker)
pf <- pfilter(ricker,Np=1000)
class(pf) #pomp
plot(pf) #cond.lik, ess, y
logLik(pf) #-138.3507
#run it again
pf <- pfilter(pf)
logLik(pf) #-139.4748
#its a Monte Carlo alg. so get a slightly different answer
#' Note that, by default, running `pfilter` on a `pfilterd.pomp` object causes the computation to be re-run with the same parameters as before.
#' Any additional arguments we add override these defaults.
pf <- pfilter(pf,Np=100)
logLik(pf) #-139.8989


## Building a custom pomp object
#' The usefulness of **pomp** in scientific research hinges on its facilities for implementing the full range of POMP models.
#' To get started building custom `pomp` models, see this [introductory tutorial](./ricker.html).

#######################
#######################
## L2 - Simulation of system of stochastic dynamic models
set.seed(594709947L)
#boarding school example
read.table("https://kingaa.github.io/sbied/stochsim/bsflu_data.txt") -> bsflu
head(bsflu)
#model movement from S to I & I to R as rbinom
sir_step <- Csnippet("
  double dN_SI = rbinom(S,1-exp(-Beta*I/N*dt));
                     double dN_IR = rbinom(I,1-exp(-gamma*dt));
                     S -= dN_SI;
                     I += dN_SI - dN_IR;
                     R += dN_IR;
                     ")
#init conditions (1 index case)
sir_init <- Csnippet("
                     S = N-1;
                     I = 1;
                     R = 0;
                     ")
#' We fold these `Csnippets`, with the data, into a `pomp` object thus:
pomp(bsflu,time="day",t0=0,rprocess=euler.sim(sir_step,delta.t=1/6),
     initializer=sir_init,paramnames=c("N","Beta","gamma"),
     statenames=c("S","I","R")) -> sir
#adding a variable H - to track the incidence
sir_step <- Csnippet("
  double dN_SI = rbinom(S,1-exp(-Beta*I/N*dt));
  double dN_IR = rbinom(I,1-exp(-gamma*dt));
  S -= dN_SI;
  I += dN_SI - dN_IR;
  R += dN_IR;
  H += dN_IR;
")
sir_init <- Csnippet("
                     S = N-1;
                     I = 1;
                     R = 0;
                     H = 0;
                     ")
pomp(sir,rprocess=euler.sim(sir_step,delta.t=1/6),initializer=sir_init,
     paramnames=c("Beta","gamma","N"),statenames=c("S","I","R","H")) -> sir
#case reports (B) as a binom process rbinom(Ht - Ht0-1, p)
#then reset B to 0 after each observation
pomp(sir,zeronames="H") -> sir
#Now, to include the observations in the model, we must write both a dmeasure and an rmeasure component:
dmeas <- Csnippet("lik = dbinom(B,H,rho,give_log);")
rmeas <- Csnippet("B = rbinom(H,rho);")
#' and put these into our `pomp` object:
sir <- pomp(sir,rmeasure=rmeas,dmeasure=dmeas,statenames="H",paramnames="rho")

## Testing the model - simulations
#some sensible params (R0=1.5, 1 day infectious, so beta=1.5, and 2600 pop, rho=rep rate)
sims <- simulate(sir,params=c(Beta=1.5,gamma=1,rho=0.9,N=2600),
                 nsim=20,as.data.frame=TRUE,include.data=TRUE)
ggplot(sims,mapping=aes(x=time,y=B,group=sim,color=sim=="data"))+
  geom_line()+guides(color=FALSE)
#fiddle to make more realistic
sims <- simulate(sir,params=c(Beta=2,gamma=0.4,rho=0.6,N=2600),
                 nsim=20,as.data.frame=TRUE,include.data=TRUE)
ggplot(sims,mapping=aes(x=time,y=B,group=sim,color=sim=="data"))+
  geom_line()+guides(color=FALSE)

## make an seir model in pomp
seir_step <- Csnippet("
  double dN_SE = rbinom(S,1-exp(-Beta*I/N*dt));
                     double dN_EI = rbinom(E,1-exp(-alpha*dt));
                     double dN_IR = rbinom(I,1-exp(-gamma*dt));
                     S -= dN_SE;
                     E += dN_SE - dN_EI;
                     I += dN_EI - dN_IR;
                     R += dN_IR;
                     H += dN_IR;
                     ")
seir_init <- Csnippet("
                     S = N-1;
                     E = 0;
                     I = 1;
                     R = 0;
                     H = 0;
                     ")
pomp(bsflu,time="day",t0=0,rprocess=euler.sim(seir_step,delta.t=1/6),
     initializer=seir_init,paramnames=c("N","Beta","alpha","gamma"),
     statenames=c("S","E","I","R","H")) -> seir
pomp(seir,zeronames="H") -> seir
dmeas <- Csnippet("lik = dbinom(B,H,rho,give_log);")
rmeas <- Csnippet("B = rbinom(H,rho);")
#' and put these into our `pomp` object:
seir <- pomp(seir,rmeasure=rmeas,dmeasure=dmeas,statenames="H",paramnames="rho")
set.seed(594709947L)
sims <- simulate(seir,params=c(Beta=6,gamma=0.4,alpha=1,rho=0.7,N=2600),
                 nsim=20,as.data.frame=TRUE,include.data=TRUE)
ggplot(sims,mapping=aes(x=time,y=B,group=sim,color=sim=="data"))+
  geom_line()+guides(color=FALSE)


  