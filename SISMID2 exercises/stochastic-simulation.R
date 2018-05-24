### R code from vignette source 'stochastic-simulation.rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: stochastic-simulation.rnw:70-83
###################################################
SI.onestep <- function (x, params) {     #function for one step of the stochastic SI epidemic
   X <- x[2]                             #the second element of x is number of susceptibles X
   Y <- x[3]                             #the third element of x is number of infecteds Y
   with(
     as.list(params),
     {
       new.Y <- Y+1                      #whenever an event occurs we increase infecteds by 1...
       new.X <- X-1                      #and descrease susceptibles by 1
       tau <- rexp(n=1,rate=params$beta*X*Y/(X+Y))  #exponential random time to next event
       c(tau=tau,X=new.X,Y=new.Y)        #store result
       }
     )
}


###################################################
### code chunk number 2: stochastic-simulation.rnw:88-97
###################################################
SI.model <- function (x, params, nstep) { #function to iterate the stochastic SI for nstep events
  output <- array(dim=c(nstep+1,3))       #set up an array to store all the results
  colnames(output) <- c("time","X","Y")   #name the variables in the array
  output[1,] <- x                         #the first record of the array is the initial condition
  for (k in 1:nstep) {                    #iterate the model for nstep events
    output[k+1,] <- x <- SI.onestep(x,params) #update x and store result
  }
  output                                  #return output
}


###################################################
### code chunk number 3: stochastic-simulation.rnw:104-121
###################################################
set.seed(38499583)                        #set the random seed so results are repeatable
nsims <- 10                               #number of simulations to run
pop.size <- 200                           #total size of the population
Y0 <- 2                                   #initial number infected
nstep <- pop.size-Y0                      #how many steps to run? until everyone infected
xstart <- c(time=0,X=(pop.size-Y0),Y=Y0)  #initial conditions
params <- list(beta=0.3) #parameters
data <- vector(mode='list',length=nsims)  #create a list called ``data'' to store all runs
for (k in 1:nsims) {                      #simulate k different runs
  data[[k]] <- as.data.frame(SI.model(xstart,params,nstep))  #main simulation step
  data[[k]]$cum.time <- cumsum(data[[k]]$time) #calculates the running sum of inter-event intervals
}
max.y<-max(data[[1]]$cum.time)               #find the maximum time any simulation ran (to set x axis)
plot(c(0,pop.size),c(0,pop.size),type='n',xlab='time',ylab='incidence',xlim=c(0,max.y)) #set up plot
for (k in 1:nsims) {                      #loop over each simulation...
  lines(Y~cum.time,data=data[[k]],col=k,type='o')  #to plot
}


###################################################
### code chunk number 4: stochastic-simulation.rnw:137-162
###################################################

SIR.onestep <- function (x, params) {  #function to calculate one step of stochastic SIR
  X <- x[2]                            #local variable for susceptibles
  Y <- x[3]                            #local variable for infecteds
  Z <- x[4]                            #local variable for recovereds
  N <- X+Y+Z                           #total population size (subject to demographic change)
  with(                                #use with as in deterministic model to simplify code
       as.list(params), 
       {
         total.rate <- mu*N+beta*X*Y/N+mu*X+mu*Y+gamma*Y+mu*Z #calculate ``total rate''
         tau <- rexp(n=1,rate=total.rate)                     #inter-event time
         new.xyz <- c(X,Y,Z) #initialize a local variable at previous state variable values
         U <- runif(1)       #uniform random deviate
         new.xyz<-c(X,Y,Z-1) #death of recovered id ``default''      
         if (U<=(mu*N+beta*X*Y/N+mu*X+gamma*Y+mu*Y)/total.rate) new.xyz<-c(X,Y-1,Z) 
			#death of infected
         if (U<=(mu*N+beta*X*Y/N+mu*X+gamma*Y)/total.rate) new.xyz<-c(X,Y-1,Z+1)    
			#recovery of infected
         if (U<=(mu*N+beta*X*Y/N+mu*X)/total.rate) new.xyz<-c(X-1,Y,Z)    #death of a susceptible
         if (U<=(mu*N+beta*X*Y/N)/total.rate) new.xyz<-c(X-1,Y+1,Z)       #transmission event
         if (U<=(mu*N/total.rate)) new.xyz<-c(X+1, Y, Z)                  #birth of susceptible
         c(tau,new.xyz) #store result
       }
       )
}


###################################################
### code chunk number 5: stochastic-simulation.rnw:167-177
###################################################

SIR.model <- function (x, params, nstep) {  #function to simulate stochastic SIR
  output <- array(dim=c(nstep+1,4))         #set up array to store results
  colnames(output) <- c("time","X","Y","Z") #name variables
  output[1,] <- x                           #first record of output is initial condition
  for (k in 1:nstep) {                      #iterate for nstep steps
    output[k+1,] <- x <- SIR.onestep(x,params)
  }
  output                                    #return output
}


###################################################
### code chunk number 6: stochastic-simulation.rnw:182-201
###################################################
set.seed(38499583)                #set seed
nsims <- 10                       #number of simulations
pop.size <- 100                   #total population size
Y0 <- 8                           #initial number infected
X0 <- round(0.98*pop.size)        #initial number suscepitlble (~98% of population)
nstep <- 1600                     #number of events to simulate
xstart <- c(time=0,X=X0,Y=Y0,Z=pop.size-X0-Y0) #initial conditions
params <- list(mu=0.00001,beta=60,gamma=365/13) #parameters
data <- vector(mode='list',length=nsims) #initialize list to store the output
for (k in 1:nsims) {              #simulate nsims times
  data[[k]] <- as.data.frame(SIR.model(xstart,params,nstep))
  data[[k]]$cum.time <- cumsum(data[[k]]$time)
}
max.time<-data[[1]]$cum.time[max(which(data[[1]]$Y>0))] #maximum time in first simulation
max.y<-1.8*max(data[[1]]$Y)       #find max infected in run 1 and increase by 80% for plot
plot(Y~cum.time,data=data[[1]],xlab='Time',ylab='Incidence',col=1,xlim=c(0,max.time),ylim=c(0,max.y))
for (k in 1:nsims) {              #add multiple epidemics to plot
  lines(Y~cum.time,data=data[[k]],col=k,type='o')
}


