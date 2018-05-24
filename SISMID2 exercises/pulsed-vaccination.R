### R code from vignette source 'pulsed-vaccination.rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: pulsed-vaccination.rnw:76-93
###################################################
require(deSolve)                          #deSolve library needed for this computing session

sir.model.open <- function (t, x, params) {    #here we begin a function with three arguments
  S <- x[1]                               #create local variable S, the first element of x
  I <- x[2]                               #create local variable I
  R <- x[3]                               #create local variable R
  with(                                   #we can simplify code using "with"
       as.list(params),                   #this argument to "with" lets us use the variable names
       {                                  #the system of rate equations
         dS <- mu*(S+I+R) - beta*S*I - mu*S
         dI <- beta*S*I - gamma*I - mu*I
         dR <- gamma*I - mu*R
         dx <- c(dS,dI,dR)                #combine results into a single vector dx
         list(dx)                         #return result as a list
       }
       )
}


###################################################
### code chunk number 2: parameters
###################################################
R0 <- 10
N <-  1					                       #population size
mu <- 0.02                             #per capita birth/death rate
gamma <- 365/10  			                 #recovery rate (in years)
beta <- R0*(gamma+mu)/N 	             #transmission rate
xstart <- c(S=0.2, I=0.001, R=1-0.2-0.001)	 #initial conditions, must sum to one
Tmax <- 120                            #integrate for 200 years after transients
params <- c(beta=beta, gamma=gamma, mu=mu)           #parameter vector
tau <- 0.1                             #size of time step
times <- seq(0, Tmax, by=tau)          #function seq returns a sequence


###################################################
### code chunk number 3: pulsed-vaccination.rnw:113-114
###################################################
out <- ode(xstart,times,sir.model.open,params, method='ode45', rtol=1e-7)  


###################################################
### code chunk number 4: pulsed-vaccination.rnw:119-124
###################################################
op <- par(fig=c(0,0.5,0,1),mar=c(4,4,1,1))                      #set graphical parameters
plot(I~time,data=out, type='l', lwd=2)                          #plot the I variable against time
par(fig=c(0.5,1,0,1),mar=c(4,1,1,1),new=T)                      #re-set graphical parameters
plot(I~S,data=out,log='xy',yaxt='n',xlab='S', type='l', lwd=2)  #plot phase portrait
par(op)                                                         #re-set graphical parameters


###################################################
### code chunk number 5: reset-parameters
###################################################
xstart <- out[which(out[,1]==50),2:4]


###################################################
### code chunk number 6: new-parameters
###################################################
pv <- 0.1                     # fraction of susceptibles vaccomated
Tv <- 4                       # number of years between pulses
vacc.events <- floor(Tmax/Tv) # number of pulses in Tmax years


###################################################
### code chunk number 7: dataframe
###################################################
data <- data.frame(S=out[which(out[,1]==50),2],
                   I=out[which(out[,1]==50),3],
                   R=out[which(out[,1]==50),4])


###################################################
### code chunk number 8: pulsed-vaccination.rnw:155-162
###################################################
for(i in 1:vacc.events){
  out <- ode(xstart, seq(tau, Tv, by=tau), sir.model.open, params, method='ode45', rtol=1e-7)
  xstart <- out[dim(out)[1],2:4]        # reset initial condition
  xstart[1] <- (1-pv)*(tail(out,1)[2])  # vaccinate susceptibles
  xstart[3] <- xstart[3]+(pv)*(tail(out,1)[2])  # move to recovered class
  data <- rbind(data,out[,2:4])         # store result
}


###################################################
### code chunk number 9: plot-vax
###################################################
data$time <- seq(50, Tmax+50, by=tau)
par(mar=c(5,4,4,4)+0.1)
plot(data$time[1:500], data$I[1:500], type='l', xlab='Time', ylab='', col='red', axes=FALSE)
axis(2, col.axis='red')
mtext(side=2, line=2.5, 'Infected', col='red')
box()
axis(1)
par(new=TRUE)
plot(data$time[1:500], data$S[1:500], type='l', xlab='', ylab='', axes=FALSE, col='black')
axis(4)
mtext(side=4, line=2.5, 'Susceptibles')


