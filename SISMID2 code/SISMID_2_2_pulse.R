## Exercise 2
## Pulsed vaccinations

require(deSolve)

sir.model.open <- function(t,x,params){
  S <- x[1]
  I <- x[2]
  R <- x[3]
  with(as.list(params),{
    dS = mu*(S+I+R) - beta*S*I - mu*S
    dI = beta*S*I - I*(mu + gamma)
    dR = gamma*I - mu*R
    dx <- c(dS,dI,dR)
    list(dx)
  }
  )
}

# define model parameters
R0 <- 10; N <- 1; Tmax=120; tau = 0.1;
params <- c(
    beta=NA,
    gamma=365/10, # recovery rate in years
    mu=0.02 # mu = per capita p.a. birth/death rate
)
params[["beta"]] <- R0*(params[["gamma"]]+params[["mu"]])/N #annual t-rate
xstart <- c(S=0.2, I=0.001, R=1-0.2-0.001)
times <- seq(0, Tmax, tau)

# solve model
out <- ode(xstart,times,sir.model.open,params, method='ode45', rtol=1e-7)

op <- par(fig=c(0,0.5,0,1),mar=c(4,4,1,1)) #set graphical parameters
plot(I~time,data=out, type='l', lwd=2) #plot the I variable against time
par(fig=c(0.5,1,0,1),mar=c(4,1,1,1),new=T) #re-set graphical parameters
plot(I~S,data=out,log='xy',yaxt='n',xlab='S', type='l', lwd=2) #plot phase portrait
par(op)


# so the system as settled to equilibrium around t=50 so use this as initialisation of pulse vaccination
xstart <- out[which(out[,1]==50),2:4] #S (roughly)= 0.1, R= 0.9 and I=0

##################################
##################################
## adding pulsed vaccinations 
pv <- 0.1 #fraction of S you can vaccinate per pulse
Tv <- 4   #period between pulses 
vacc.events <- floor(Tmax/Tv) #number of pulses in Tmax years = 30

data <- data.frame(S=out[which(out[,1]==50),2],
                   I=out[which(out[,1]==50),3],
                   R=out[which(out[,1]==50),4]
                   )
# now use terminal conditions of transient period as initial conditions
# move forward Tv units of time
# vaccinate pv (move from S to R)
# repeat until Tmax

for(i in 1:vacc.events){
  out <- ode(xstart,seq(tau, Tv, by=tau),sir.model.open,params, method='ode45', rtol=1e-7)
  xstart <- out[dim(out)[1],2:4] #reset initial conditions
  xstart[1] <- (1-pv)*(tail(out,1)[2])  #move the vaccinated out of S
  xstart[3] <- xstart[3]+(pv)*(tail(out,1)[2]) #move the vaccinated into R
  data <- rbind(data,out[,2:4]) # store result
}

par(mfrow=c(1,1))
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


#Exercise 1. Why might this periodicity be of interest? What might be the consequences of \peaks" and
# troughs" for public health?
# peaks of lots of cases could be heavier disease burden on health system?
# epidemic troughs can also be utilised for elimination campaigns

#By modifying the value of pv, see if you can locate the vaccination threshold. Does this
# agree with the analytically predicted threshold
mu <- params[["mu"]]
crit <- function(pv,mu,Tv,R0){
  ((mu*Tv-pv)*(exp(mu*Tv)-1)+mu*pv*Tv)/(mu*Tv*(pv-1+exp(mu*Tv)))-(1/R0)
}
  
plot(seq(0,0.6,by=0.01), crit(seq(0,0.6,by=0.01), mu, Tv, R0), type='l', xlab='p_v', ylab='Elimination > abline(h=0, lty=2)')
abline(h=0, lty=2)

x <- seq(0,0.6,0.01)
thresh <- crit(x, mu, Tv, R0)

thresholds <- x[thresh<0]
thresholds[1] # = 0.54

# or if you're more smart:
print(vax <- uniroot(crit, lower=0, upper=0.6, mu=mu, Tv=Tv, R0=R0))
  #root = 0.5312612 = vacc threshold

# Exercise 3. One thing we might be interested in knowing is how our vaccination strategy affects mean
#disease prevalance. Devise a strategy to study this question and produce a plot illustrating the result.

# cycle through a range of plausible pv values: (run exercise 1 again first)
# define model parameters
R0 <- 10; N <- 1; Tmax=120; tau = 0.1;
params <- c(
  beta=NA,
  gamma=365/10, # recovery rate in years
  mu=0.02 # mu = per capita p.a. birth/death rate
)
params[["beta"]] <- R0*(params[["gamma"]]+params[["mu"]])/N #annual t-rate
xstart <- c(S=0.2, I=0.001, R=1-0.2-0.001)
times <- seq(0, Tmax, tau)
# solve model
out <- ode(xstart,times,sir.model.open,params, method='ode45', rtol=1e-7)
# so the system as settled to equilibrium around t=50 so use this as initialisation of pulse vaccination
data.critical <- data.frame(pv=numeric(0), I=numeric(0))
for (pv in seq(0,0.6,0.01)){
xstart <- out[which(out[,1]==50),2:4] #S (roughly)= 0.1, R= 0.9 and I=0
Tv <- 4   #period between pulses 
vacc.events <- floor(Tmax/Tv) #number of pulses in Tmax years = 30

data <- data.frame(S=out[which(out[,1]==50),2],
                   I=out[which(out[,1]==50),3],
                   R=out[which(out[,1]==50),4]
)
for(i in 1:vacc.events){
  out2 <- ode(xstart,seq(tau, Tv, by=tau),sir.model.open,params, method='ode45', rtol=1e-7)
  xstart <- out2[dim(out2)[1],2:4] #reset initial conditions
  xstart[1] <- (1-pv)*(tail(out2,1)[2])  #move the vaccinated out2 of S
  xstart[3] <- xstart[3]+(pv)*(tail(out2,1)[2]) #move the vaccinated into R
  data <- rbind(data,out2[,2:4]) # store result
}
data.critical <- rbind(data.critical, data.frame(pv, I=mean(tail(data$I,1000))))
}
head(data.critical)
tail(data.critical)

plot(data.critical$pv,data.critical$I, type='l', xlab='pv', ylab='Mean infected')
vax <- uniroot(crit, lower=0, upper=0.6, mu=mu, Tv=Tv, R0=R0)
abline(v=vax$root, lty=2)

#### make the fancy plot sometime...

