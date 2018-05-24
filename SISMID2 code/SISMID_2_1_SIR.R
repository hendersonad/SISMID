## Exercise 1 
## Basic SIR model

# demographically closed SIR model 

require(deSolve)

# Define the ODE model
sir.model.closed <- function(t, x, params){ #t=time, x=init.state, params=parameters
  S <- x[1]
  I <- x[2]
  R <- x[3]
  with(
    as.list(params),
    {
      dS <- -beta*S*I
      dI <- beta*S*I - gamma*I
      dR <- gamma*I
      dx <- c(dS,dI,dR)
      list(dx)
    }
  )
}

# Set initial conditions
#   time length
times <- seq(0,120,5)
#   list of parameters
params <- c(
  beta=0.3,
  gamma=1/7
)
#   Initial conditions - PROPORTIONS not NUMBERS
xstart <- c(
  S=9999/10000,
  I=1/10000,
  R=0
)

# solve the ODEs 
out <- as.data.frame(ode(xstart,times,sir.model.closed,params))

# plot the results
op <- par(fig=c(0,0.5,0,1),mar=c(4,4,1,1)) # set graph parameters
plot(data=out,I~time,type='b')                  ## Plot infected over time 
par(fig=c(0.5,1,0,1),mar=c(4,1,1,1),new=T) # set new graph parameters
plot(data=out,I~S,type='b',yaxt='n',xlab='S')   ## Plot infected vs susceptible ("phase portrait")
par(op)

########################################################################
########################################################################
## Explore dynamics over different beta and gamma values
beta_seq <- c(0.1, 1, 10, 100)
gamma_seq <- c(1, 5, 10, 50)
# time series
par(mfrow=c(4,4))
for(beta in beta_seq){
  params[["beta"]] <- beta
  
  for(gamma in gamma_seq){
  params[["gamma"]] <- 1/gamma
  
  out <- as.data.frame(ode(xstart,times,sir.model.closed,params))
  plot(data=out,I~time,type='b',xlab='S', main=paste0("gamma=",gamma, " beta=",beta))
  }
}

# phase portrait
par(mfrow=c(4,4))
for(beta in beta_seq){
  params[["beta"]] <- beta
  
  for(gamma in gamma_seq){
  params[["gamma"]] <- 1/gamma
  
  out <- as.data.frame(ode(xstart,times,sir.model.closed,params))
  plot(data=out,I~S,type='b',xlab='S', main=paste0("gamma=",gamma, " beta=",beta))
  }
}


## Explore dynamics over different initial state values
par(mfrow=c(1,3))
params <- c(beta=0.3, gamma=1/7)
xstart <- xstart <- c( # 1 infected
  S=9999/10000,
  I=1/10000,
  R=0
)
out <- as.data.frame(ode(xstart,times,sir.model.closed,params))
plot(data=out,I~time,type='b',xlab='time', main="1 infected")


xstart <- xstart <- c( # 100 infected
  S=9990/10000,
  I=100/10000,
  R=0
)
out <- as.data.frame(ode(xstart,times,sir.model.closed,params))
plot(data=out,I~time,type='b',xlab='time', main="100 infected")


xstart <- xstart <- c( # 1000 infected
  S=9000/10000,
  I=1000/10000,
  R=0
)
out <- as.data.frame(ode(xstart,times,sir.model.closed,params))
plot(data=out,I~time,type='b',xlab='time', main="10000 infected")


########################################################################
########################################################################
# Change to demographically open model
  
  sir.model.open <- function(t, x, params){ #t=time, x=init.state, params=parameters
    S <- x[1]
    I <- x[2]
    R <- x[3]
    with(
      as.list(params),
      {
        dS <- -beta*S*I + mu*(S+I+R) - mu*S
        dI <- beta*S*I - gamma*I - mu*I
        dR <- gamma*I - mu*R
        dx <- c(dS,dI,dR)
        list(dx)
      }
    )
  }
  times <- seq(0,120, 1)
  xstart <- xstart <- c( # 100 infected
    S=9999/10000,
    I=1/10000,
    R=0
  )
  # birth rate 1/(365*60) <- assumes lifespan of 60 years and 1 offspring per person
  params1 <- c(beta=0.3, gamma=1/7, mu=1/52)
  params2 <- c(beta=0.3, gamma=1/7, mu=1/520)
  params3 <- c(beta=0.3, gamma=1/7, mu=1/5200)
  
  out1 <- as.data.frame(ode(xstart,times,sir.model.open,params1))
  out2 <- as.data.frame(ode(xstart,times,sir.model.open,params2))
  out3 <- as.data.frame(ode(xstart,times,sir.model.open,params3))
  
  op <- par(fig=c(0,0.5,0,1),mar=c(4,4,1,1)) # set graph parameters
  plot(data=out1,I~time,type='l', col=4, ylim=c(0, 0.3), lwd=2)
  lines(data=out2,I~time, type='l', col=2, lwd=2)
  lines(data=out3,I~time, type='l', col=3, lwd=2)
  par(fig=c(0.5,1,0,1),mar=c(4,1,1,1),new=T) # set new graph parameters
  plot(data=out1,I~S,type='l', col=4, log='x',yaxt='n',xlab='S', ylim=c(0.01, 0.3), xlim=c(0.01, 1), lwd=2)
  lines(data=out2,I~S, type='l', col=2, lwd=2)
  lines(data=out3,I~S, type='l', col=3, lwd=2)
  legend('topleft', col=c('red','blue','green'), lty=1,
          legend=c('short lifespan', 'intermediate lifespan', 'long lifespan'))
  par(op)

########################################################################
########################################################################
# Change to SEIR model
seir.model.open <- function(t, x, params){ 
  S <- x[1]
  E <- x[2]
  I <- x[3]
  R <- x[4]
  with(
    as.list(params),
    {
      dS <- -beta*S*I + mu*(S+E+I+R) - mu*S
      dE <- beta*S*I - sigma*E - mu*E
      dI <- sigma*E - gamma*I - mu*I
      dR <- gamma*I - mu*R
      dx <- c(dS,dE,dI,dR)
      list(dx)
    }
  )
}
  
seir.model.closed <- function(t, x, params){ 
  S <- x[1]
  E <- x[2]
  I <- x[3]
  R <- x[4]
  with(
    as.list(params),
    {
      dS <- -beta*S*I 
      dE <- beta*S*I - sigma*E
      dI <- sigma*E - gamma*I
      dR <- gamma*I 
      dx <- c(dS,dE,dI,dR)
      list(dx)
    }
  )
}
times <- seq(0,240, 1)
xstart <- xstart <- c( # 100 infected
  S=9999/10000,
  E=0,
  I=1/10000,
  R=0
)

params1 <- c(beta=0.3, gamma=1/7, sigma=1/0.0001, mu=1/52)
params2 <- c(beta=0.3, gamma=1/7, sigma=1/7, mu=1/52)
params3 <- c(beta=0.3, gamma=1/7, sigma=1/14, mu=1/52)

out1 <- as.data.frame(ode(xstart,times,seir.model.closed,params1))
out2 <- as.data.frame(ode(xstart,times,seir.model.closed,params2))
out3 <- as.data.frame(ode(xstart,times,seir.model.closed,params3))

op <- par(fig=c(0,0.5,0,1),mar=c(4,4,1,1)) # set graph parameters
plot(data=out1,I~time,type='l', col=4, ylim=c(0, 0.3), lwd=2)
lines(data=out2,I~time, type='l', col=2, lwd=2)
lines(data=out3,I~time, type='l', col=3, lwd=2)
par(fig=c(0.5,1,0,1),mar=c(4,1,1,1),new=T) # set new graph parameters
plot(data=out1,I~S,type='l', col=4, log='x',yaxt='n',xlab='S', ylim=c(0.01, 0.3), xlim=c(0.01, 1), lwd=2)
lines(data=out2,I~S, type='l', col=2, lwd=2)
lines(data=out3,I~S, type='l', col=3, lwd=2)
legend('topleft', col=c('red','blue','green'), lty=1,
       legend=c('short lifespan', 'intermediate lifespan', 'long lifespan'))
par(op)