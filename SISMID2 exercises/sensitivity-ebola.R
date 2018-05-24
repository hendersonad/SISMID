### R code from vignette source 'sensitivity-ebola.rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: sensitivity-ebola.rnw:86-121
###################################################
legrand <- function(t, x, params){

  # states
  S <- x[1]
  E <- x[2]
  I <- x[3]
  H <- x[4]
  F <- x[5]
  R <- x[6]
  N <- S+E+I+H+F+R
  
  # calculated parameters
  gammaih = params$gammai*params$gammah/(params$gammah - params$gammai)
  gammadh = params$gammad*params$gammah/(params$gammah - params$gammad)
  # rates
  r1 <- (params$betaI*S*I + params$betaH*S*H + params$betaF*S*F)/N
  r2 <- params$alpha*E
  r3 <- params$gammah*params$theta1*I
  r4 <- gammadh*params$delta2*H
  r5 <- params$gammaf*F
  r6 <- params$gammai*(1-params$theta1)*(1-params$delta1)*I
  r7 <- params$delta1*(1-params$theta1)*params$gammad*I
  r8 <- gammaih*(1-params$delta2)*H
  
  # derivatives
  dS <- -r1
  dE <- r1 - r2
  dI <- r2 - r3 - r7 - r6
  dH <- r3 - r4 - r8
  dF <- r4 + r7 - r5
  dR <- r5 + r6 + r8
    
  #output
  out <- list(c(dS, dE, dI, dH, dF, dR))
}


###################################################
### code chunk number 2: sensitivity-ebola.rnw:125-143
###################################################
require(deSolve)

times <- seq(0, 52, by=1)  #solve for 52 weeks
params <- list(betaI=2.532,
               betaH=0.012,
               betaF=0.462,
               alpha=7/12,
               gammah=7/4.2,
               gammai=7/10,
               gammaf=7/2,
               gammad = 7/8,
               theta1=0.65,
               delta1=0.47,
               delta2=0.42)
pop.size <- 470000
I0 <- 9
xstart <- c(S=pop.size-I0, E=0, I=I0 , H=0, F=0, R=0)
out <- as.data.frame(ode(xstart, times, legrand, params))


###################################################
### code chunk number 3: sensitivity-ebola.rnw:147-156
###################################################
plot.legrand <- function(out){
  par(mfrow=c(3,2))
  plot(out$time, out$S, xlab='Time', ylab='Susceptible', type='l', col='steelblue4', lwd=3)
  plot(out$time, out$E, xlab='Time', ylab='Exposed', type='l', col='steelblue4', lwd=3)
  plot(out$time, out$I, xlab='Time', ylab='Infectious', type='l', col='steelblue4', lwd=3)
  plot(out$time, out$H, xlab='Time', ylab='Hospitalized', type='l', col='steelblue4', lwd=3)
  plot(out$time, out$F, xlab='Time', ylab='Funeral', type='l', col='steelblue4', lwd=3)
  plot(out$time, out$R, xlab='Time', ylab='Removed', type='l', col='steelblue4', lwd=3)
}


###################################################
### code chunk number 4: sensitivity-ebola.rnw:159-160
###################################################
plot.legrand(out)


###################################################
### code chunk number 5: sensitivity-ebola.rnw:169-207
###################################################

# model of Legrand et al (2007) with interventions
legrand2 <- function(t, x, params){

  # states
  S <- x[1]
  E <- x[2]
  I <- x[3]
  H <- x[4]
  F <- x[5]
  R <- x[6]
  N <- S+E+I+H+F+R
  
  # calculated parameters
  gammaih = params$gammai*params$gammah/(params$gammah - params$gammai)
  gammadh = params$gammad*params$gammah/(params$gammah - params$gammad)
  
  # rates
  r1 <- (params$betaI*S*I*(t<T) + params$betaI*S*I*(1-params$z)*(t>=params$T) + params$betaH*S*H*(t<params$T) + params$betaF*S*F*(t<params$T))/N
  r2 <- params$alpha*E
  r3 <- params$gammah*params$theta1*I
  r4 <- gammadh*params$delta2*H
  r5 <- params$gammaf*F
  r6 <- params$gammai*(1-params$theta1)*(1-params$delta1)*I
  r7 <- params$delta1*(1-params$theta1)*params$gammad*I
  r8 <- gammaih*(1-params$delta2)*H
  
  # derivatives
  dS <- -r1
  dE <- r1 - r2
  dI <- r2 - r3 - r7 - r6
  dH <- r3 - r4 - r8
  dF <- r4 + r7 - r5
  dR <- r5 + r6 + r8
    
  #output
  out <- list(c(dS, dE, dI, dH, dF, dR))
}


###################################################
### code chunk number 6: sensitivity-ebola.rnw:211-228
###################################################
params2 <- list(betaI=2.532,
               betaH=0.012,
               betaF=0.462,
               alpha=7/12,
               gammah=7/4.2,
               gammai=7/10,
               gammaf=7/2,
               gammad = 7/8,
               theta1=0.65,
               delta1=0.47,
               delta2=0.42,
               T=7,
               z=0.88)
pop.size <- 470000
I0 <- 9
xstart <- c(S=pop.size-I0, E=0, I=I0 , H=0, F=0, R=0)
out2 <- as.data.frame(ode(xstart, times, legrand2, params2))


###################################################
### code chunk number 7: sensitivity-ebola.rnw:232-233
###################################################
plot.legrand(out2)


###################################################
### code chunk number 8: sensitivity-ebola.rnw:240-241
###################################################
get.size <- function(out) tail(out$R,1)


###################################################
### code chunk number 9: sensitivity-ebola.rnw:247-256
###################################################
output0 <- rep(NA, 15)
intervention.times <- seq(1,15)
for(T in intervention.times){
  params2$T <- T
  out3 <- as.data.frame(ode(xstart, times, legrand2, params2))
  output0[T] <- get.size(out3)  
}

plot(intervention.times, output0, log='y', xlab='Intervention time', ylab='Epidemic size', las=1)


###################################################
### code chunk number 10: sensitivity-ebola.rnw:264-268
###################################################
require(lhs)            #add the lhs library
h <- 1000               #choose number of points
set.seed(6242015)
lhs<-maximinLHS(h,12)   #simulate


###################################################
### code chunk number 11: sensitivity-ebola.rnw:273-297
###################################################
betaI.min <- 1
betaI.max <- 4
betaH.min <- 0.01
betaH.max <- 0.5
betaF.min <- 0.1
betaF.max <- 4
alpha.min <- 7/21
alpha.max <- 7/2
gammah.min <- 7/2
gammah.max <- 7/1
gammai.min <- 7/21
gammai.max <- 7/2
gammaf.min <- 7/7
gammaf.max <- 7/1
gammad.min <- 7/21
gammad.max <- 7/2
theta1.min <- 0
theta1.max <- 1
delta1.min <- 0
delta1.max <- 1
delta2.min <- 0
delta2.max <- 1
z.min <- 0
z.max <- 1


###################################################
### code chunk number 12: sensitivity-ebola.rnw:301-314
###################################################
params.set <- cbind(
  betaI = lhs[,1]*(betaI.max-betaI.min)+betaI.min,
  betaH = lhs[,2]*(betaH.max-betaH.min)+betaH.min,
  betaF = lhs[,3]*(betaF.max-betaF.min)+betaF.min,
  alpha = lhs[,4]*(alpha.max-alpha.min)+alpha.min,
  gammah = lhs[,5]*(gammah.max-gammah.min)+gammah.min,
  gammai = lhs[,6]*(gammai.max-gammai.min)+gammai.min,
  gammaf = lhs[,7]*(gammaf.max-gammaf.min)+gammaf.min,
  gammad = lhs[,8]*(gammad.max-gammad.min)+gammad.min,  
  theta1 = lhs[,9]*(theta1.max-theta1.min)+theta1.min,  
  delta1 = lhs[,10]*(delta1.max-delta1.min)+delta1.min,  
  delta2 = lhs[,11]*(delta2.max-delta2.min)+delta2.min,    
  z = lhs[,12]*(z.max-z.min)+z.min)  


###################################################
### code chunk number 13: sensitivity-ebola.rnw:319-320
###################################################
levels <- 15


###################################################
### code chunk number 14: sensitivity-ebola.rnw:325-326
###################################################
h2 <-250


###################################################
### code chunk number 15: sensitivity-ebola.rnw:331-345
###################################################
j <- 1  
data <- data.frame(matrix(rep(NA,levels*h2*14),nrow=levels*h2))
for(i in 1:h2){
  for (T in intervention.times){
    
    data[j,1:13] <- params <- as.list(c(params.set[i,], T=T))
    out <- as.data.frame(ode(xstart, times, legrand2, params))
    data[j,14] <- get.size(out)
    j <- j+1
    
  }
}
names(data) <- c(names(params),'outbreak.size')
save(data, file='data.Rdata')


###################################################
### code chunk number 16: sensitivity-ebola.rnw:351-356
###################################################
load('data.Rdata')
plot(intervention.times, output0, type='l', lwd=3, ylim=c(10,2e6), log='y',
            xlab='Intervention times',
            ylab='Epidemic size')
points(data$T, data$outbreak.size, pch=19, cex=0.3, col='blue')


###################################################
### code chunk number 17: sensitivity-ebola.rnw:363-366
###################################################
boxplot(data$outbreak.size~data$T, ylim=c(10,2e6), border='blue', log='y', pch='.')
par(new=TRUE)
plot(intervention.times, output0, type='l', lwd=3, ylim=c(10,2e6), log='y', xlab='Intervention times', ylab='Epidemic size', axes=FALSE)


###################################################
### code chunk number 18: sensitivity-ebola.rnw:382-386
###################################################
library(sensitivity)
bonferroni.alpha <- 0.05/12
prcc <- pcc(data[,1:12], data[,14], nboot = 1000, rank=TRUE, conf=1-bonferroni.alpha)
save(prcc, file='prcc.Rdata')


###################################################
### code chunk number 19: sensitivity-ebola.rnw:392-405
###################################################
library(sensitivity)
load('prcc.Rdata')
summary <- print(prcc)
par(mar=c(7,4,4,2)+0.1)
plot(summary$original, main='Partial rank correlation coefficients', ylim=c(-1,1),
     xlab='', ylab='Coefficient',
     axes=FALSE)
axis(2)
axis(1, at=seq(1:12), labels=row.names(summary), las=2)
mtext(text='Parameter', side=1, line=4.5)
box()
for(i in 1:12) lines(c(i,i),c(summary[i,4], summary[i,5]))
abline(h=0)


