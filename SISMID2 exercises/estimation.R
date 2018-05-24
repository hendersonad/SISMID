### R code from vignette source 'estimation.rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: estimation.rnw:56-59
###################################################
load('data.RData')     #load the data and plot flu cases
plot(flu,type='b',log='y',main='Epidemic in a British boarding school', cex.main=0.85,
 xlab='Day', ylab='Active influenza cases')


###################################################
### code chunk number 2: estimation.rnw:101-105
###################################################
model<-lm(log(flu[1:4])~day[1:4],data=flu);  #fit a linear model
summary(model)         #summary statistics for fit model
slope<-coef(model)[2]  #extract slope parameter
slope                 #print to screen


###################################################
### code chunk number 3: estimation.rnw:127-131
###################################################
load('data.RData')
niamey[5,3]<-0  #replace a "NA"
niamey<-data.frame(biweek=rep(seq(1,16),3),site=c(rep(1,16),rep(2,16),rep(3,16)),
                   cases=c(niamey[,1],niamey[,2],niamey[,3])) #define "biweeks"


###################################################
### code chunk number 4: estimation.rnw:135-139
###################################################
plot(niamey$biweek,niamey$cases,type='p',col=niamey$site,xlab='Biweek',ylab='Cases')
lines(niamey$biweek[niamey$site==1],niamey$cases[niamey$site==1])
lines(niamey$biweek[niamey$site==2],niamey$cases[niamey$site==2],col=2)
lines(niamey$biweek[niamey$site==3],niamey$cases[niamey$site==3],col=3)


###################################################
### code chunk number 5: estimation.rnw:147-155
###################################################
closed.sir.model <- function (t, x, params) {  #SIR model equations
  S <- x[1]
  I <- x[2]
  beta <- params
  dS <- -beta*S*I
  dI <- beta*S*I-(365/13)*I
  list(c(dS,dI))
}


###################################################
### code chunk number 6: estimation.rnw:162-174
###################################################

sse.sir <- function(params0,data,site){  #function to calculate squared errors
  data<-data[data$site==site,]    #working dataset, based on site
  t <- data[,1]*14/365            #time in biweeks
  cases <- data[,3]               #number of cases
  beta <- exp(params0[1])            #parameter beta
  S0 <- exp(params0[2])           #initial susceptibles
  I0 <- exp(params0[3])           #initial infected        
  out <- as.data.frame(ode(c(S=S0,I=I0),times=t,closed.sir.model,beta,hmax=1/120))
  sse<-sum((out$I-cases)^2)       #sum of squared errors
}



###################################################
### code chunk number 7: estimation.rnw:183-192
###################################################
library(deSolve)   #differential equation library
params0<-c(-3.2,7.3,-2.6)  #initial guess

fit1 <- optim(params0,sse.sir,data=niamey,site=1) #fit
exp(fit1$par)  #back-transform parameters
fit2 <- optim(params0,sse.sir,data=niamey,site=2) #fit
exp(fit2$par)  #back-transform parameters
fit3 <- optim(params0,sse.sir,data=niamey,site=3) #fit
exp(fit3$par)  #back-transform parameters


###################################################
### code chunk number 8: estimation.rnw:197-219
###################################################

par(mfrow=c(3,1))   #set up plotting area for multiple panels
plot(cases~biweek,data=subset(niamey,site==1),type='b',col='blue', pch=21) #plot site 1
t <- subset(niamey,site==1)[,1]*14/365
mod.pred<-as.data.frame(ode(c(S=exp(fit1$par[2]),I=exp(fit1$par[3])),times=t,
                              closed.sir.model,exp(fit1$par[1]),hmax=1/120))
                              #obtain model predictions
lines(mod.pred$I~subset(niamey,site==1)[,1]) #and plot as a line

plot(cases~biweek,data=subset(niamey,site==2),type='b',col=site) #site 2
t <- subset(niamey,site==2)[,1]*14/365
mod.pred<-as.data.frame(ode(c(S=exp(fit2$par[2]),I=exp(fit2$par[3])),times=t,
                              closed.sir.model,exp(fit2$par[1]),hmax=1/120))
lines(mod.pred$I~subset(niamey,site==2)[,1])


plot(cases~biweek,data=subset(niamey,site==3),type='b',col=site) #site 3
t <- subset(niamey,site==3)[,1]*14/365
mod.pred<-as.data.frame(ode(c(S=exp(fit3$par[2]),I=exp(fit3$par[3])),times=t,
                              closed.sir.model,exp(fit3$par[1]),hmax=1/120))
lines(mod.pred$I~subset(niamey,site==3)[,1])



