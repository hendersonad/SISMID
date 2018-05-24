## Exercise 3
## Estimating R0

require(deSolve)
setwd('/Users/AHenderson/OneDrive - London School of Hygiene and Tropical Medicine/PhD/SISMID/SISMID2 code')

test <- load('data/data.RData')     #load the data and plot flu cases
plot(flu,type='b',log='y',main='Epidemic in a British boarding school', cex.main=0.85,
     xlab='Day', ylab='Active influenza cases')


## R0 from final outbreak size
final.size.est <- function(R0, final.size){1 - final.size - exp(-final.size*R0)}

Z <- 512/764
x <- seq(0,3,0.01)
y <- final.size.est(x, Z)
plot(x,y, type='l')
abline(h=0, lty=2)
abline(v=1.655, lty=2)

####################################
## exercise 1
# plot relationship between R0 and final outbreak size for R(Inf) between 0 and 1

final.epi.size <- seq(0.001, 0.999, 0.001)
R0 <- c()
for(z in final.epi.size){
  R0 <- c(R0, log(1-z)/-z)
}
plot(final.epi.size, R0, type='l', xlab='Epidemic size')

####################################
####################################
####################################
## Linear approximation - first stage of epidemic log transform tells us R0
model <- lm(log(flu[1:4])~day[1:4],data=flu)
summary(model) 
slope <- coef(model)[2]
slope #slope = 1.09491

# so R0=slope/gamma + 1
## assume recovery rate = 2.5 days so gamma=1/2.5= 0.4 
R0 <- slope/(0.4) + 1 # = 3.7

####################################
## exercise 2
#   confined once symptomatic and no furhter infection. Assume confined within 24 hours
# so gamma  = 1/1 = 1
(R0.2 <- slope/(1) + 1 ) 
  # R0 = 2.1
# or within 12 hours gamma = 1/0.5 = 2
(R0.3 <- slope/(2) + 1 ) 
  # R0 = 1.5

####################################
## exercise 3
# measles in Niamey: infectious period is 2 weeks or 14/365= 0.0384 years
# estimate R0
head(niamey) # three communities 
plot(seq(0,15*2,by=2),niamey$V1,type='b',log='y',main='Measles in Niamey', cex.main=0.85,
     xlab='Epi.week', ylab='Measles cases')
par(new=T)
plot(seq(0,15*2,by=2),niamey$V2,type='b', log='y', yaxt='n', xlab='',ylab='',col=2)
par(new=T)
plot(seq(0,15*2,by=2),niamey$V3,type='b', log='y', yaxt='n', xlab='',ylab='', col=3)

#focus on community 1 and trajectory over 1st 18 data points
plot(seq(0,15*2,by=2),niamey$V1,type='b',log='y',main='Measles in Niamey', cex.main=0.85,
     xlab='Epi.week', ylab='Measles cases')

model <- lm(data=niamey, log(niamey[1:10,1])~seq(0,18,by=2))
slope <- coef(model)[2]

gamma <- 1/2 # infectious period = 2 weeks
R0 <- slope/(1/2) + 1 # 1.44

####################################
## exercise 4
slope <- NULL
se <- NULL
for (i in 3:18){
model <- lm(log(niamey[1:i,1])~seq(0,(i-1)*2,by=2))
slope <- c(slope, as.numeric(coef(model)[2]))
se <- c(se, summary(model)$coefficients[4])
}

R0 <- slope/0.5 + 1
plot(R0,se)
#so as you include more data R0 becomes more precise but more bias introduced

####################################
####################################
####################################
## Least squares estimation

#replace NA with zero
load('data/data.Rdata')
niamey[5,3] <- 0
niamey<-data.frame(biweek=rep(seq(1,16),3),site=c(rep(1,16),rep(2,16),rep(3,16)),
                   cases=c(niamey[,1],niamey[,2],niamey[,3])) #define "biweeks": data from wide to long format

plot(niamey$biweek,niamey$cases,type='p',col=niamey$site,xlab='Biweek',ylab='Cases')
lines(niamey$biweek[niamey$site==1],niamey$cases[niamey$site==1])
lines(niamey$biweek[niamey$site==2],niamey$cases[niamey$site==2],col=2)
lines(niamey$biweek[niamey$site==3],niamey$cases[niamey$site==3],col=3)

# the model
closed.sir.model1 <- function (t, x, params) {  #SIR model equations
  S <- x[1]
  I <- x[2]
  beta <- params[1]
  #gamma <- params[2]
  dS <- -beta*S*I
  dI <- beta*S*I-(365/13)*I # Removal rate is hard wired as 13 days
  list(c(dS,dI))
}

# the function
sse.func <- function(params0, data, site){
  data <- data[data$site==site,]
  t <- data[,1]*14/365 # time in biweeks
  cases <- data[,3]     # cases for that site
  beta <- exp(params0[1])  #parameter beta
  S0 <- exp(params0[2])     # parameter S0
  I0 <- exp(params0[3])     # parameter I0
  out <- as.data.frame(ode(c(S=S0, I=I0), times=t, closed.sir.model1, 
                    parms=c(beta=beta), hmax=1/120))
  sse <- sum((out$I - cases)^2)
}
# sneaky maths trick - we'll chuck the whole number list at this function
#so parameterise in terms of alternative vars log(beta) log(S0)...
#these map to model params on range 0 -> Inf - the range tha is biologicall meaningful

# the optimisation
require(deSolve)
params0 <- c(-3.2, # 0.04
              7.3,
              -2.6 #0.07
            )
fit1 <- optim(params0,sse.func, data=niamey, site=1, hessian=T)
exp(fit1$par)
fit2 <- optim(params0,sse.func, data=niamey, site=2, hessian=T)
exp(fit2$par)
fit3 <- optim(params0,sse.func, data=niamey, site=3, hessian=T)
exp(fit3$par)

# plot results
par(mfrow=c(2,2))
#site1
plot(cases~biweek,data=subset(niamey,site==1),type='p',pch=21, col=4) #plot site 1
t <- subset(niamey,site==1)[,1]*14/365 # time as biweeks
mod.pred<-as.data.frame(ode(c(S=exp(fit1$par[2]),I=exp(fit1$par[3])),times=t,closed.sir.model,
                              c(exp(fit1$par[1]),exp(fit1$par[4])),hmax=1/120)) # get model predicitons for site1
lines(mod.pred$I~subset(niamey,site==1)[,1]) #and plot as a line
#site2
plot(cases~biweek,data=subset(niamey,site==2),type='p',pch=21, col=2) #plot site 1
t <- subset(niamey,site==1)[,1]*14/365 # time as biweeks
mod.pred<-as.data.frame(ode(c(S=exp(fit2$par[2]),I=exp(fit2$par[3])),times=t,closed.sir.model,
                              c(exp(fit2$par[1]),exp(fit2$par[4])),hmax=1/120)) # get model predicitons for site1
lines(mod.pred$I~subset(niamey,site==2)[,1]) #and plot as a line
#site3
plot(cases~biweek,data=subset(niamey,site==3),type='p',pch=21, col=3) #plot site 1
t <- subset(niamey,site==1)[,1]*14/365 # time as biweeks
mod.pred<-as.data.frame(ode(c(S=exp(fit3$par[2]),I=exp(fit3$par[3])),times=t,closed.sir.model,
                              c(exp(fit3$par[1]),exp(fit3$par[4])),hmax=1/120)) # get model predicitons for site1
lines(mod.pred$I~subset(niamey,site==3)[,1]) #and plot as a line

####################################
## exercise 5
 #modify the code to estimate gamma at the same time...
closed.sir.model2 <- function (t, x, params) {  #SIR model equations
  S <- x[1]
  I <- x[2]
  beta <- params[1]
  gamma <- params[2]
  dS <- -beta*S*I
  dI <- beta*S*I-(gamma)*I # Removal rate is hard wired as 13 days
  list(c(dS,dI))
}
sse.func <- function(params0, data, site){
  data <- data[data$site==site,]
  t <- data[,1]*14/365 # time in biweeks
  cases <- data[,3]     # cases for that site
  beta <- exp(params0[1])  #parameter beta
  S0 <- exp(params0[2])     # parameter S0
  I0 <- exp(params0[3])     # parameter I0
  gamma <- exp(params0[4])     # parameter gamma
  out <- as.data.frame(ode(c(S=S0, I=I0), times=t, closed.sir.model2, 
                           parms=c(beta=beta, gamma=gamma), hmax=1/120))
  sse <- sum((out$I - cases)^2)
}
params0 <- c(-3.2, # 0.04
             7.3,
             -2.6, #0.07
             log(13/365)
)
fit4 <- optim(params0,sse.func, data=niamey, site=1, hessian=T)
exp(fit4$par)
fit5 <- optim(params0,sse.func, data=niamey, site=2, hessian=T)
exp(fit5$par)
fit6 <- optim(params0,sse.func, data=niamey, site=3, hessian=T)
exp(fit6$par)

par(mfrow=c(2,2))
#site1
plot(cases~biweek,data=subset(niamey,site==1),type='p',pch=21, col=4) #plot site 1
t <- subset(niamey,site==1)[,1]*14/365 # time as biweeks
mod.pred<-as.data.frame(ode(c(S=exp(fit1$par[2]),I=exp(fit1$par[3])),times=t,closed.sir.model1,
                            c(exp(fit1$par[1]),exp(fit1$par[4])),hmax=1/120)) # get model predicitons for site1
mod.pred2<-as.data.frame(ode(c(S=exp(fit4$par[2]),I=exp(fit4$par[3])),times=t,closed.sir.model2,
                            c(exp(fit4$par[1]),exp(fit4$par[4])),hmax=1/120)) # get model predicitons for site1
lines(mod.pred$I~subset(niamey,site==1)[,1]) #and plot as a line
lines(mod.pred2$I~subset(niamey,site==1)[,1],col=2) #and plot as a line
#site2
plot(cases~biweek,data=subset(niamey,site==2),type='p',pch=21, col=5) #plot site 1
t <- subset(niamey,site==3)[,1]*14/365 # time as biweeks
mod.pred<-as.data.frame(ode(c(S=exp(fit2$par[2]),I=exp(fit2$par[3])),times=t,closed.sir.model1,
                            c(exp(fit2$par[1]),exp(fit2$par[4])),hmax=1/120)) # get model predicitons for site1
mod.pred2<-as.data.frame(ode(c(S=exp(fit5$par[2]),I=exp(fit5$par[3])),times=t,closed.sir.model2,
                            c(exp(fit5$par[1]),exp(fit5$par[4])),hmax=1/120)) # get model predicitons for site1
lines(mod.pred$I~subset(niamey,site==2)[,1]) #and plot as a line
lines(mod.pred2$I~subset(niamey,site==2)[,1],col=2) #and plot as a line
#site3
plot(cases~biweek,data=subset(niamey,site==3),type='p',pch=21, col=3) #plot site 1
t <- subset(niamey,site==3)[,1]*14/365 # time as biweeks
mod.pred<-as.data.frame(ode(c(S=exp(fit3$par[2]),I=exp(fit3$par[3])),times=t,closed.sir.model1,
                            c(exp(fit3$par[1]),exp(fit3$par[4])),hmax=1/120)) # get model predicitons for site1
mod.pred2<-as.data.frame(ode(c(S=exp(fit6$par[2]),I=exp(fit6$par[3])),times=t,closed.sir.model2,
                            c(exp(fit6$par[1]),exp(fit6$par[4])),hmax=1/120)) # get model predicitons for site1
lines(mod.pred$I~subset(niamey,site==3)[,1]) #and plot as a line
lines(mod.pred2$I~subset(niamey,site==3)[,1],col=2) #and plot as a line

## basically they are very similar but a bit different if estimate gamma as well....

####################################
## exercise 6
# What happens if one or both of the other unknowns (I0 and S0) is fixed instead of gamma

#let I0 float and be fitted but gamma is fixes
load('data/data.Rdata')
niamey[5,3] <- 0
niamey<-data.frame(biweek=rep(seq(1,16),3),site=c(rep(1,16),rep(2,16),rep(3,16)),
                   cases=c(niamey[,1],niamey[,2],niamey[,3])) #define "biweeks": data from wide to long format

plot(niamey$biweek,niamey$cases,type='p',col=niamey$site,xlab='Biweek',ylab='Cases')
lines(niamey$biweek[niamey$site==1],niamey$cases[niamey$site==1])
lines(niamey$biweek[niamey$site==2],niamey$cases[niamey$site==2],col=2)
lines(niamey$biweek[niamey$site==3],niamey$cases[niamey$site==3],col=3)

# the model
closed.sir.model1 <- function (t, x, params) {  #SIR model equations
  S <- x[1]
  I <- x[2]
  beta <- params[1]
  gamma <- params[2]
  dS <- -beta*S*I
  dI <- beta*S*I-(gamma)*I # Removal rate is hard wired as 13 days
  list(c(dS,dI))
}
# the function
sse.func <- function(params0, data, site){
  data <- data[data$site==site,]
  t <- data[,1]*14/365 # time in biweeks
  cases <- data[,3]     # cases for that site
  beta <- exp(params0[1])  #parameter beta
  S0 <- 12500     # parameter S0
  I0 <- exp(params0[2])     # parameter I0
  gamma <- exp(params0[3])     # parameter I0
  out <- as.data.frame(ode(c(S=S0, I=I0), times=t, closed.sir.model1, 
                           parms=c(beta=beta, gamma=gamma), hmax=1/120))
  sse <- sum((out$I - cases)^2)
}
# the optimisation
params0 <- c(-3.2, # 0.04
             -2.6, #0.07
             log(13/365)
)

fit1 <- optim(params0,sse.func, data=niamey, site=1, hessian=T)
exp(fit1$par)
fit2 <- optim(params0,sse.func, data=niamey, site=2, hessian=T)
exp(fit2$par)
fit3 <- optim(params0,sse.func, data=niamey, site=3, hessian=T)
exp(fit3$par)

exp(fit1$par) # let I0 be estimated
exp(fit4$par)
exp(fit2$par) # let I0 be estimated
exp(fit5$par)
exp(fit3$par) # let I0 be estimated
exp(fit6$par)
