### R code from vignette source 'structured-models.rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: structured-models.rnw:54-55
###################################################
require(deSolve)


###################################################
### code chunk number 2: structured-models.rnw:63-65
###################################################
load('data.RData')    #load the data
plot(measles$Time,measles$Cases,type='l', xlab='Year',ylab='Cases')  #plot cases over time


###################################################
### code chunk number 3: structured-models.rnw:72-78
###################################################
#plot measles cases in subsequent two year intervals using the constructed variable "TwoYEAR"
plot(measles$TwoYear,measles$Cases,type='p',pch=20,col='grey',xlab='Time (years)',ylab='Cases')
#fit a smooth line using loess -- notice data must be ordered for loess to fit properly
smooth.cases<-loess(measles$Cases[order(measles$TwoYear)]~measles$TwoYear[order(measles$TwoYear)],
                    span=0.3)
lines(smooth.cases$x,smooth.cases$fitted,lwd=2)  #add smooth fit


###################################################
### code chunk number 4: structured-models.rnw:86-99
###################################################
age.model<-function(t,x,parms){  #a function to return derivatives of age structured model
 S<-x[1:4]     #S are the first four elements of x
 E<-x[5:8]     #E are the next four elements of x
 I<-x[9:12]    #I are the last four elements ofx
 dx<-vector(length=12)   #a vector to store the derivatives
  for(a in 1:4){   #loop over age classes
   tmp <- (parms$beta[a,]%*%I)*S[a]   #temporary variable with infection rate
   dx[a] <- parms$nu[a]*55/75 - tmp - parms$mu[a]*S[a]                 #dS
   dx[a+4] <- tmp - parms$sigma*E[a] - parms$mu[a]*E[a]                #dE
   dx[a+8] <- parms$sigma*E[a] - parms$gamma*I[a] - parms$mu[a]*I[a]   #dI
 }  
return(list(dx))  #return the result
}


###################################################
### code chunk number 5: structured-models.rnw:104-112
###################################################
y0<-c(0.05, 0.01, 0.01, 0.008, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001)
#initialize state variables

#a list of model parameters
parms<-list(beta=matrix(c(2.089, 2.089, 2.086, 2.037, 2.089, 9.336, 2.086, 2.037, 2.086, 2.086,
            2.086, 2.037, 2.037, 2.037, 2.037,2.037),nrow=4,byrow=TRUE),
            sigma=1/8, gamma=1/5, nu=c(1/(55*365),0,0,0),  mu=c(1/(55*365),0,0,0))
parms


###################################################
### code chunk number 6: structured-models.rnw:117-146
###################################################
n=c(6,4,10,55)/75   #number of years in each age class
maxTime <- 100*365  #number of days in 100 years
T0=0                #initial time
S=c()               #initialize S
E=c()               #initialize E
I=c()               #initialize E
T=c()               #initialize T, a vector to hold times
while(T0<maxTime){  #loop over times
  y=lsoda(y0,c(T0, T0+365),age.model,parms) #solve diff'l equation for each time
  T=rbind(T, y[2,1])      #store results
  S=rbind(S, y[2,2:5])
  E=rbind(E, y[2,6:9])
  I=rbind(I, y[2,10:13])
  #Now do the yearly movements
  #Note use of "tail" to pull off the last value in a vector
  y0[1]=tail(y,1)[2]-tail(y,1)[2]/6
  y0[2]=tail(y,1)[3]+tail(y,1)[2]/6 - tail(y,1)[3]/4
  y0[3]=tail(y,1)[4]+tail(y,1)[3]/4 - tail(y,1)[4]/10
  y0[4]=tail(y,1)[5]+tail(y,1)[4]/10
  y0[5]=tail(y,1)[6]-tail(y,1)[6]/6
  y0[6]=tail(y,1)[7]+tail(y,1)[6]/6 - tail(y,1)[7]/4
  y0[7]=tail(y,1)[8]+tail(y,1)[7]/4 - tail(y,1)[8]/10
  y0[8]=tail(y,1)[9]+tail(y,1)[8]/10
  y0[9]=tail(y,1)[10]-tail(y,1)[10]/6
  y0[10]=tail(y,1)[11]+tail(y,1)[10]/6 - tail(y,1)[11]/4
  y0[11]=tail(y,1)[12]+tail(y,1)[11]/4 - tail(y,1)[12]/10
  y0[12]=tail(y,1)[13]+tail(y,1)[12]/10
  T0=tail(T,1)
}


###################################################
### code chunk number 7: structured-models.rnw:148-163
###################################################

#plot
par(mfrow=c(2,1))    #set up plotting region
plot(T,S[,1],type='l',xlim=c(0,45000),ylim=c(0,0.06),xlab='Time (days)',
     ylab='Proportion susceptible') #plot susceptibles in youngest age class
lines(T,S[,2],col='blue')           #susceptibles in second age class
lines(T,S[,3],col='red')            #susceptibles in third age class
lines(T,S[,4],col='green')          #susceptibles in oldest age class
legend(x='topright',legend=c('<5','6-9','10-19','20+'),  #add legent
     col=c('black','blue','red','green'),lty=1,bty='n')
plot(T,I[,1],type='l',log='y',xlim=c(0,45000),xlab='Time (days)', #plot infected
     ylab='Proportion infected')
lines(T,I[,2],col='blue')
lines(T,I[,3],col='red')
lines(T,I[,4],col='green')


