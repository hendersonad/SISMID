### R code from vignette source 'sensitivity.rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: sensitivity.rnw:91-95
###################################################
lambdaS <- function(betaSU,YSU,betaST,YST,pSU,YRU,pST,YRT,N){
  (betaSU*YSU+betaST*YST+pSU*betaSU*YRU+pST*betaST*YRT)/N
}
lambdaR <- function(betaRU,YRU,betaRT,YRT,N) (betaRU*YRU+betaRT*YRT)/N


###################################################
### code chunk number 2: sensitivity.rnw:113-137
###################################################
hiv.model.sf <- function(t, x, params){
  X <- x[1]
  Y <- x[2]
  YSU <- x[3]
  YST <- x[4]
  YRU <- x[5]
  YRT <- x[6]
  N <- X+YSU+YST+YRU+YRT
  betaST <- params$alpha*params$betaSU
  betaRU <- params$alpha2*params$betaSU
  betaRT <- params$alpha3*betaRU
  sigmaS <- params$FS*(params$nuSU+params$mu)/(1-params$FS)
  sigmaR <- params$FR*(params$nuRU+params$mu+params$q)/(1-params$FR)
  lambdaS <- as.numeric(lambdaS(params$betaSU,YSU,betaST,YST,params$pSU,YRU,params$pST,YRT,N))
  lambdaR <- as.numeric(lambdaR(betaRU,YRU,betaRT,YRT,N))
  dX <- params$pi - (params$c*(lambdaS+lambdaR)+params$mu)*X
  dY <- (params$c*(lambdaS+lambdaR)+params$mu)*X
  dYSU <- lambdaS*params$c*X + params$q*YRU + params$gS*YST - YSU*(sigmaS+params$nuSU+params$mu)
  dYST <- YSU*sigmaS - YST*(params$gS + params$r + params$nuST + params$mu)  
  dYRU <- X*params$c*lambdaR + YRT*params$gR - YRU*(params$q + params$e*sigmaR 
                                                    + params$nuRU + params$mu) 
  dYRT <- YRU*params$e*sigmaR + YST*params$r - YRT*(params$gR + params$nuRT + params$mu)  
  list(c(dX, dY, dYSU, dYST, dYRU, dYRT))
}


###################################################
### code chunk number 3: sensitivity.rnw:142-168
###################################################
require(deSolve)
times <- seq(0, 20, by=1)  #solve for 20 years
params <- list(betaSU=0.1, 
               alpha=0.1+(0.5-0.1)/2, 
               pSU=0.5, 
               pST=0.5, 
               FS=0.7, 
               FR=0.7, 
               alpha2=0.5, 
               alpha3=0.5, 
               pi=2133, 
               c=1.7, 
               mu=1/30, 
               q=52/6, 
               gS=0.05, 
               nuSU=1/12, 
               r=0.10, 
               nuST=1/27, 
               gR=0.1, 
               e=0.5, 
               nuRU=1/27, 
               nuRT=1/27)
pop.size <- 800000*0.5*0.15
xstart <- c(X=pop.size*0.7, Y=pop.size*0.3, YSU=pop.size*0.3*0.5*0.999,
            YST=pop.size*0.3*0.5*0.85, YRU=pop.size*0.3*0.5*0.001, YRT=pop.size*0.3*0.5*0.15)
out <- as.data.frame(lsoda(xstart, times, hiv.model.sf, params))


###################################################
### code chunk number 4: sensitivity.rnw:173-178
###################################################
par(mfrow=c(2,2))
plot(out$time, out$YSU, xlab='Time', ylab='HIV+ Sensitive/Untreated')
plot(out$time, out$YST, xlab='Time', ylab='HIV+ Sensitive/Treated')
plot(out$time, out$YRU, xlab='Time', ylab='HIV+ Resistant/Untreated')
plot(out$time, out$YRT, xlab='Time', ylab='HIV+ Resistant/Treated')


###################################################
### code chunk number 5: sensitivity.rnw:185-200
###################################################
infections.prevented.baseline<-c()
for(F in seq(0,0.98,by=0.02)){
  params <- list(betaSU=0.1, alpha=0.1+(0.5-0.1)/2 , pSU=0.5, pST=0.5, FS=F, FR=F, alpha2=0.5 , alpha3=0.5 , pi=2133, c=1.7, mu=1/30, q=52/6, gS=0.05, nuSU=1/12, r=0.10, nuST=1/27, gR=0.1, e=0.5, nuRU=1/27 , nuRT=1/27)
  out <- as.data.frame(lsoda(xstart, times, hiv.model.sf, params))
  new.infections.art <- sum(diff(out$Y))

  params <- list(betaSU=0.1, alpha=0.1+(0.5-0.1)/2 , pSU=0.5, pST=0.5, FS=0, FR=0, alpha2=0.5 , alpha3=0.5 , pi=2133, c=1.7, mu=1/30, q=52/6, gS=0.05, nuSU=1/12, r=0.10, nuST=1/27, gR=0.1, e=0.5, nuRU=1/27 , nuRT=1/27)
  out <- as.data.frame(lsoda(xstart, times, hiv.model.sf, params))
  new.infections.no.art <- sum(diff(out$Y))

  infections.prevented.baseline <- c(infections.prevented.baseline, new.infections.no.art-new.infections.art)
}
plot(seq(0,0.98,by=0.02), infections.prevented.baseline, type='l', lwd=3, ylim=c(-200,20000),
            xlab='Fraction of cases treated',
            ylab='Infections prevented')


###################################################
### code chunk number 6: sensitivity.rnw:206-209
###################################################
require(lhs)            #add the lhs library
h <- 500                #choose number of points
lhs<-maximinLHS(h,18)   #simulate


###################################################
### code chunk number 7: sensitivity.rnw:214-250
###################################################
betaSU.min <- 0.0
betaSU.max <- 0.2
alpha.min <- 0.2
alpha.max <- 0.4
pSU.min <- 0.4
pSU.max <- 0.6
pST.min <- 0.4
pST.max <- 0.6
alpha2.min <- 0.4 
alpha2.max <- 0.6
alpha3.min <- 0.4
alpha3.max <- 0.6
pi.min <- 2133-250
pi.max <- 2133+250
c.min <- 1.2
c.max <- 2.2
mu.min <- 0.0333-0.01
mu.max <- 0.0333+0.01
q.min <- 8.667-1
q.max <- 8.667+1
gS.min <- 0.04
gS.max <- 0.06
nuSU.min <- 0.08333-0.01
nuSU.max <- 0.08333+0.01
r.min <- 0.05
r.max <- 0.15
nuST.min <- 0.037-0.01
nuST.max <- 0.037+0.01
gR.min <- 0.05
gR.max <- 0.15
e.min <- 0.4
e.max <- 0.6
nuRU.min <- 0.037-0.01
nuRU.max <- 0.037+0.01
nuRT.min <- 0.037-0.01
nuRT.max <- 0.037+0.01


###################################################
### code chunk number 8: sensitivity.rnw:254-273
###################################################
params.set <- cbind(
  betaSU = lhs[,1]*(betaSU.max-betaSU.min)+betaSU.min,
  alpha = lhs[,2]*(alpha.max-alpha.min)+alpha.min,
  pSU = lhs[,3]*(pSU.max-pSU.min)+pSU.min,
  pST = lhs[,4]*(pST.max-pST.min)+pST.min,
  alpha2 = lhs[,5]*(alpha2.max-alpha2.min)+alpha2.min,
  alpha3 = lhs[,6]*(alpha3.max-alpha3.min)+alpha3.min,
  pi = lhs[,7]*(pi.max-pi.min)+pi.min,
  c = lhs[,8]*(c.max-c.min)+c.min,
  mu = lhs[,9]*(mu.max-mu.min)+mu.min,
  q = lhs[,10]*(q.max-q.min)+q.min,
  gS = lhs[,11]*(gS.max-gS.min)+gS.min,
  nuSU = lhs[,12]*(nuSU.max-nuSU.min)+nuSU.min,
  r = lhs[,13]*(r.max-r.min)+r.min,
  nuST = lhs[,14]*(nuST.max-nuST.min)+nuST.min,
  gR = lhs[,15]*(gR.max-gR.min)+gR.min,
  e = lhs[,16]*(e.max-e.min)+e.min,
  nuRU = lhs[,17]*(nuRU.max-nuRU.min)+nuRU.min,
  nuRT = lhs[,18]*(nuRT.max-nuRT.min)+nuRT.min)


###################################################
### code chunk number 9: sensitivity.rnw:280-281
###################################################
l <- 19


###################################################
### code chunk number 10: sensitivity.rnw:286-287
###################################################
h2 <-250


###################################################
### code chunk number 11: sensitivity.rnw:292-314
###################################################
j <- 1  
data <- data.frame(matrix(rep(NA,l*h2*21),nrow=l*h2))
for(i in 1:h2){
  for (F in seq(0,0.9,length=l)){
    xstart <- c(X=pop.size*0.7, Y=pop.size*0.3, YSU=pop.size*0.3*0.5*0.999,
            YST=pop.size*0.3*0.5*0.85, YRU=pop.size*0.3*0.5*0.001, YRT=pop.size*0.3*0.5*0.15)

    params <- as.list(c(params.set[i,], FS=0, FR=0))
    out <- as.data.frame(lsoda(xstart, times, hiv.model.sf, params))
    new.infections.no.art <- sum(diff(out$Y))
    
    params <- as.list(c(params.set[i,], FS=F, FR=F))
    out <- as.data.frame(lsoda(xstart, times, hiv.model.sf, params))
    new.infections.art <- sum(diff(out$Y))
       
    infections.prevented <- new.infections.no.art-new.infections.art
    data[j,1:20] <- params
    data[j,21] <- infections.prevented
    j <- j+1
  }
}
names(data) <- c(names(params),'Infections.Prevented')


###################################################
### code chunk number 12: sensitivity.rnw:319-323
###################################################
plot(seq(0,0.98,by=0.02), infections.prevented.baseline, type='l', lwd=3, ylim=c(-200,20000),
            xlab='Fraction of cases treated',
            ylab='Infections prevented')
points(data$FS, data$Infections.Prevented, pch=19, cex=0.3, col='blue')


###################################################
### code chunk number 13: sensitivity.rnw:328-331
###################################################
boxplot(data$Infections.Prevented~data$FS, ylim=c(-200,20000), border='blue')
par(new=TRUE)
plot(seq(0,0.98,by=0.02), infections.prevented.baseline, type='l', lwd=3, ylim=c(-200,20000), xlab='', ylab='Infections prevented', axes=FALSE)


###################################################
### code chunk number 14: sensitivity-package
###################################################
require(sensitivity)
prcc <- pcc(data[,1:20], data[,21], nboot = 100, rank=TRUE)


###################################################
### code chunk number 15: plot-prcc
###################################################
plot(prcc)


