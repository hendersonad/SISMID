### R code from vignette source 'social-distancing.rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: packages
###################################################
require(deSolve)

#define a color palette
pulp <- c(rgb(101,48,47, maxColorValue=255),    #brown
          rgb(210,202,203, maxColorValue=255),  #gray
          rgb(211,141,101, maxColorValue=255),  #peach
          rgb(223,45,39, maxColorValue=255),    #red
          rgb(250,208,10, maxColorValue=255),   #yellow
          rgb(16,16,18, maxColorValue=255),     #black
          rgb(76,106,147, maxColorValue=255)    #blue
          )


###################################################
### code chunk number 2: sir-closed
###################################################
sir.model.closed <- function (t, x, params) {    #here we begin a function with three arguments
  S <- x[1]                               #create local variable S, the first element of x
  I <- x[2]                               #create local variable I
  R <- x[3]                               #create local variable R
  with(                                   #we can simplify code using "with"
       as.list(params),                   #this argument to "with" lets us use the variable names
       {                                  #the system of rate equations
         dS <- -beta*S*I
         dI <- beta*S*I-gamma*I
         dR <- gamma*I
         dx <- c(dS,dI,dR)                #combine results into a single vector dx
         list(dx)                         #return result as a list
       }
       )
}


###################################################
### code chunk number 3: parameters
###################################################
R0 <- 1.8
N <-  1  				              #population size
gamma <- 7/2.6  			        #recovery rate (in years)
beta <- R0*(gamma)/N 	        #transmission rate
phi <- 0.3333
D <- 12
xstart <- c(S=0, I=0, R=0)	  #initial conditions, must sum to one
xstart[2] <- 1/58000000
xstart[1] <- 1-xstart[2]
xstart[3] <- 1-sum(xstart)
Tmax <- 52                    #integrate for 200 years after transients
params <- c(beta=beta, gamma=gamma)   #parameter vector
tau <- 0.1                    #size of time step
times <- seq(0, Tmax, by=tau) #function seq returns a sequence


###################################################
### code chunk number 4: social-distancing.rnw:125-126
###################################################
out <- as.data.frame(ode(xstart,times,sir.model.closed,params, method='ode45', rtol=1e-7))


###################################################
### code chunk number 5: plot-1
###################################################
plot(out$time,
     out$I,
     type='l',
     lwd=2,
     xlab='Time',
     ylab='Infected',
     col=pulp[7])  #plot the I variable against time


###################################################
### code chunk number 6: sd-initiated
###################################################
Tvec <- c(3, 5, 6, 7)


###################################################
### code chunk number 7: color (eval = FALSE)
###################################################
## require(EBImage)   # available on Bioconductor
## require(rPlotter)  # https://github.com/woobe/rPlotter
## 
## ## Define a list of images from Tarantino's movies
## lst_tar <- list(
##   reservoir_dogs = "http://filmhash.files.wordpress.com/2011/06/reservoir-dogs-051.jpg",
##   pulp_fiction = "http://www.scoutlondon.com/blog/wp-content/uploads/2012/05/Pulp-Fiction.jpg",
##   kill_bill = "http://www.moviegoods.com/Assets/product_images/1010/477803.1010.A.jpg",
##   django = "http://www.comingsoon.net/nextraimages/djangounchainednewposter.jpg"
##   )
## 
## 
## ## Create palette for each image and save them all into one PNG
## png("example_tarantino.png", width = 1000, height = 1000, res = 150)
## par(mfrow = c(4,4))
## 
## for (n_tar in 1:length(lst_tar)) {
##   tmp_url <- unlist(lst_tar[n_tar])
##   if (n_tar %% 2 != 0) display(readImage(tmp_url), method = "raster")
##   set.seed(1234)
##   pie(rep(1, 3), col = extract_colours(tmp_url, 3))
##   pie(rep(1, 5), col = extract_colours(tmp_url, 5))
##   if (n_tar %% 2 == 0) display(readImage(tmp_url), method = "raster")
## }


###################################################
### code chunk number 8: solve
###################################################

# plot baseline again
plot(out$time,
     out$I,
     type='l',
     lwd=2,
     xlab='Time',
     ylab='Infected',
     col=pulp[7])  #plot the I variable against time

# create a list to store information for the legend
legend.info <- data.frame(T=numeric(0),col=character(0),stringsAsFactors=FALSE)      

for(i in 1:length(Tvec)){
  T <- Tvec[i]

  out1 <- ode(xstart,                               # start at initial condition
             seq(0,T,tau),                          # solve from 0 to T
             sir.model.closed,params,
             method='ode45',
             rtol=1e-7)  
  
  out2 <- ode(tail(out1,1)[2:4],                    # start at end of last solution
             seq(T,T+D,tau),                        # solve from T to T+D
             sir.model.closed,
             c(beta=beta*(1-phi), gamma=gamma),     # change beta
             method='ode45',
             rtol=1e-7)
  
  out3 <- ode(tail(out2,1)[2:4],                    # start at end of last solution
             seq(T+D, Tmax, tau),                   # solve from T+D to Tmax
             sir.model.closed,
             params,                                # reset parameters
             method='ode45',
             rtol=1e-7) 
  
  data <- as.data.frame(rbind(out1, out2, out3))
  
  lines(data$time,
        data$I,
        col=pulp[i],
        lwd=2)
  
  legend.info[i,] <- c(T=T, col=as.character(pulp[i]))
}

legend('topright',
       legend=paste('T=',legend.info$T, sep=''),
       lty=1,
       lwd=3,
       col=legend.info$col,
       bty='n')


