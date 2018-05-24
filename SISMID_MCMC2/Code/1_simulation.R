################################################################################
# R-code to simulate from a Markov and Non-Markov Stochastic epidemic models
################################################################################


# This function assumes 1 initial infective and N-1 initially susceptibles
# Per-person infection rate is beta/N

simSIR.Markov <- function(N, beta, gamma) {

  # initial number of infectives and susceptibles;
  I <- 1
  S <- N-1;

  # recording time;
  t <- 0;
  times <- c(t);

  # a vector which records the type of event (1=infection, 2=removal)
  type <- c(1); 

  while (I > 0) {  #stop the model once disease extinct
 
    # time to next event;
    t <- t + rexp(1, (beta/N)*I*S + gamma*I); # proposed time to next event
    times <- append(times, t);      # add this proposed time to the 'times' vecor
    
    if (runif(1) < beta*S/(beta*S + N*gamma)) { #add infection with propbability = pr(infection)
      # infection
      I <- I+1;
      S <- S-1;
      type <- append(type, 1); ## type=1 is an infection
    }
    else {
      #removal
      # S <- S 
      I <- I-1
      type <- append(type, 2); ## type=2 is a removal
    }
  }
  
  # record the final size , i.e. the number of initially susceptlbles who contracted the disease sometime during the epidemic.

  # record the times of events (infections/removals) as well as the type
    
  res <- list("t"=times, "type"=type); ##RESULTS: the new times vector and type i.e. what happened
  res
}

## with final size calc
simSIR.Markov.fsize <- function(N, beta, gamma) {

  # initial number of infectives and susceptibles;
  I <- 1
  S <- N-1;

  # recording time;
  t <- 0;
  times <- c(t);

  # a vector which records the type of event (1=infection, 2=removal)
  type <- c(1); 

  while (I > 0) {  #stop the model once disease extinct
 
    # time to next event;
    t <- t + rexp(1, (beta/N)*I*S + gamma*I); # proposed time to next event
    times <- append(times, t);      # add this proposed time to the 'times' vecor
    
    if (runif(1) < beta*S/(beta*S + N*gamma)) { #add infection with propbability = pr(infection)
      # infection
      I <- I+1;
      S <- S-1;
      type <- append(type, 1); ## type=1 is an infection
    }
    else {
      #removal
      # S <- S 
      I <- I-1
      type <- append(type, 2); ## type=2 is a removal
    }
  }
  
  # record the final size , i.e. the number of initially susceptlbles who contracted the disease sometime during the epidemic.
  final.size=sum(type[type==1])
	duration=max(times)
  # record the times of events (infections/removals) as well as the type
    
  res <- list("size"=final.size, "duration"=duration); ##RESULTS: the new times vector and type i.e. what happened
  res
}



## alternative algorithm for the SIR Mchain
simSIR.Markov.alternative <- function(N, beta, gamma) {

  # initial number of infectives and susceptibles;
  I <- 1
  S <- N-1;

  # recording time;
  t <- 0;
  times <- c(t);

  # a vector which records the type of event (1=infection, 2=removal)
  type <- c(1);

  while (I > 0) {

    ############################################
    # simulate times to the next possible events
    ############################################
    
    # time to next infection
    if (S > 0) { # if susceptible hosts still exist then draw random possible time to infection
      t.next.infection <- t +  rexp(1, (beta/N)*I*S)
    }else {
      t.next.infection <- Inf;
    }
    
    # ALSO generate a possibe time to next removal    
    t.next.removal <- t + rexp(1, gamma*I)


    # check which of the two events happens first
    if (t.next.infection < t.next.removal) {
      # infection occurs
      I <- I+1;
      S <- S-1;
      type <- append(type, 1); ## type 1 so infection 
      times <- append(times, t.next.infection);
      t <- t.next.infection
    } else {
      #removal occurs
      I <- I-1
      times <- append(times, t.next.removal);
      type <- append(type, 2); ## type 2 so removal
      t <- t.next.removal
    }
  }
  
  # record the final size , i.e. the number of initially susceptlbles who contracted the disease sometime during the epidemic.
  final.size=sum(type[type==1])
  duration=max(times)
  
  # record the times of events (infections/removals) as well as the type
  
  res <- list("t"=times, "type"=type, "labels"=labels, "size"=final.size, "duration"=duration);
  res
}


simSIR.Non.Markov.constant <- function(N, beta, k) {

  # initial number of infectives and susceptibles;
  I <- 1
  S <- N-1;

  # recording time;
  t <- 0;
  times <- c(t);

  # create a vector containing the removal times of all the current infectives.
  r <- k 

  # a vector which records the type of event (1=infection, 2=removal)
  type <- c(1);
  
  # a counter for labelling the individuals
  lambda <- 1;

  # a vector to store the labels
  labels <- c(1);

  while (I > 0) {

    ############################################
    # simulate times to the next possible events
    ############################################
    
    # time to next infection
    if (S > 0) {
      T  <- rexp(1, (beta/N)*I*S)
    } else {
      T <- Inf;
    }

    # time to next removal
    R <- min(r, na.rm=TRUE);
    
    # check which of the two events happens first
    if (t + T < R) {
      # infection occurs
      I <- I+1;
      S <- S-1;
      r <- append(r, t + T + k) ## k = constant infection period
      type <- append(type, 1);
      times <- append(times, t + T);

      lambda <- lambda + 1;
      labels <- append(labels, lambda)
      t <- t + T
    } else {
      #removal occurs
      I <- I-1
      type <- append(type, 2);
      index.min.r <- which(min(r, na.rm=TRUE)==r)
      r[index.min.r] <- NA
      labels <- append(labels, index.min.r)
      times <- append(times, R);
      t <- R

      # update the vector of 
    }
  }
  
  # record the final size , i.e. the number of initially susceptlbles who contracted the disease sometime during the epidemic.
  final.size=sum(type[type==1])
  duration=max(times)
  
  # record the times of events (infections/removals) as well as the type
  
  res <- list("t"=times, "type"=type, "labels"=labels, "size"=final.size, "duration"=duration);
  res
}



simSIR.Non.Markov.gamma <- function(N, beta, gamma, delta) {

  # initial number of infectives and susceptibles;
  I <- 1
  S <- N-1;

  # recording time;
  t <- 0;
  times <- c(t);

  # create a vector containing the removal times of all the current infectives.
  k <- rgamma(1, gamma, delta)
  r <- k 

  # a vector which records the type of event (1=infection, 2=removal)
  type <- c(1);

  # a counter for labelling the individuals
  lambda <- 1;
  
  # a vector to store the labels
  labels <- c(1);

  while (I > 0) {

    ############################################
    # simulate times to the next possible events
    ############################################
    
    # time to next infection
    if (S > 0) {
      T  <- rexp(1, (beta/N)*I*S)
    }
    else {
      T <- Inf;
    }

    # time to next removal
    R <- min(r, na.rm=TRUE);
    
    # check which of the two events happens first
    if (t + T < R) {
      # infection occurs
      I <- I+1;
      S <- S-1;
      k <- rgamma(1, gamma, delta)
      r <- append(r, t + T + k)

      lambda <- lambda + 1;
      labels <- append(labels, lambda)
      type <- append(type, 1);
      times <- append(times, t + T);
      t <- t + T
    }
    else {
      #removal occurs
      I <- I-1
      type <- append(type, 2);
      index.min.r <- which(min(r, na.rm=TRUE)==r)
      r[index.min.r] <- NA
      labels <- append(labels, index.min.r)
      times <- append(times, R);
      t <- R      
    }
  }
  
  # record the final size , i.e. the number of initially susceptlbles who contracted the disease sometime during the epidemic.
  final.size=sum(type[type==1])
  duration=max(times)
  
  # record the times of events (infections/removals) as well as the type
  
  res <- list("t"=times, "type"=type, "labels"=labels, "size"=final.size, "duration"=duration);
  res
}

################################################################################
# Call the functions
################################################################################
################################################################################
##' 1
set.seed(05081991)
sir1 <- simSIR.Markov(21, 0.9, 1) 
sir2 <- simSIR.Markov(21, 2, 1) 
sir3 <- simSIR.Markov(21, 4, 1) 
print(c(max(sir1$t),max(sir2$t),max(sir3$t)))
plot(sir3$t, sir3$type,  type='l', col=1, main='length of outbreaks', xlab = 'time', ylab='type: 1=infection, 2=removal')
lines(sir2$t,sir2$type,   type='l', col=2)
lines(sir1$t,sir1$type,   type='l', col=3)
## transmission rate increase therefore infection survives longer

################################################################################
##' 2 - record final size
set.seed(05081991)
sir1 <- simSIR.Markov.fsize(21, 0.9, 1) 
sir2 <- simSIR.Markov.fsize(21, 2, 1) 
sir3 <- simSIR.Markov.fsize(21, 4, 1) 
print(c(max(sir1$size),max(sir2$size),max(sir3$size))) ## final size: epidemics 3 infects everyone
print(c(max(sir1$duration),max(sir2$duration),max(sir3$duration))) ## duration: epidemics 3 longest


################################################################################
##' 3 - simulation based estimate of final size
par(mfrow=c(2,2))
set.seed(05081991)
sir_f.size <- NULL
for(beta in c(0.9, 2, 4, 8)){
for(i in 1:1000){
  sir1 <- simSIR.Markov.fsize(21, beta, 1) 
  sir_f.size <- c(sir_f.size, sir1$size)
}
hist(sir_f.size, main=paste0('R0 = ', beta))
}

## Now with constant gamma (not assumed exponential dist.)
set.seed(05081991)
sir_f.size <- NULL
for(beta in c(0.9, 2, 4, 8)){
  for(i in 1:1000){
    sir1 <- simSIR.Non.Markov.constant(21, beta, 1) 
    sir_f.size <- c(sir_f.size, sir1$size)
  }
  hist(sir_f.size, col=2, main=paste0('R0 = ', beta/1, " gamma = 1"))
}

## Now with gamma dist infectious period
set.seed(05081991)
sir_f.size <- NULL
for(beta in c(0.9, 2, 4, 8)){
  for(i in 1:1000){
    sir1 <- simSIR.Non.Markov.gamma(21, beta, 1, 0.5) 
    sir_f.size <- c(sir_f.size, sir1$size)
  }
  hist(sir_f.size, col=3, main=paste0('Beta = ', beta/1, " gamma ~ (1,0.5)"))
}

set.seed(05081991)
sir_f.size1 <- NULL
sir_f.size2 <- NULL
sir_f.size3 <- NULL
sir_f.size4 <- NULL
beta = 2
  for(i in 1:1000){
    sir1 <- simSIR.Markov.fsize(21, beta, 1) 
    sir2 <- simSIR.Markov.alternative(21, beta, 1) 
    sir3 <- simSIR.Non.Markov.constant(21, beta, 1) 
    sir4 <- simSIR.Non.Markov.gamma(21, beta, 1, 0.5) 
      sir_f.size1 <- c(sir_f.size1, sir1$size)
      sir_f.size2 <- c(sir_f.size2, sir2$size)
      sir_f.size3 <- c(sir_f.size3, sir3$size)
      sir_f.size4 <- c(sir_f.size4, sir4$size)
  }
hist(sir_f.size1, main=paste0(" gamma ~ Exp(1)"), xlab="Size of epidemic")
hist(sir_f.size2, col=2, main=paste0(" gamma ~ Exp(1) ALT"), xlab="Size of epidemic")
hist(sir_f.size3, col=3, main=paste0(" gamma = 1"), xlab="Size of epidemic")
hist(sir_f.size4, col=4, main=paste0(" gamma ~ (1,0.5)"), xlab="Size of epidemic")


################################################################################
##' 4 dist. of duration
set.seed(05081991)
sir_dur1 <- NULL
sir_dur2 <- NULL
sir_dur3 <- NULL
sir_dur4 <- NULL
beta = 2
  for(i in 1:1000){
    sir1 <- simSIR.Markov.fsize(21, beta, 1) 
    sir2 <- simSIR.Markov.alternative(21, beta, 1) 
    sir3 <- simSIR.Non.Markov.constant(21, beta, 1) 
    sir4 <- simSIR.Non.Markov.gamma(21, beta, 1, 0.5) 
      sir_dur1 <- c(sir_dur1, sir1$duration)
      sir_dur2 <- c(sir_dur2, sir2$duration)
      sir_dur3 <- c(sir_dur3, sir3$duration)
      sir_dur4 <- c(sir_dur4, sir4$duration)
  }
hist(sir_dur1, main=paste0(" gamma ~ Exp(1)"), xlab="Duration of outbreak")
hist(sir_dur2, col=2, main=paste0(" gamma ~ Exp(1) ALT"), xlab="Duration of outbreak")
hist(sir_dur3, col=3, main=paste0(" gamma = 1"), xlab="Duration of outbreak")
hist(sir_dur4, col=4, main=paste0(" gamma ~ (1,0.5)"), xlab="Duration of outbreak")


################################################################################
##' 5 - gamma with WEIBULL dist
##' 
simSIR.Non.Markov.weibull <- function(N, beta, gamma, delta) {

# initial number of infectives and susceptibles;
I <- 1
S <- N-1;

# recording time;
t <- 0;
times <- c(t);

# create a vector containing the removal times of all the current infectives.
k <- rweibull(1, gamma, delta)
r <- k 

# a vector which records the type of event (1=infection, 2=removal)
type <- c(1);

# a counter for labelling the individuals
lambda <- 1;

# a vector to store the labels
labels <- c(1);

while (I > 0) {
  
  ############################################
  # simulate times to the next possible events
  ############################################
  
  # time to next infection
  if (S > 0) {
    T  <- rexp(1, (beta/N)*I*S)
  }
  else {
    T <- Inf;
  }
  
  # time to next removal
  R <- min(r, na.rm=TRUE);
  
  # check which of the two events happens first
  if (t + T < R) {
    # infection occurs
    I <- I+1;
    S <- S-1;
    k <- rweibull(1, gamma, delta)
    r <- append(r, t + T + k)
    
    lambda <- lambda + 1;
    labels <- append(labels, lambda)
    type <- append(type, 1);
    times <- append(times, t + T);
    t <- t + T
  }
  else {
    #removal occurs
    I <- I-1
    type <- append(type, 2);
    index.min.r <- which(min(r, na.rm=TRUE)==r)
    r[index.min.r] <- NA
    labels <- append(labels, index.min.r)
    times <- append(times, R);
    t <- R      
  }
}

# record the final size , i.e. the number of initially susceptlbles who contracted the disease sometime during the epidemic.
final.size=sum(type[type==1])
duration=max(times)

# record the times of events (infections/removals) as well as the type

res <- list("t"=times, "type"=type, "labels"=labels, "size"=final.size, "duration"=duration);
res
}

set.seed(05081991)
sir_f.size <- NULL
for(beta in c(0.9, 2, 4, 8)){
  for(i in 1:1000){
    sir1 <- simSIR.Non.Markov.weibull(21, beta, 1, 0.5) 
    sir_f.size <- c(sir_f.size, sir1$size)
  }
  hist(sir_f.size, col=5, main=paste0('beta = ', beta, " gamma ~ WEIBULL(1,0.5)"))
}

set.seed(05081991)
sir_f.size1 <- NULL
sir_f.size2 <- NULL
sir_f.size3 <- NULL
sir_f.size4 <- NULL
beta = 2
for(i in 1:1000){
  sir1 <- simSIR.Markov.fsize(21, beta, 1) 
  sir2 <- simSIR.Non.Markov.weibull(21, beta, 1, 0.5) 
  sir3 <- simSIR.Non.Markov.constant(21, beta, 1) 
  sir4 <- simSIR.Non.Markov.gamma(21, beta, 1, 0.5) 
  sir_f.size1 <- c(sir_f.size1, sir1$size)
  sir_f.size2 <- c(sir_f.size2, sir2$size)
  sir_f.size3 <- c(sir_f.size3, sir3$size)
  sir_f.size4 <- c(sir_f.size4, sir4$size)
}
hist(sir_f.size1, main=paste0(" gamma ~ Exp(1)"), xlab="Size of epidemic")
hist(sir_f.size2, col=2, main=paste0(" gamma ~ Weibull(1,0.5)"), xlab="Size of epidemic")
hist(sir_f.size3, col=3, main=paste0(" gamma = 1"), xlab="Size of epidemic")
hist(sir_f.size4, col=4, main=paste0(" gamma ~ (1,0.5)"), xlab="Size of epidemic")


################################################################################
##' 6 - SEIR model
##' 

simSEIR.Non.Markov.weibull <- function(N, beta, alpha, gamma, delta) {

# initial number of infectives and susceptibles;
E <- 0
I <- 1
S <- N-1;

# recording time;
t <- 0;
times <- c(t);

# create a vector containing the removal times of all the current infectives.
k <- rweibull(1, gamma, delta)
r <- k 
e <- 0

# a vector which records the type of event (1=infection, 2=removal)
type <- c(1);

# a counter for labelling the individuals
lambda <- 1;

# a vector to store the labels
labels <- c(1);

while (I > 0) {
  
  ############################################
  # simulate times to the next possible events
  ############################################
  
  # time to next infection
  if (S > 0) {
    T  <- rexp(1, (beta/N)*I*S)
  }
  else {
    T <- Inf;
  }
  
  # time to next removal
  E <- min(e, na.rm=T)
  R <- min(r, na.rm=TRUE);
  
  # check which of the two events happens first
  if (t + T + alpha < R) {
    # infection occurs
    I <- I+1;
    S <- S-1;
    k <- rweibull(1, gamma, delta)
    r <- append(r, t + T + alpha + k)
    
    lambda <- lambda + 1;
    labels <- append(labels, lambda)
    type <- append(type, 1);
    times <- append(times, t + T + alpha);
    t <- t + T + alpha
  }
  else {
    #removal occurs
    I <- I-1
    type <- append(type, 2);
    index.min.r <- which(min(r, na.rm=TRUE)==r)
    r[index.min.r] <- NA
    labels <- append(labels, index.min.r)
    times <- append(times, R);
    t <- R      
  }
  ## Include E state change.....
}

# record the final size , i.e. the number of initially susceptlbles who contracted the disease sometime during the epidemic.
final.size=sum(type[type==1])
duration=max(times)

# record the times of events (infections/removals) as well as the type

res <- list("t"=times, "type"=type, "labels"=labels, "size"=final.size, "duration"=duration);
res
}

par(mfrow=c(2,2))
set.seed(05081991)
sir_f.size <- NULL
sir_dur <- NULL
for(beta in c(0.9, 4)){
  for(i in 1:1000){
    sir1 <- simSEIR.Non.Markov.weibull(21, beta, alpha=1, gamma=1, delta=0.5) 
    sir_f.size <- c(sir_f.size, sir1$size)
    sir_dur <- c(sir_dur, sir1$duration)
  }
  hist(sir_f.size, col=5, main=paste0('FINAL SIZE: beta = ', beta), xlab='Final size')
  hist(sir_dur, col=5, main=paste0('DURATION: beta = ', beta), xlab='Duration')
}

