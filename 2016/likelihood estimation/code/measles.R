#' ---
#' title: "Case study: Measles in large and small towns"
#' author: "Aaron A. King"
#' output:
#'   html_document:
#'     toc: yes
#'     toc_depth: 3
#' bibliography: ../sbied.bib
#' csl: ../ecology.csl
#' 
#' ---
#' 
#' \newcommand\prob[1]{\mathbb{P}\left[{#1}\right]}
#' \newcommand\expect[1]{\mathbb{E}\left[{#1}\right]}
#' \newcommand\var[1]{\mathrm{Var}\left[{#1}\right]}
#' \newcommand\dist[2]{\mathrm{#1}\left(#2\right)}
#' \newcommand\dlta[1]{{\Delta}{#1}}
#' \newcommand\scinot[2]{$#1 \times 10^{#2}$\xspace}
#' \newcommand{\mortality}{m}
#' \newcommand{\birth}{b}
#' \newcommand{\loglik}{\ell}
#' \newcommand{\immigration}{\iota}
#' \newcommand{\amplitude}{a}
#' \newcommand{\cohort}{c}
#' \newcommand{\R}{\textsf{R}}
#' \newcommand{\Rzero}{{R_0}}
#' 
#' [Licensed under the Creative Commons Attribution-NonCommercial license](http://creativecommons.org/licenses/by-nc/4.0/).
#' Please share and remix noncommercially, mentioning its origin.  
#' ![CC-BY_NC](../graphics/cc-by-nc.png)
#' 
#' Produced in **R** version `r getRversion()` using **pomp** version `r packageVersion("pomp")`.
#' 
#' 
#' ----------------------------------
#' 
#' ## Objectives
#' 
#' 1. To display a published case study using plug-and-play methods with non-trivial model complexities
#' 1. To show how extra-demographic stochasticity can be modeled
#' 1. To demonstrate the use of covariates in **pomp**
#' 1. To demonstrate the use of profile likelihood in scientific inference
#' 1. To discuss the interpretation of parameter estimates
#' 1. To emphasize the potential need for extra sources of stochasticity in modeling
#' 
#' ## Measles revisited
#' 
#' ### Motivation: challenges in inference from disease dynamics
#' 
#' - Understanding, forecasting, managing epidemiological systems increasingly depends on models.
#' - Dynamic models can be used to test causal hypotheses.
#' - Real epidemiological systems:
#'     - are nonlinear
#'     - are stochastic
#'     - are nonstationary
#'     - evolve in continuous time
#'     - have hidden variables
#'     - can be measured only with (large) error
#' - Dynamics of infectious disease outbreaks illustrate this well.
#' 
#' - Measles is the paradigm for a nonlinear ecological system that can be well described by low-dimensional nonlinear dynamics.
#' - A tradition of careful modeling studies have proposed and found evidence for a number of specific mechanisms, including
#'     - a high value of $R_0$ (c. 15--20)
#'     - seasonality in transmission rates associated with school terms
#'     - a birth cohort effect
#'     - response to changing birth rates
#'     - under-reporting
#'     - fadeouts and reintroductions that scale with city size
#'     - spatial traveling waves
#' - Much of this evidence has been amassed from fitting models to data, using a variety of methods.
#' - See @Rohani2010 for a review of some of the high points.
#' 
#' 
#' ### Outline
#' 
#' - We revisit a classic measles data set, weekly case reports in 954 urban centers in England and Wales during the pre-vaccine era (1950--1963).
#' - We examine questions regarding:
#'     - measles extinction and recolonization
#'     - transmission rates
#'     - seasonality
#'     - resupply of susceptibles
#' - We use a model that 
#'     1. expresses our current understanding of measles dynamics
#'     1. includes a long list of mechanisms that have been proposed and demonstrated in the literature
#'     1. cannot be fit by existing likelihood-based methods
#' - We examine data from large and small towns using the same model, something no existing methods have been able to do.
#' - We ask: does our perspective on this disease change when we expect the models to explain the data in detail?
#' - What bigger lessons can we learn regarding inference for dynamical systems?
#' 
#' ### He, Ionides, & King, *J. R. Soc. Interface* (2010)
#' 
#' #### Data sets
#' 
#' - Twenty towns, including
#'     - 10 largest
#'     - 10 smaller, chosen at random
#' - Population sizes: 2k--3.4M
#' - Weekly case reports, 1950--1963
#' - Annual birth records and population sizes, 1944--1963
#' 
#' 
#' 
#' 
#' #### Continuous-time Markov process model
#' 
#' Diagram of the model:
#' 
#' 
#' ![](./model_diagram.png)
#' 
#' - $B(t) = \text{birth rate, from data}$
#' - $N(t) = \text{population size, from data}$
#' - Overdispersed binomial measurement model: $\mathrm{cases}_t\,\vert\,\dlta{N}_{IR}=z_t \sim \dist{Normal}{\rho\,z_t,\rho\,(1-\rho)\,z_t+(\psi\,\rho\,z_t)^2}$
#' 
#' - Entry into susceptible class:
#' $$\mu_{BS}(t) = (1-c)\,B(t-\tau)+c\,\delta(t-t_0)\,\int_{t-1}^{t}\,B(t-\tau-s)\,ds$$
#' 
#' - Force of infection
#' $$\mu_{SE}(t) = \tfrac{\beta(t)}{N(t)}\,(I+\iota)\,\zeta(t)$$
#' 
#' - $c = \text{cohort effect}$  
#' - $\tau = \text{school-entry delay}$  
#' - $\beta(t) = \text{school-term transmission} = \begin{cases}\beta_0\,(1+a) &\text{during term}\\\beta_0\,(1-a) &\text{during vacation}\end{cases}$  
#' - $\iota = \text{imported infections}$
#' - $\zeta(t) = \text{Gamma white noise with intensity}\,\sigma_{SE}$ [@He2010,@bhadra11]
#' 
#' 
#' ## Model implementation
#' 
#' We'll load the packages we'll need, and set the random seed, to allow reproducibility.
#' Note that we'll be making heavy use of the data-munging methods in packages **plyr** and **reshape2**, a [tutorial introduction to which is provided here](https://kingaa.github.io/R_Tutorial/munging.html).
#' Also, we'll be using **ggplot2** for plotting: see [this brief tutorial](https://kingaa.github.io/R_Tutorial/viz.html#a-more-systematic-approach-the-grammar-of-graphics).
#' Finally, we'll use the convenient **magrittr** syntax, which is explained [here](https://kingaa.github.io/R_Tutorial/munging.html#the-magrittr-syntax).
#' 
## ----prelims,cache=FALSE-------------------------------------------------
library(pomp)
library(plyr)
library(reshape2)
library(magrittr)
library(ggplot2)
theme_set(theme_bw())
stopifnot(packageVersion("pomp")>="1.12")
set.seed(594709947L)

#' 
#' ### Data and covariates
#' 
#' Now we'll load the data and covariates.
#' The data are measles reports from 20 cities in England and Wales.
#' We also have information on the population sizes and birth-rates in these cities;
#' we'll treat these variables as covariates.
#' 
## ----load-data-----------------------------------------------------------
daturl <- "https://kingaa.github.io/pomp/vignettes/twentycities.rda"
datfile <- file.path(tempdir(),"twentycities.rda")
download.file(daturl,destfile=datfile,mode="wb")
load(datfile)

#' 
#' We select the data for London and pre-process the measles and demography data.
#' 
## ----plot-data-----------------------------------------------------------
measles %>% 
  mutate(year=as.integer(format(date,"%Y"))) %>%
  subset(town=="London" & year>=1950 & year<1964) %>%
  mutate(time=(julian(date,origin=as.Date("1950-01-01")))/365.25+1950) %>%
  subset(time>1950 & time<1964, select=c(time,cases)) -> dat
demog %>% subset(town=="London",select=-town) -> demogLondon

#' 
#' We plot the data and covariates.
#' 
## ----data-plot-----------------------------------------------------------
dat %>% ggplot(aes(x=time,y=cases))+geom_line()
demogLondon %>% melt(id="year") %>%
  ggplot(aes(x=year,y=value))+geom_point()+
  facet_wrap(~variable,ncol=1,scales="free_y")

#' 
#' Now, we smooth the covariates.
#' Note that we delay the entry of newborns into the susceptible pool.
#' 
## ----prep-covariates-----------------------------------------------------
demogLondon %>% 
  summarize(
    time=seq(from=min(year),to=max(year),by=1/12),
    pop=predict(smooth.spline(x=year,y=pop),x=time)$y,
    birthrate=predict(smooth.spline(x=year+0.5,y=births),x=time-4)$y
  ) -> covar

## ----covarplot-----------------------------------------------------------
plot(pop~time,data=covar,type='l')
points(pop~year,data=demogLondon)
plot(birthrate~time,data=covar,type='l')
points(births~year,data=demogLondon)
plot(birthrate~I(time-4),data=covar,type='l')
points(births~I(year+0.5),data=demogLondon)

#' 
#' ### The partially observed Markov process model
#' 
#' #### The (unobserved) process model
#' 
#' We propose a variant of the SEIR model as an explanation for these data.
#' This is a compartmental model that, diagrammatically, looks as follows.
#' 
## ----seir-diagram,echo=FALSE,cache=FALSE,eval=FALSE----------------------
## library(DiagrammeR)
## DiagrammeR("digraph SEIR {
##   graph [rankdir=TD, overlap=false, fontsize = 10]
##   node[shape=egg, label='B'] b;
##   subgraph {
##     rank=same;
##     node[shape=oval, label='S'] S;
##     node[shape=oval, label='E'] E;
##     node[shape=oval, label='I'] I;
##     node[shape=oval, label='R'] R;
##     S->E E->I I->R
##   }
##   node[shape=diamond, label='dead'] d;
##   b->S
##   {S E I R}->d
##    }",type="grViz",engine="dot",height=300,width=800)

#' 
#' ![model diagram](./model_diagram.png)
#' 
#' $B = \text{births}$  
#' $S = \text{susceptibles}$  
#' $E = \text{exposed, incubating}$  
#' $I = \text{infectious}$  
#' $R = \text{recovered}$  
#' 
#' We require a simulator for this model.
#' The following code implements a simulator.
#' 
#' Notable complexities include:
#' 
#' 1. Incorporation of the known birthrate.
#' 1. The birth-cohort effect: a specified fraction (`cohort`) of the cohort enter the susceptible pool aall at once.
#' 1. Seasonality in the transmission rate: during school terms, the transmission rate is higher than it is during holidays.
#' 1. Extra-demographic stochasticity in the form of a Gamma white-noise term acting multiplicatively on the force of infection.
#' 1. Demographic stochasticity implmented using Euler-multinomial distributions.
#' 
## ----rprocess------------------------------------------------------------
rproc <- Csnippet("
  double beta, br, seas, foi, dw, births;
  double rate[6], trans[6];
  
  // cohort effect
  if (fabs(t-floor(t)-251.0/365.0) < 0.5*dt) 
    br = cohort*birthrate/dt + (1-cohort)*birthrate;
  else 
  	br = (1.0-cohort)*birthrate;

  // term-time seasonality
  t = (t-floor(t))*365.25;
  if ((t>=7&&t<=100) || (t>=115&&t<=199) || (t>=252&&t<=300) || (t>=308&&t<=356))
      seas = 1.0+amplitude*0.2411/0.7589;
    else
      seas = 1.0-amplitude;

  // transmission rate
  beta = R0*(gamma+mu)*seas;
  // expected force of infection
  foi = beta*pow(I+iota,alpha)/pop;
  // white noise (extrademographic stochasticity)
  dw = rgammawn(sigmaSE,dt);

  rate[0] = foi*dw/dt;  // stochastic force of infection
  rate[1] = mu;			    // natural S death
  rate[2] = sigma;		  // rate of ending of latent stage
  rate[3] = mu;			    // natural E death
  rate[4] = gamma;		  // recovery
  rate[5] = mu;			    // natural I death

  // Poisson births
  births = rpois(br*dt);
  
  // transitions between classes
  reulermultinom(2,S,&rate[0],dt,&trans[0]);
  reulermultinom(2,E,&rate[2],dt,&trans[2]);
  reulermultinom(2,I,&rate[4],dt,&trans[4]);

  S += births   - trans[0] - trans[1];
  E += trans[0] - trans[2] - trans[3];
  I += trans[2] - trans[4] - trans[5];
  R = pop - S - E - I;
  W += (dw - dt)/sigmaSE;  // standardized i.i.d. white noise
  C += trans[4];           // true incidence
")

#' 
#' In the above, $C$ represents the true incidence, i.e., the number of new infections occurring over an interval.
#' Since recognized measles infections are quarantined, we argue that most infection occurs before case recognition so that true incidence is a measure of the number of individuals progressing from the I to the R compartment in a given interval.
#' 
#' We complete the process model definition by specifying the distribution of initial unobserved states.
#' The following codes assume that the fraction of the population in each of the four compartments is known.
#' 
## ----initializer---------------------------------------------------------
initlz <- Csnippet("
  double m = pop/(S_0+E_0+I_0+R_0);
  S = nearbyint(m*S_0);
  E = nearbyint(m*E_0);
  I = nearbyint(m*I_0);
  R = nearbyint(m*R_0);
  W = 0;
  C = 0;
")

#' 
#' 
#' #### The measurement model
#' 
#' We'll model both under-reporting and measurement error.
#' We want $\mathbb{E}[\text{cases}|C] = \rho\,C$, where $C$ is the true incidence and $0<\rho<1$ is the reporting efficiency.
#' We'll also assume that $\mathrm{Var}[\text{cases}|C] = \rho\,(1-\rho)\,C + (\psi\,\rho\,C)^2$, where $\psi$ quantifies overdispersion.
#' Note that when $\psi=0$, the variance-mean relation is that of the binomial distribution.
#' To be specific, we'll choose
#' $\text{cases|C} \sim f(\cdot|\rho,\psi,C)$,
#' where $$f(c|\rho,\psi,C) = \Phi(c+\tfrac{1}{2},\rho\,C,\rho\,(1-\rho)\,C+(\psi\,\rho\,C)^2)-\Phi(c-\tfrac{1}{2},\rho\,C,\rho\,(1-\rho)\,C+(\psi\,\rho\,C)^2),$$
#' where $\Phi(x,\mu,\sigma^2)$ is the c.d.f. of the normal distribution with mean $\mu$ and variance $\sigma^2$.
#' 
#' The following computes $\mathbb{P}[\text{cases}|C]$.
#' 
## ----dmeasure------------------------------------------------------------
dmeas <- Csnippet("
  double m = rho*C;
  double v = m*(1.0-rho+psi*psi*m);
  double tol = 1.0e-18;
  if (cases > 0.0) {
    lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)-pnorm(cases-0.5,m,sqrt(v)+tol,1,0)+tol;
  } else {
    lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)+tol;
  }
")

#' 
#' The following codes simulate $\text{cases} | C$.
#' 
## ----rmeasure------------------------------------------------------------
rmeas <- Csnippet("
  double m = rho*C;
  double v = m*(1.0-rho+psi*psi*m);
  double tol = 1.0e-18;
  cases = rnorm(m,sqrt(v)+tol);
  if (cases > 0.0) {
    cases = nearbyint(cases);
  } else {
    cases = 0.0;
  }
")

#' 
#' #### Constructing the `pomp` object
#' 
#' We put all the model components together with the data in a call to `pomp`:
#' 
## ----pomp-construction---------------------------------------------------
dat %>% 
  pomp(t0=with(dat,2*time[1]-time[2]),
       time="time",
       rprocess=euler.sim(rproc,delta.t=1/365.25),
       initializer=initlz,
       dmeasure=dmeas,
       rmeasure=rmeas,
       covar=covar,
       tcovar="time",
       zeronames=c("C","W"),
       statenames=c("S","E","I","R","C","W"),
       paramnames=c("R0","mu","sigma","gamma","alpha","iota",
                    "rho","sigmaSE","psi","cohort","amplitude",
                    "S_0","E_0","I_0","R_0")
  ) -> m1

#' 
#' The following codes plot the data and covariates together.
#' 
## ----plot-pomp-----------------------------------------------------------
m1 %>% as.data.frame() %>% 
  melt(id="time") %>%
  ggplot(aes(x=time,y=value))+
  geom_line()+
  facet_grid(variable~.,scales="free_y")

#' 
#' @He2010 estimated the parameters of this model.
#' The full set is included in the *R* code accompanying this document, where they are read into a data frame called `mles`.
#' 
## ----mles,include=FALSE--------------------------------------------------
read.csv(text="
town,loglik,loglik.sd,mu,delay,sigma,gamma,rho,R0,amplitude,alpha,iota,cohort,psi,S_0,E_0,I_0,R_0,sigmaSE
Bedwellty,-1125.1,0.14,0.02,4,57.9,146,0.311,24.7,0.16,0.937,0.0396,0.351,0.951,0.0396,2.64e-05,2.45e-05,0.96,0.0611
Birmingham,-3239.3,1.55,0.02,4,45.6,32.9,0.544,43.4,0.428,1.01,0.343,0.331,0.178,0.0264,8.96e-05,0.000335,0.973,0.0611
Bradford,-2586.6,0.68,0.02,4,45.6,129,0.599,32.1,0.236,0.991,0.244,0.297,0.19,0.0365,7.41e-06,4.59e-06,0.964,0.0451
Bristol,-2681.6,0.5,0.02,4,64.3,82.6,0.626,26.8,0.203,1.01,0.441,0.344,0.201,0.0358,9.62e-06,5.37e-06,0.964,0.0392
Cardiff,-2364.9,0.73,0.02,4,39,143,0.602,34.4,0.223,0.996,0.141,0.267,0.27,0.0317,1.01e-05,9.21e-06,0.968,0.0539
Consett,-1362.9,0.73,0.02,4,42.6,172,0.65,35.9,0.2,1.01,0.0731,0.31,0.406,0.0322,1.83e-05,1.97e-05,0.968,0.0712
Dalton.in.Furness,-726.1,0.3,0.02,4,73.6,257,0.455,28.3,0.203,0.989,0.0386,0.421,0.818,0.0387,2.23e-05,2.36e-05,0.961,0.0779
Halesworth,-318.6,0.51,0.02,4,49.6,210,0.754,33.1,0.381,0.948,0.00912,0.547,0.641,0.0526,1.99e-05,2.82e-05,0.947,0.0748
Hastings,-1583.7,0.21,0.02,4,56.3,74.1,0.695,34.2,0.299,1,0.186,0.329,0.396,0.0233,5.61e-06,3.4e-06,0.977,0.0955
Hull,-2729.4,0.39,0.02,4,42.1,73.9,0.582,38.9,0.221,0.968,0.142,0.275,0.256,0.0371,1.2e-05,1.13e-05,0.963,0.0636
Leeds,-2918.6,0.23,0.02,4,40.7,35.1,0.666,47.8,0.267,1,1.25,0.592,0.167,0.0262,6.04e-05,3e-05,0.974,0.0778
Lees,-548.1,1.1,0.02,4,45.6,244,0.612,29.7,0.153,0.968,0.0311,0.648,0.681,0.0477,2.66e-05,2.08e-05,0.952,0.0802
Liverpool,-3403.1,0.34,0.02,4,49.4,39.3,0.494,48.1,0.305,0.978,0.263,0.191,0.136,0.0286,0.000184,0.00124,0.97,0.0533
London,-3804.9,0.16,0.02,4,28.9,30.4,0.488,56.8,0.554,0.976,2.9,0.557,0.116,0.0297,5.17e-05,5.14e-05,0.97,0.0878
Manchester,-3250.9,0.66,0.02,4,34.4,56.8,0.55,32.9,0.29,0.965,0.59,0.362,0.161,0.0489,2.41e-05,3.38e-05,0.951,0.0551
Mold,-296.5,0.25,0.02,4,67.4,301,0.131,21.4,0.271,1.04,0.0145,0.436,2.87,0.064,2.61e-05,2.27e-05,0.936,0.0544
Northwich,-1195.1,2.25,0.02,4,45.6,147,0.795,30.1,0.423,0.948,0.0602,0.236,0.402,0.0213,1.32e-05,1.58e-05,0.979,0.0857
Nottingham,-2703.5,0.53,0.02,4,70.2,115,0.609,22.6,0.157,0.982,0.17,0.34,0.258,0.05,1.36e-05,1.41e-05,0.95,0.038
Oswestry,-696.1,0.49,0.02,4,37.3,168,0.631,52.9,0.339,1.04,0.0298,0.263,0.476,0.0218,1.56e-05,1.61e-05,0.978,0.0699
Sheffield,-2810.7,0.21,0.02,4,54.3,62.2,0.649,33.1,0.313,1.02,0.853,0.225,0.175,0.0291,6.04e-05,8.86e-05,0.971,0.0428
",stringsAsFactors=FALSE) -> mles

## ----mle-----------------------------------------------------------------
mles %>% subset(town=="London") -> mle
paramnames <- c("R0","mu","sigma","gamma","alpha","iota",
                "rho","sigmaSE","psi","cohort","amplitude",
                "S_0","E_0","I_0","R_0")
mle %>% extract(paramnames) %>% unlist() -> theta
mle %>% subset(select=-c(S_0,E_0,I_0,R_0)) %>%
  knitr::kable(row.names=FALSE)

#' 
#' Verify that we get the same likelihood as @He2010.
#' 
## ----pfilter1------------------------------------------------------------
library(foreach)
library(doParallel)

registerDoParallel()

set.seed(998468235L,kind="L'Ecuyer")

foreach(i=1:4,
        .packages="pomp",
        .options.multicore=list(set.seed=TRUE)
) %dopar% {
  pfilter(m1,Np=10000,params=theta)
} -> pfs
logmeanexp(sapply(pfs,logLik),se=TRUE)

#' 
#' Simulations at the MLE:
#' 
## ----sims1,fig.height=8--------------------------------------------------
m1 %>% 
  simulate(params=theta,nsim=9,as.data.frame=TRUE,include.data=TRUE) %>%
  ggplot(aes(x=time,y=cases,group=sim,color=(sim=="data")))+
  guides(color=FALSE)+
  geom_line()+facet_wrap(~sim,ncol=2)

#' 
#' #### Parameter transformations
#' 
#' The parameters are constrained to be positive, and some of them are constrained to lie between $0$ and $1$.
#' We can turn the likelihood maximization problem into an unconstrained maximization problem by transforming the parameters.
#' The following Csnippets implement such a transformation and its inverse.
#' 
## ----transforms----------------------------------------------------------
toEst <- Csnippet("
  Tmu = log(mu);
  Tsigma = log(sigma);
  Tgamma = log(gamma);
  Talpha = log(alpha);
  Tiota = log(iota);
  Trho = logit(rho);
  Tcohort = logit(cohort);
  Tamplitude = logit(amplitude);
  TsigmaSE = log(sigmaSE);
  Tpsi = log(psi);
  TR0 = log(R0);
  to_log_barycentric (&TS_0, &S_0, 4);
")

fromEst <- Csnippet("
  Tmu = exp(mu);
  Tsigma = exp(sigma);
  Tgamma = exp(gamma);
  Talpha = exp(alpha);
  Tiota = exp(iota);
  Trho = expit(rho);
  Tcohort = expit(cohort);
  Tamplitude = expit(amplitude);
  TsigmaSE = exp(sigmaSE);
  Tpsi = exp(psi);
  TR0 = exp(R0);
  from_log_barycentric (&TS_0, &S_0, 4);
")


pomp(m1,toEstimationScale=toEst,
     fromEstimationScale=fromEst,
     statenames=c("S","E","I","R","C","W"),
     paramnames=c("R0","mu","sigma","gamma","alpha","iota",
                  "rho","sigmaSE","psi","cohort","amplitude",
                  "S_0","E_0","I_0","R_0")) -> m1

#' 
#' ### Construction of a likelihood profile
#' 
#' [The linked document shows how a likelihood profile can be constructed using IF2](./profile.html).
#' 
#' ## Results from @He2010
#' 
#' ### Fitting procedures
#' 
#' - A large Latin-hypercube design was used to initiate searches.
#' - Iterated filtering was used to maximize the likelihood.
#' - We obtained point estimates of all parameters for 20 cities.
#' - We constructed profile likelihoods to quantify uncertainty in London and Hastings.
#' 
#' 
#' #### Imported infections
#' 
#' $$\text{force of infection} = \mu_{SE}=\frac{\beta(t)}{N(t)}\,(I+\iota)\,\zeta(t)$$
#' 
#' 
#' 
#' #### Seasonality
#' 
#' 
#' 
#' ### Notable findings
#' 
#' #### Cohort effect
#' 
#' 
#' 
#' #### Birth delay
#' 
#' 
#' #### Reporting rate
#' 
#' 
#' #### Predicted vs observed critical community size
#' 
#' 
#' ### Problematic results
#' 
#' #### $R_0$
#' 
#' - Recall that $R_0$ is the basic reproduction number: a measure of how communicable an infection is.
#' - Existing estimates of $R_0$ (c. 15--20) come from two sources:
#' - serology surveys
#' - models fit to data using feature-based methods
#' 
#' 
#' #### Parameter estimates
#' 
#' 
#' $r = \mathrm{cor}(\log{\hat\theta},\log{N_{1950}})$
#' 
#' #### Extrademographic stochasticity
#' 
#' $$\mu_{SE}=\frac{\beta(t)}{N(t)}\,(I+\iota)\,\zeta(t)$$
#' 
#' 
#' ### Questions
#' 
#' 1. What does it mean that parameter estimates from the fitting disagree with estimates from other data?
#' 1. How can one interpret the correlation between infectious period and city size in the parameter estimates?
#' 1. How do we interpret the need for extrademographic stochasticity in this model?
#' 
#' Simulations at the MLE:
#' 
## ----sims2---------------------------------------------------------------
m1 %>% 
  simulate(params=theta,nsim=100,as.data.frame=TRUE,include.data=TRUE) %>%
  subset(select=c(time,sim,cases)) %>%
  mutate(data=sim=="data") %>%
  ddply(~time+data,summarize,
        p=c(0.05,0.5,0.95),q=quantile(cases,prob=p,names=FALSE)) %>%
  mutate(p=mapvalues(p,from=c(0.05,0.5,0.95),to=c("lo","med","hi")),
         data=mapvalues(data,from=c(TRUE,FALSE),to=c("data","simulation"))) %>%
  dcast(time+data~p,value.var='q') %>%
  ggplot(aes(x=time,y=med,color=data,fill=data,ymin=lo,ymax=hi))+
  geom_ribbon(alpha=0.2)

#' 
#' ## Exercises
#' 
#' ### Exercise: Reformulate the model
#' 
#' Modify the @He2010 model to remove the cohort effect.
#' Run simulations and compute likelihoods to convince yourself that the resulting codes agree with the original ones for `cohort = 0`.
#' 
#' Now modify the transmission seasonality to use a sinusoidal form.
#' How many parameters must you use?
#' Fixing the other parameters at their MLE values, compute and visualize a profile likelihood over these parameters.
#' 
#' ### Exercise: Extrademographic stochasticity
#' 
#' Set the extrademographic stochasticity parameter $\sigma_{SE}=0$, set $\alpha=1$, and fix $\rho$ and $\iota$ at their MLE values, then maximize the likelihood over the remaining parameters. 
#' How do your results compare with those at the MLE? 
#' Compare likelihoods but also use simulations to diagnose differences between the models.
#' 
#' --------------------------
#' 
#' ## [Back to course homepage](../index.html)
#' ## [**R** codes for this document](https://raw.githubusercontent.com/kingaa/sbied/master/measles/measles.R)
#' ## [Profile likelihood computation for this example](./profile.html)
#' 
#' ----------------------
#' 
#' ## References
