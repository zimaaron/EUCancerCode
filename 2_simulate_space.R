## this file simulates some space only data that is like the woman age
## 50-54 dataset for a single year
##
## AOZ Nov 2018

## first, load the real dataset so we can use it to inform our realistics simualtion

###########
## NOTES ##
###########

## 1) age 50-54 is age group 11 


## TODO!
## simulate in space?
space <- TRUE

###########
## SETUP ##
###########

## setup working environment, 
## load raw data and prepare two lists:
# 1) data.list: list of all data objects passed into TMB
# 2) init.list: list of initial values for all params
source('~/Desktop/EUCancer/EUCancerCode/0_prep_cancer.R')

##########################
## USE REALISTIC VALUES ##
##########################
real.yr <- 2008

## -- get pops -- ##
pops  <- aggregate(pop.i ~ cname, nat, mean)
pops  <- rbind(pops,aggregate(pop.i ~ cname, reg, mean))

## -- assign data types to countries -- ##
types <- subset(nat, year == real.yr & age == 1)[, .(cname, type)]
table(types$type)

sim.dat.setup <- merge(types, pops, by = 'cname', all.x = TRUE, all.y = FALSE)

## type 2 country registry coverage
type2Morty <- subset(type2Mort, year == real.yr & a == 11)
regDaty    <- subset(regDat, year == real.yr & a == 11)
type2cov   <- data.table(cname = type2Morty$cname, reg.cov = regDaty$pop.m / type2Morty$pop.m)

## set pop.m for type II countries
sim.dat.setup[, pop.m := -1.0]
for(mm in 1:nrow(type2cov)){
  c.row <- which(sim.dat.setup$cname == type2cov$cname[mm])
  sim.dat.setup[c.row, pop.m := pop.i * type2cov[mm, reg.cov]]
}

##############
## SIMULATE ##
##############

## function to simulate - includes filling in missing populations
## TODO mocve more things into arguments
sim.space.data <- function(a.i  = -6.5, ## incidence intercept
                           a.mi = -1.0, ## mi.ratio  intercept
                           ctry.re = TRUE, 
                           v.i.sd  = 0.5, ## incidence country RE SD
                           v.mi.sd = 0.5, ## mi.ratio  country RE SD 
                           space = TRUE,
                           s.i.sd  = 0.5, ## incidence country SPATIAL RE SD
                           s.mi.sd = 0.5, ## mi.ratio  country SPATIAL RE SD 
                           seed = NULL,
                           dat.template = sim.dat.setup, 
                           plot = TRUE
                           ){

  ## global intercepts
  a.i  <<- a.i  ## set globally so we can compare estiamted vs true
  a.mi <<- a.mi ## set globally so we can compare estiamted vs true

  sim.dat <- copy(dat.template)

  if(!is.null(seed)){
    set.seed(seed)
  }

  ## TODO add on missing populations
  if(mean(is.na(sim.dat$pop.i)) < 1){
    for(ii in which(is.na(sim.dat$pop.i))){
      sim.dat[ii, pop.i := runif(n = 1, min = min(sim.dat$pop.i, na.rm = T),
                                 max = max(sim.dat$pop.i, na.rm = T))]
    }
  }

  ## country offsets
  log_tau_v.i  <<- log(1 / v.i.sd ^ 2)  ## set these globaly so we have the true value for comparison
  log_tau_v.mi <<- log(1 / v.mi.sd ^ 2) ## set these globaly so we have the true value for comparison
  v.i  <- rnorm(n = Countries, mean = 0, sd = v.i.sd)
  v.mi <- rnorm(n = Countries, mean = 0, sd = v.mi.sd)
  ## make them mean zero
  v.i <- v.i - mean(v.i)
  v.mi <- v.mi - mean(v.mi)

  
  ## ~~ spatial rand eff simulation ~~ ##

  ## Lemma 3.1 in Besag & Kooperberg 1995 ##
  ricar <- function( nsim, omega, D, W ){
    require(MASS)
    n <- dim(W)[1]
    Q <- solve(D) %*% (diag(n) - W) / omega ^ 2
    Ainvt <- cbind( diag(1, n - 1), cbind(rep(0, n - 1)) )
    Ainv <- t(Ainvt)
    Qy <- Ainvt %*% Q %*% Ainv 
    Vy <- solve(Qy)
    y <- mvrnorm(n = 1, mu = rep(0, n - 1 ), Sigma = Vy)
    x <- Ainv %*% y
    return(x)
  }

  ## generate the ICAR REs ##
  W <- matrix(NA, ncol = Countries, nrow = Countries)
  A <- as.matrix(Amat)
  m <- apply(Amat, 1, sum)
  D <- diag( 1 / m )

  for(i in 1:Countries){
    W[i,] <- A[i, ] / m[i]
  }

  log_tau_s.i  <<- log(1 / s.i.sd ^ 2)  ## set these globaly so we have the true value for comparison
  log_tau_s.mi <<- log(1 / s.mi.sd ^ 2) ## set these globaly so we have the true value for comparison
  s.i  <- ricar(1, s.i.sd , D, W)
  s.mi <- ricar(1, s.mi.sd, D, W)

  ## -- make incidence rates and mortality-incidence ratios -- ##
  p.c <- a.i  ## incidence intercept 
  r.c <- a.mi ## mi.ratio  intercept
  
  if(ctry.re){
    message('\n-- Simulating p.c and r.c with country random effects\n')
    p.c <- p.c + v.i
    r.c <- r.c + v.mi
  }
  
  if(space){
    message('\n-- Simulating p.c and r.c with spatial effects\n')
    p.c <- p.c + s.i
    r.c <- r.c + s.i
  }

  ## convert to appropriate space with link functions
  p.c <- exp(p.c)
  r.c <- exp(r.c) / (1 + exp(r.c))
  
  sim.dat[, p.c := p.c]
  sim.dat[, r.c := r.c]

  ## -- now we can simulate data for each of the country types -- ##
  for(cc in 1:nrow(sim.dat)){
    
    ## -- simulate type 1 countries: national mortality and registry data -- ##
    if(sim.dat[cc, type] == 1){
      
      ## (r)emainder country cases and deaths.
      ## since there is no (l)ocal info, the (r)emainder is the national estimates
      sim.dat[cc, r.cases  := rpois(n = 1, lambda = pop.i * p.c)]
      sim.dat[cc, r.deaths := rbinom(n = 1, size = r.cases, prob = r.c)]
    }  

    ## -- simulate type 2 countries: local mortality and registry data & complete national mortality -- ##
    if(sim.dat[cc, type] == 2){
      
      ## simulate the (l)ocal info
      sim.dat[cc, l.cases  := rpois(n = 1, lambda = pop.m * p.c)]
      sim.dat[cc, l.deaths := rbinom(n = 1, size = l.cases, prob = r.c)]

      ## and simulate the (r)emainder for the country
      sim.dat[cc, r.deaths := rpois(n = 1, lambda = (pop.i - pop.m) * p.c * r.c)]

    }  

    ## -- simulate type 3 countries: national mortality data -- ##
    if(sim.dat[cc, type] == 3){

      ## since there is no (l)ocal info, the (r)emainder is the national estimates
      sim.dat[cc, r.deaths := rpois(n = 1, lambda = pop.i * p.c * r.c)]
      
    }  
    
    
    ## -- simulate type 4 countries: national mortality and registry data -- ##
    if(sim.dat[cc, type] == 4){

      ## no data    
    }  

  }


  ## check the simulated data
  if(plot){
    par(mfrow = c(1, 2))
    plot(sim.dat$p.c, sim.dat$r.cases / sim.dat$pop.i, main = 'Incidence',
         ylab = "MLE estimate", xlab = 'Truth')
    v <- sim.dat$r.cases / sim.dat$pop.i
    if(space | ctry.re){
      abline(lm(v ~ p.c, data = sim.dat))
      abline(a = 0, b = 1, col = 2)
      legend("topleft", legend = c('y=x', 'lm()'), col = 2:1, lwd = 2)
    }else{
      abline(h = sim.dat$p.c[1])
      legend("topleft", legend = 'truth', lwd = 2)
    }
    
    plot(sim.dat$r.c, sim.dat$r.deaths / sim.dat$r.cases, main = 'Mortality',
         ylab = "MLE estimate", xlab = 'Truth')
    v <- sim.dat$r.deaths / sim.dat$r.cases
    if(space | ctry.re){
      abline(lm(v ~ r.c, data = sim.dat))
      abline(a = 0, b = 1, col = 2)
      legend("topleft", legend = c('y=x', 'lm()'), col = 2:1, lwd = 2)
    }else{
      abline(h = sim.dat$r.c[1])
      legend("topleft", legend = 'truth', lwd = 2)
    }
    
  }
  
  return(sim.dat)
}

## ----------------------------------------------------- ##
## -- Peform simulation and plot the daa vs the truth -- ##
## ----------------------------------------------------- ##
with.ctry.re <- FALSE
with.space   <- FALSE

sim.dat <- sim.space.data(space = with.space, ctry.re = with.ctry.re, seed = 1988)


## -- restructure to feed into TMB -- ##
Data <- list(N        = nrow(sim.dat), ## total countries in dataset
             Ctype    = sim.dat$type , ## vector of types
             pop_nat  = sim.dat$pop.i, ## natl pop
             pop_reg  = sim.dat$pop.m, ## regi pop
             r_cases  = sim.dat$r.cases, ## remainder cases
             r_deaths = sim.dat$r.deaths,## remainder deaths
             l_cases  = sim.dat$l.cases, ## local cases
             l_deaths = sim.dat$l.deaths,## local deaths
             K        = K, ## 
             I        = I, ## identity matrix
             options = c(1, ## include priors?
                         0, ## model with country REs?
                         0, ## model with spatial effect?
                         0) ## adreport?
             )

## -- set initial values for all parameters -- ##
Init <- list(aI  = 0,
             aMI = 0)

Init.ctry.re <- list(vI  = rep(0, nrow(sim.dat)),
                     vMI = rep(0, nrow(sim.dat)),
                     log_tau_vI  = 0,
                     log_tau_vMI = 0)

Init.space <- list(uI  = rep(0, nrow(sim.dat)),
                   uMI = rep(0, nrow(sim.dat)),
                   lambdaI  = 0.5,
                   lambdaMI = 0.5, 
                   log_tau_uI  = 0,
                   log_tau_uMI = 0)

## even if we do not want to fit these params, TMB still needs them in the list
Init <- c(Init, Init.ctry.re, Init.space)

## -- set which parameters are random effefts -- ##
Rand <- NULL

## -- set which params will not be needed in this model run (due to selected options) -- ##
ADmap_list <- list()

## -- make a vector of true (non-random effects) params -- ##
true.par <- c(a.i, a.mi)

## add to the base model obbjects things related to country random effects
if(with.ctry.re){
  Data[['options']][2] <- 1
  Rand <- c(Rand, "vI", "vMI")
  true.par <- c(true.par, c(log_tau_v.i, log_tau_v.mi))
}else{
  tmp.list <- Init.ctry.re
  for(i in 1:length(tmp.list)){
    tmp.list[[i]] <- rep(factor(NA), length = length(tmp.list[[i]]))
  }
  ADmap_list <- c(ADmap_list, tmp.list)
}

## add to the base model obbjects things related to spatial random effect
if(with.space){
  Data[['options']][3] <- 1
  Rand <- c(Rand, "uI", "uMI")
  lambdaI <- lambdaMI <- 0.99 ## since we simulate from ICAR which is really lambda=1
  true.par <- c(true.par, c(lambdaI, lambdaMI,
                            log_tau_s.i, log_tau_s.mi))
}else{
  tmp.list <- Init.space
  for(i in 1:length(tmp.list)){
    tmp.list[[i]] <- rep(factor(NA), length = length(tmp.list[[i]]))
  }
  ADmap_list <- c(ADmap_list, tmp.list)
}

## -- set lower and upper bounds on params -- ##
L <- rep(-Inf, length(unlist(Init)))
names(L) <- names(unlist(Init))
U <- rep(+Inf, length(unlist(Init)))
names(U) <- names(unlist(Init))
L[which(grepl('lambda', names(L)))] <- 1e-3
U[which(grepl('lambda', names(U)))] <- 1 - 1e-3

## -- set nlminb control params -- ##
Control <- list(eval.max = 500, iter.max = 500)

## -- now, all objects are prepped: compile TMB -- ##
TMB::compile('~/Desktop/EUCancer/EUCancerCode/cancer_space_sim_model.cpp')
dyn.load(dynlib('~/Desktop/EUCancer/EUCancerCode/cancer_space_sim_model'))

## -- make autodiff likelihood functions -- ##
obj <- MakeADFun(data = Data,
                 parameters = Init,
                 random = Rand,
                 map = ADmap_list, 
                 hessian = TRUE,
                 DLL = 'cancer_space_sim_model')
## obj$env$inner.control$tol10 <- 0 ## helps avoid early exiting

## run optimization
opt0 <- nlminb(start     = obj$par,
               objective = obj$fn,
               gradient  = obj$gr,
               lower     = L, 
               upper     = U, 
               control = Control
               )

opt0
obj$fn(true.par)


## run optimization
## opt0 <- optim(par = obj$par,
##               fn  = obj$fn,
##               gr  = obj$gr,
##               lower     = L, 
##               upper     = U,
##               method = 'L-BFGS-B'#, 
##               #control = Control
##               )

comp <- data.frame(truth = c(a.i, a.mi, log_tau_v.i, log_tau_v.mi), 
                   est = opt0$par[1:4])

print(comp)
message(sprintf('Sum of u terms is: %.9f', sum(opt0$par[2:(1 + ngrp)])))

## run SD report
SD0 <- TMB::sdreport(obj, getJointPrecision = TRUE)

lower <- SD0$par.fixed - 2 * sqrt(diag(SD0$cov.fixed))
upper <- SD0$par.fixed + 2 * sqrt(diag(SD0$cov.fixed))
mean(c(a.i, a.mi) > lower & c(a.i, a.mi) < upper)


#############
## SCRATCH ##
#############

## ## simulate poisson and look at MLEs
## par(mfrow = c(4, 4))
## for(i in 1:16){
##   pop <- 10000
##   n.sim <- 12
##   true.rate <- runif(n = n.sim, min = 0, max = .003)
##   sim.count <- rpois(n = n.sim, lambda = true.rate * pop)
##   plot(true.rate, sim.count / pop);abline(a = 0, b = 1, col = 'red', lwd = 2) ;abline(lm(sim.count / pop ~ true.rate))
## }













