## this file simulates some space only data that is like the woman age
## 50-54 dataset for a single year
##
## AOZ Nov 2018

## first, load the real dataset so we can use it to inform our realistics simualtion

###########
## NOTES ##
###########

## 1) age 50-54 is age group 11 


###########
## SETUP ##
###########

## setup working environment, 
## load raw data and prepare two lists:
# 1) data.list: list of all data objects passed into TMB
# 2) init.list: list of initial values for all params
source('~/Desktop/EUCancer/EUCancerCode/0_prep_cancer.R')

## load any extra packages
require(TMB)
require(Matrix)


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

sim.dat <- merge(types, pops, by = 'cname', all.x = TRUE, all.y = FALSE)

## type 2 country registry coverage
type2Morty <- subset(type2Mort, year == real.yr & a == 11)
regDaty    <- subset(regDat, year == real.yr & a == 11)
type2cov   <- data.table(cname = type2Morty$cname, reg.cov = regDaty$pop.m / type2Morty$pop.m)

## set pop.m for type II countries
sim.dat[, pop.m := -1.0]
for(mm in 1:nrow(type2cov)){
  c.row <- which(sim.dat$cname == type2cov$cname[mm])
  sim.dat[c.row, pop.m := pop.i * type2cov[mm, reg.cov]]
}

##############
## SIMULATE ##
##############
set.seed(1989)

## TODO add on missing populations
if(mean(is.na(sim.dat$pop.i)) < 1){
  for(ii in which(is.na(sim.dat$pop.i))){
    sim.dat[ii, pop.i := runif(n = 1, min = min(sim.dat$pop.i, na.rm = T),
                               max = max(sim.dat$pop.i, na.rm = T))]
 }
}

## global intercepts
a.i  <- -6.5
a.mi <- -1.0

## country offsets
v.i.sd  <- 0.5
v.mi.sd <- 0.5
log_tau_v.i  <- log(1 / v.i.sd ^ 2)
log_tau_v.mi <- log(1 / v.mi.sd ^ 2)
v.i  <- rnorm(n = Countries, mean = 0, sd = v.i.sd)
v.mi <- rnorm(n = Countries, mean = 0, sd = v.mi.sd)

## ~~ spatial rand eff ~~ ##

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


## generate the ICAR REs ####
W <- matrix(NA, ncol = Countries, nrow = Countries)
A <- as.matrix(Amat)
m <- apply(Amat,1,sum)
D <- diag( 1 / m )

for(i in 1:Countries){
  #i<-2
  W[i,]<-A[i,] / m[i]
}

s.i.sd  <- 0.5
s.mi.sd <- 0.5
log_tau_s.i  <- log(1 / s.i.sd ^ 2)
log_tau_s.mi <- log(1 / s.mi.sd ^ 2)
s.i  <- ricar(1, s.i.sd , D, W)
s.mi <- ricar(1, s.mi.sd, D, W)

## -- make incidence rates -- ##
p.c <- exp(a.i + v.i + s.i)
sim.dat[, inc.rates := p.c]

## -- make mortality rates -- ##
r.c <- exp(a.mi + v.mi + s.mi) / (1 + exp(a.mi + v.mi + s.mi))
sim.dat[, mi.ratios := r.c]

## -- simulate data -- ##4
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

## -- restructure to feed into TMB -- ##
Data <- list(N        = nrow(sim.dat), ## total countries in dataset
             Ctype    = sim.dat$type, ## vector of types
             pop_nat  = sim.dat$pop.i, ## natl pop
             pop_reg  = sim.dat$pop.m, ## regi pop
             r_cases  = sim.dat$r.cases, ## remainder cases
             r_deaths = sim.dat$r.deaths,## remainder deaths
             l_cases  = sim.dat$l.cases, ## local cases
             l_deaths = sim.dat$l.deaths,## local deaths
             K        = K, ## 
             I        = I, ## identity matrix
             options = c(1, ## include priors?             
                         0) ## adreport?
             )

## -- set initial values for all parameters -- ##
Init <- list(aI  = 0,
             aMI = 0,
             vI  = rep(0, nrow(sim.dat)),
             vMI = rep(0, nrow(sim.dat)),
             log_tau_vI  = 0,
             log_tau_vMI = 0,
             uI  = rep(0, nrow(sim.dat)),
             uMI = rep(0, nrow(sim.dat)),
             lambdaI  = 0.5,
             lambdaMI = 0.5, 
             log_tau_uI  = 0,
             log_tau_uMI = 0            
             )

Init <- list(aI  = -6,
             aMI = -1,
             vI  = rep(0, nrow(sim.dat)),
             vMI = rep(0, nrow(sim.dat)),
             log_tau_vI  = 1.5,
             log_tau_vMI = 1.5,
             uI  = rep(0, nrow(sim.dat)),
             uMI = rep(0, nrow(sim.dat)),
             lambdaI  = 0.99,
             lambdaMI = 0.99, 
             log_tau_uI  = 1.5,
             log_tau_uMI = 1.5            
             )


## -- set lower and upper bounds on params -- ##
L <- rep(-Inf, length(unlist(Init)))
names(L) <- names(unlist(Init))
U <- rep(+Inf, length(unlist(Init)))
names(U) <- names(unlist(Init))
L[which(grepl('lambda', names(L)))] <- 1e-3
U[which(grepl('lambda', names(U)))] <- 1 - 1e-3


L <- rep(-Inf, length(L))

## -- set which params are random -- ##
Rand <- c("vI", "vMI", "uI", "uMI")

## -- set nlminb control params -- ##
Control <- list(eval.max = 500, iter.max = 500)

## -- set vector of true fixed effects -- ##
lambdaI = lambdaMI = 0.99 ## since we simulate from ICAR
true.par <- c(a.i, a.mi,
              log_tau_v.i, log_tau_v.mi,
              lambdaI, lambdaMI,
              log_tau_s.i, log_tau_s.mi)

## -- compile TMB -- ##
TMB::compile('~/Desktop/EUCancer/EUCancerCode/cancer_space_sim_model.cpp')
dyn.load(dynlib('~/Desktop/EUCancer/EUCancerCode/cancer_space_sim_model'))

## -- make autodiff likelihood functions -- ##
obj <- MakeADFun(data = Data,
                 parameters = Init,
                 random = Rand, 
                 hessian = TRUE,
                 DLL = 'cancer_space_sim_model')
obj$env$inner.control$tol10 <- 0

## run optimization
opt0 <- nlminb(start     = obj$par,
               objective = obj$fn,
               gradient  = obj$gr,
               lower     = L, 
               upper     = U, 
               control = Control
               )


## run optimization
opt0 <- optim(par = obj$par,
              fn  = obj$fn,
              gr  = obj$gr,
              lower     = L, 
              upper     = U,
              method = 'L-BFGS-B'#, 
              #control = Control
              )

comp <- data.frame(truth = c(a, u, log(sobs)),
                   est = opt0$par[1:(ngrp + 2)])
rownames(comp) <- c('a', paste0('u', 1:ngrp), 'logSobs')
print(comp)
message(sprintf('Sum of u terms is: %.9f', sum(opt0$par[2:(1 + ngrp)])))

## run SD report
SD0 <- TMB::sdreport(obj, getJointPrecision = TRUE)
