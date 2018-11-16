## this script is used to launch the eu cancer model
## aoz oct 2018

###########
## NOTES ##
###########

## DATA TYPES
# type 1: national inc and mort
# type 2: subnat inc & mort and natl mort
# type 3: only natl mort
# type 4: no data


## TODO - I copied existing code notation which is not as written below!
## NOTATION - consistent with paper
# a : age group
# t : time
# c : country

# .l : local area in country
# .r : remainder population not covered in .l local registries

# n.l : population at a,t,c (total in all local registries)
# y.l : incident reported cases at a,t,c (total of all registries)
# z.l : mortality (reported deaths) at a,t,c (total of all registries)
# n.r : remainder population at a,t,c (pop not covered by registries)
# y.r : incident reported cases at a,t,c in remainder pop (not covered by registries)
# z.r : mortality (reported deaths) at a,t,c in remainder pop (not covered by registries)

# y = y.l + y.r : all reported cases in a,c,t
# z = z.l + z.r : all reported deaths in a,c,t
# n = n.l + n.r : total pop in a,c,t

# p : Prob(reported incidence | a,c,t) 
# q : Prob(reported mortality | a,c,t) 
# r : Prob(reported mortality | reported incidence, a,c,t)



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

########################
## PREPARE TO RUN TMB ##
########################


# K and I are passed to the TMB in the data list
K <- Matrix(diag(apply(Amat, 1, sum)) - Amat, sparse = TRUE) ## TMB expects sparse mat
I <- Matrix(diag(Countries), sparse = TRUE) ## TMB expects sparse mat

data.list <- list(NY=NY, 
                  popY=popY,
                  Y=Y,
                  Z=Z,
                  Yc=Yc,
                  Ya=Ya,
                  Yt=Yt,

                  Ndeaths=Ndeaths,
                  popDeaths=popDeaths,
                  Deaths=Deaths,
                  Deaths_c=Deaths_c,
                  Deaths_a=Deaths_a,
                  Deaths_t=Deaths_t,

                  Countries=Countries,
                  Ages=Ages,
                  K=K,
                  I=I,

                  ## ## for 2012 projections
                  ## N2012=N2012,
                  ## iarc_c=iarc_c,
                  ## iarc_t=iarc_t,
                  ## iarc_a=iarc_a,
                  
                  ## Nrates=Nrates,
                  index_c=index_c,
                  index_t=index_t,
                  index_a=index_a,

                  options = c(1) ## LOGICAL: include priors?
                  )

param.list <- list(
  betI = (0.0),
  betMI=(0.0),

  bI = rep(0, Countries),
  bMI = rep(0, Countries),
  betaI = rep(0, Countries),
  betaMI = rep(0, Countries),
  
  gammaI = rep(-9.4, Ages),
  gammaMI = rep(-0.8, Ages),
  deltaI = matrix(rep(0.0, Ages * Countries), ncol=Ages),
  deltaMI = matrix(rep(0.0, Ages * Countries), ncol=Ages),
  log_tau_bi = log(100),
  log_tau_bmi = log(100),
  log_tau_betai = log(100),
  log_tau_betami = log(100),
  log_tau_delti = log(100),
  log_tau_deltmi = log(100),
  
  log_tau_gami = log(100),
  log_tau_gammi = log(100),
  lambdaI = 0.5, ## TODO must be between (0,1)
  lambdaMI = 0.5 ## TODO must be between (0,1)
  
)


## TODO bounds for params
## L <- c(-Inf, rep(-Inf, ngrp), -Inf)
## U <- c(Inf, rep(Inf, ngrp), Inf)

## load TMB and compile model
TMB::compile('~/Desktop/EUCancer/EUCancerCode/cancer_model.cpp')
dyn.load(dynlib('~/Desktop/EUCancer/EUCancerCode/cancer_model'))

## make autodiff likelihood functions
obj <- MakeADFun(data = data.list,
                 parameters = param.list, 
                 hessian = TRUE,
                 DLL = 'cancer_model')

## run optimization
opt0 <- nlminb(start     = obj$par,
               objective = obj$fn,
               gradient  = obj$gr,
               lower     = L, 
               upper     = U)

comp <- data.frame(truth = c(a, u, log(sobs)),
                   est = opt0$par[1:(ngrp + 2)])
rownames(comp) <- c('a', paste0('u', 1:ngrp), 'logSobs')
print(comp)
message(sprintf('Sum of u terms is: %.9f', sum(opt0$par[2:(1 + ngrp)])))

## run SD report
SD0 <- TMB::sdreport(obj, getJointPrecision = TRUE)
