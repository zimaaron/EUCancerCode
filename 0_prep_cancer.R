## this script preps things for the eu cancer model
## aoz - oct 2018

###########
## NOTES ##
###########

## DATA TYPES
# type 1: national inc and mort
# type 2: subnat inc & mort and natl mort
# type 3: only natl mort
# type 4: no data


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

## set paths
setwd('~/Desktop/EUCancer/')
data.dir <- paste0(getwd(), '/Data/')
code.dir <- paste0(getwd(), '/EUCancerCode/')

## load packages
require(data.table)
require(dplyr)

#######################
## LOAD DATA OBJECTS ##
#######################

load(paste0(data.dir, 'AdjacencyMatrix.rdata')) ## Amat
amat    <- Amat
bii     <- fread(paste0(data.dir,'BreastIarcIhme.csv'))
nat     <- fread(paste0(data.dir,'C50_national_1990_2010_Aggregated_20160515.csv'))
reg     <- fread(paste0(data.dir,'C50_registry_1990_2010_Aggregated_20160515.csv'))
qual    <- fread(paste0(data.dir, 'Quality_info.csv'))

## take a look at the datasets
if(FALSE){

  str(bii)
  head(bii)

  str(nat)
  nat
  table(nat[,cname], nat[,year])
  unique(nat[, type])
  ## complete across country years
  ## NOTE: the reason there are 17 entries per country/year is the 17 age groups
  ## NOTE: data types 1-4 exist in nat
  
  str(reg)
  reg
  table(reg[,cname], reg[,year])
  reg[cname == 'Bosnia', ]
  unique(reg[, type])
  ## note complete - only where there exists data
  ## only data type 2 exist in reg

  qual
  
}


###################################
###################################
##                               ##
## PROCESS AND PREP DATA OBJECTS ##
##                               ##
###################################
###################################


#-- create K (structure) matrix --#
k.u <- as.matrix((-1)*amat)
m <- apply(amat, 1, sum)
diag(k.u) <- m
qr(k.u)$rank # 39

# K and I are passed to the TMB in the data list
K <- diag(apply(Amat, 1, sum)) - Amat
I <- diag(Countries)


# -- reorder the national data -- #
natsub <- nat %>% arrange(id)

# -- reorder the registry data -- #
regsub <- reg %>% arrange(id)


#################################################
## create the Y matrix, type 1 and 2 countries ##
#################################################

# -- Type 1 Incidence: National Incidence data available -- #
type1 <- natsub %>% filter(!is.na(cases))

# make sure there seem to be the rigth number of types
with(type1, table(type))
with(natsub, table(type))
with(type1, table(as.character(cname)))

# -- Type 2 Incidence: Registry Incidence data available -- #
# reorder the columns of the registry data#
regsub <- regsub %>% select(names(type1))

# -- Combining Type 1 & 2 Incidence -- #
type12 <- rbind(type1, regsub) %>% arrange(id)
with(type12, table(type))


# -- creating vectors to be used in Stan code -- #
Y <- type12$cases
popY <- type12$pop.m

# - used for indices - #
Ya <- as.numeric(as.factor(type12$a))
Yt <- type12$t
Yc <- type12$c

# - Deaths to accompany cases - #
Z <- type12$deaths

# - the number of Type 1 & 2 Countries - #
NY <- length(Y)


# ---- Occasionally there are more deaths than cases (due to lag in mortality)
# ---- For these cases we currently set deaths=cases.
# ---- This will likely need a different solution for other cancers,
# ---- but it is pretty rare in breast cancer

# if Y < Z, replace Z with Y
table(Y < Z)
table(Ya[Y < Z])
table(Yc[Y < Z])
table(Yt[Y < Z])
if(FALSE){
  plot(Y[which(Y < Z)],Z[which(Y < Z)])
  abline(0, 1, lty=2)
}

summary((Z[which(Y < Z)] - Y[which(Y < Z)]) / Z[which(Y < Z)])
length(Z[which(Y < Z)])
table(Ya[which(Y < Z)])

# If there are more deaths than cases, then increase the number of cases
Y[which(Y < Z)] <- Z[which(Y < Z)]
table(Y < Z)


###########################################################################################################
## Extract Mortality Data for Type 2 (remainder not coverged by registry, nat-reg) and Type 3 (national) ##
###########################################################################################################


# ---- Deaths for type 2 ----- #
type2Mort <- natsub %>% filter(type==2) %>% arrange(id)

#
dim(type2Mort)
dim(regsub)

# sometimes there is registry data in country/years for which there is complete
# national data. Thus the difference in dimension here.  In those cases I limit
# the registry data to only those for which there is no national incidence data
# as we are better off modeling just national when it is available.  We assume it is free
# from the possible bias that may exist in registries.  However, in exploratory
# analyses the registries and national rates seem to match up nicely.

# -- grab the appropriate registry data -- #
regDat <- regsub %>% filter(id %in% type2Mort$id) %>% arrange(id)
dim(regDat)
dim(type2Mort)


# -- isolate just the 'remainder' mortality, that which is not covered by registries -- #
type2Deaths <- type2Mort$deaths - regDat$deaths

# -- remainder population -- #
type2pop <- type2Mort$pop.m - regDat$pop.m

# -- needed for indices in the MCMC -- #
type2Deaths_a <- as.numeric(as.factor(regDat$a))
type2Deaths_c <- regDat$c
type2Deaths_t <- regDat$t

# ---- Type 3 Deaths: National Mortality (with no available incidence) --- #
type3Mort <- natsub %>% filter(type==3) %>% arrange(id)

type3Deaths <- type3Mort$deaths
type3pop <- type3Mort$pop.m

# -- indices for MCMC -- #
type3Deaths_a <- as.numeric(as.factor(type3Mort$a))
type3Deaths_c <- type3Mort$c
type3Deaths_t <- type3Mort$t

# --- Combine Type 2 and 3 Deaths to be used in the Mortality Model in the MCMC --- #
Deaths <- c(type2Deaths, type3Deaths)
Deaths[Deaths < 0] <- 0
popDeaths <- c(type2pop, type3pop)
Deaths_a  <- c(type2Deaths_a, type3Deaths_a)
Deaths_t  <- c(type2Deaths_t, type3Deaths_t)
Deaths_c  <- c(type2Deaths_c, type3Deaths_c)

# -- used for the loop -- #
Ndeaths <- length(Deaths)

# --- Finally for the random effects how many factors are there for country/year/age --- #
# For the Random Effects #
Countries <- with(natsub, max(cnum))
Years     <- with(natsub, max(t)) # go out to 2012 to compare with IARC
Ages      <- with(natsub, length(unique(a)))

index <- data.frame(year=rep(1:Years, each = Countries * Ages))
index$country <- rep(1:Countries, each=Ages, times=Years)
index$age     <- rep(1:Ages, times = Years * Countries)
index <- index %>% arrange(country, year, age)

Nrates  <- nrow(index)
index_a <- index$age
index_c <- index$country
index_t <- index$year

#################################################
#################################################
##                                             ##
## COMBINE PREPPED OBJECTS INTO A LIST FOR TMB ##
##                                             ##
#################################################
#################################################

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
                  
                  Nrates=Nrates,
                  index_c=index_c,
                  index_t=index_t,

                  index_a=index_a)

init.list <- list(
  list(
    betI = (0.0),betMI=(0.0),

    bI = rep(0, Countries),
    bMI = rep(0, Countries),
    betaI = rep(0, Countries),
    betaMI = rep(0, Countries),

    gammaI = rep(-9.4, Ages),
    gammaMI = rep(-0.8, Ages),
    deltaI = matrix(rep(0.0, Ages*Countries),ncol=Ages),
    deltaMI = matrix(rep(0.0, Ages*Countries),ncol=Ages),
    tau_bi = 100,
    tau_bmi = 100,
    tau_betai = 100,
    tau_betami = 100,
    tau_delti = 100,
    tau_deltmi = 100,

    tau_gami = 100,
    tau_gammi = 100,
    lambdaI = 0.5,
    lambdaMI = 0.5
  )
) 

