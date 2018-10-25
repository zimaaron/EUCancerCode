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

########################
## PREPARE TO RUN TMB ##
########################
