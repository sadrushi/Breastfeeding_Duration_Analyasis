library(tidyverse)
library(meta)
library(readxl)

settings.meta('RevMan5')

###############################################################################

Maternal_age <- read_excel("Datasets/Maternal age.xlsx", 
                           sheet = "Final data")

## run function metagen to get estimates using generic inverse variance method  
mR <- metagen( 
  TE = log(AdjustedHR), lower = log(CILower), upper = log(CIHigher),  
  sm = "HR", common=F, random=T, 
  studlab = `Article No`, data = Maternal_age,
  method.tau = "DL",   ## method to calculate Tau 
  method.random.ci = "classic",  ## method to calculate estimator's CI 
) 

summary(mR) 

###############################################################################

Education.primary <- read_excel("Datasets/Education.xlsx", 
                                sheet = "Final data - primary")

## run function metagen to get estimates using generic inverse variance method  
mR <- metagen( 
  TE = log(Hazardratio), lower = log(CILower), upper = log(CIHigher),  
  sm = "HR", common=F, random=T, 
  studlab = `Article No`, data = Education.primary,
  method.tau = "DL",   ## method to calculate Tau 
  method.random.ci = "classic",  ## method to calculate estimator's CI 
) 

summary(mR) 

Education.secondary <- read_excel("Datasets/Education.xlsx", 
                                  sheet = "Final data - secondary")

## run function metagen to get estimates using generic inverse variance method  
mR <- metagen( 
  TE = log(`Hazard ratio`), lower = log(`CI Lower`), upper = log(`CI Higher`),  
  sm = "HR", common=F, random=T, 
  studlab = `Article No`, data = Education.secondary,
  method.tau = "DL",   ## method to calculate Tau 
  method.random.ci = "classic",  ## method to calculate estimator's CI 
) 

summary(mR) 

###############################################################################

Employment <- read_excel("Datasets/Employment.xlsx", 
                         sheet = "Final data")

## run function metagen to get estimates using generic inverse variance method  
mR <- metagen( 
  TE = log(`Hazard ratio`), lower = log(`CI Lower`), upper = log(`CI Higher`),  
  sm = "HR", common=F, random=T, 
  studlab = `Article No`, data = Employment,
  method.tau = "DL",   ## method to calculate Tau 
  method.random.ci = "classic",  ## method to calculate estimator's CI 
) 

summary(mR) 

###############################################################################

Delivery_mode <- read_excel("Datasets/Delivery mode.xlsx", 
                            sheet = "Final data")

## run function metagen to get estimates using generic inverse variance method  
mR <- metagen( 
  TE = log(Hazardratio), lower = log(CILower), upper = log(CIHigher),  
  sm = "HR", common=F, random=T, 
  studlab = `Article No`, data = Delivery_mode,
  method.tau = "DL",   ## method to calculate Tau 
  method.random.ci = "classic",  ## method to calculate estimator's CI 
) 

summary(mR) 

###############################################################################

Wealth_status.poor <- read_excel("Datasets/Wealth status.xlsx", 
                                 sheet = "Final data - poor")

## run function metagen to get estimates using generic inverse variance method  
mR <- metagen( 
  TE = log(`Hazard ratio`), lower = log(`Lower CI`), upper = log(`Higher CI`),  
  sm = "HR", common=F, random=T,
  studlab = `Article No`, data = Wealth_status.poor,
  method.tau = "DL",   ## method to calculate Tau 
  method.random.ci = "classic",  ## method to calculate estimator's CI 
) 

summary(mR) 

Wealth_status.middle <- read_excel("Datasets/Wealth status.xlsx", 
                                   sheet = "Final data - middle")

## run function metagen to get estimates using generic inverse variance method  
mR <- metagen( 
  TE = log(`Hazard ratio`), lower = log(`Lower CI`), upper = log(`Higher CI`),  
  sm = "HR", common=F, random=T, 
  studlab = `Article No`, data = Wealth_status.middle,
  method.tau = "DL",   ## method to calculate Tau 
  method.random.ci = "classic",  ## method to calculate estimator's CI 
) 

summary(mR) 

###############################################################################

Child_s_Gender <- read_excel("Datasets/Child's Gender.xlsx", 
                             sheet = "Final data")

## run function metagen to get estimates using generic inverse variance method  
mR <- metagen( 
  TE = log(Hazardratio), lower = log(CILower), upper = log(CIHigher),  
  sm = "HR", common=F, random=T, 
  studlab = `Article No`, data = Child_s_Gender,
  method.tau = "DL",   ## method to calculate Tau 
  method.random.ci = "classic",  ## method to calculate estimator's CI 
) 

summary(mR) 

###############################################################################

Birth_order.second <- read_excel("Datasets/Birth order.xlsx", 
                                 sheet = "Final data - second")

## run function metagen to get estimates using generic inverse variance method  
mR <- metagen( 
  TE = log(`Hazard ratio`), lower = log(`CI Lower`), upper = log(`CI Higher`),  
  sm = "HR", common=F, random=T,
  studlab = `Article No`, data = Birth_order.second,
  method.tau = "DL",   ## method to calculate Tau 
  method.random.ci = "classic",  ## method to calculate estimator's CI 
) 

summary(mR) 

Birth_order.third <- read_excel("Datasets/Birth order.xlsx", 
                                sheet = "Final data - third")

## run function metagen to get estimates using generic inverse variance method  
mR <- metagen( 
  TE = log(`Hazard ratio`), lower = log(`CI Lower`), upper = log(`CI Higher`),  
  sm = "HR", common=F, random=T,
  studlab = `Article No`, data = Birth_order.third,
  method.tau = "DL",   ## method to calculate Tau 
  method.random.ci = "classic",  ## method to calculate estimator's CI 
) 

summary(mR) 


###############################################################################

Birth_weight <- read_excel("Datasets/Birth weight.xlsx", 
                           sheet = "Final data")

## run function metagen to get estimates using generic inverse variance method  
mR <- metagen( 
  TE = log(AdjustedHR), lower = log(`AdjustedCI L`), upper = log(`AdjustedCI H`),  
  sm = "HR", common=F, random=T, 
  studlab = `Article No`, data = Birth_weight,
  method.tau = "DL",   ## method to calculate Tau 
  method.random.ci = "classic",  ## method to calculate estimator's CI 
) 

summary(mR) 

###############################################################################

Age_at_feeding <- read_excel("Datasets/Age at feeding.xlsx", 
                             sheet = "Final data")

## run function metagen to get estimates using generic inverse variance method  
mR <- metagen( 
  TE = log(`Hazard ratio`), lower = log(`CI Lower`), upper = log(`CI Upper`),  
  sm = "HR", common=F, random=T,
  studlab = `Article No`, data = Age_at_feeding,
  method.tau = "DL",   ## method to calculate Tau 
  method.random.ci = "classic",  ## method to calculate estimator's CI 
) 

summary(mR) 

