# loading libraries
library(tidyverse)
library(gridExtra)
library(DescTools)
library(readxl)

# Importing EDA_Data sets
EDA_Data <- read_excel("Datasets/Data for BF duration analysis.xlsx")


#---------------------------------------------------------

# Maternal age at delivery
table(EDA_Data$Age.at.delivery)
round(prop.table(table(EDA_Data$Age.at.delivery))*100, digits = 1)

# Maternal educational level
table(EDA_Data$Maternal.educational.level)
round(prop.table(table(EDA_Data$Maternal.educational.level))*100, digits = 1)

# Maternal employment status
table(EDA_Data$Maternal.employment)
round(prop.table(table(EDA_Data$Maternal.employment))*100, digits = 1)

# Wealth status
table(EDA_Data$Wealth.status)
round(prop.table(table(EDA_Data$Wealth.status))*100, digits = 1)

# Delivery mode
table(EDA_Data$Delivery.mode)
round(prop.table(table(EDA_Data$Delivery.mode))*100, digits = 1)

# Child's Birth order
table(EDA_Data$Birth.order)
round(prop.table(table(EDA_Data$Birth.order))*100, digits = 1)

# Child's gender
table(EDA_Data$Child.gender)
round(prop.table(table(EDA_Data$Child.gender))*100, digits = 1)

# Child's age at first supplementary feeding
table(EDA_Data$Age.at.first.feeding)
round(prop.table(table(EDA_Data$Age.at.first.feeding))*100, digits = 1)

# Child's birth weight
table(EDA_Data$Child.birth.weight)
round(prop.table(table(EDA_Data$Child.birth.weight))*100, digits = 1)

# Child's age in years
EDA_Data <- EDA_Data %>% mutate(Child.age.years = round(Child.age/12,1))
summary(EDA_Data$Child.age.years)
MeanCI(EDA_Data$Child.age.years, level = 0.95)
