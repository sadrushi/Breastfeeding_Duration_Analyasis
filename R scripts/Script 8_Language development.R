library(tidyverse)
library(car)

# Loading dataset
cordata <- read_csv("Datasets/Data for language development.csv")

# Factorizing variables
cordata <- cordata %>%
  mutate(Age.at.delivery = factor(Age.at.delivery, levels = c("18-30","31-45")),
         Maternal.educational.level = factor(Maternal.educational.level, levels = c("Up to ordinary level","Up to advanced level","Higher education")),
         Maternal.employment = factor(Maternal.employment, levels = c("Unemployed","Employed")),
         Wealth.status = factor(Wealth.status, levels = c("Lower","Middle","Higher")),
         Child.gender = factor(Child.gender, levels = c("Female","Male")),
         Birth.order = factor(Birth.order, levels = c("First","Second","Third or above")),
         Child.birth.weight = factor(Child.birth.weight, levels = c("Low","Normal")))

# Defining reference levels
cordata$Age.at.delivery <- relevel(cordata$Age.at.delivery, ref = "31-45")
cordata$Maternal.educational.level <- relevel(cordata$Maternal.educational.level,ref = "Higher education")
cordata$Maternal.employment <- relevel(cordata$Maternal.employment,ref = "Employed")
cordata$Wealth.status <- relevel(cordata$Wealth.status,ref = "Higher")
cordata$Child.gender <- relevel(cordata$Child.gender,ref = "Female")
cordata$Birth.order <- relevel(cordata$Birth.order,ref = "First")
cordata$Child.birth.weight <- relevel(cordata$Child.birth.weight, ref = "Low")

# calculating 75th percentile of age of achievement
coo <- quantile(cordata$Makes.cooing.sounds, 0.75)
mimic <- quantile(cordata$Mimics.sounds, 0.75)
pro <- quantile(cordata$Pronounce.individual.words, 0.75)
sen <- quantile(cordata$Make.sentences, 0.75)


# categorizing variables
cordata <- cordata %>% mutate(Breastfeeding.group = ifelse(EBF.Time >=6 & 
                                                             Breastfeeding.time >= 24, "group 1","group 2"))

cordata$Breastfeeding.group <- as.factor(cordata$Breastfeeding.group)
cordata$Breastfeeding.group <- relevel(cordata$Breastfeeding.group, ref = "group 2")


cordata <- cordata %>% mutate(coo2 = ifelse(Makes.cooing.sounds <= coo, 1, 0))
cordata <- cordata %>% mutate(mimic2 = ifelse(Mimics.sounds <= mimic, 1, 0))
cordata <- cordata %>% mutate(indword2 = ifelse(Pronounce.individual.words <= 
                                                  pro, 1, 0))
cordata <- cordata %>% mutate(sentecs2 = ifelse(Make.sentences <= sen, 1, 0))


# Defining domain

cordata <- cordata %>% mutate(Language = ifelse(coo2 == 1 & mimic2 == 1 & 
                                                  indword2 == 1 & sentecs2 == 1,
                                                1, 0))


# Fitting the model
full.language <- glm(Language~Breastfeeding.group+Age.at.delivery+
                       Maternal.educational.level+Wealth.status+
                       Maternal.employment+ Child.birth.weight +
                       Child.gender + Birth.order,
                     data = cordata, family = binomial)

summary(full.language)

x <- exp(cbind(OR = coef(full.language), confint(full.language, level = 0.9)))

round(x, 3)
