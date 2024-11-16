library(tidyverse)
library(car)

# Loading data
cordata <- read_csv("Datasets/Data for gross motor development.csv")

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
roll <- quantile(cordata$Rolls.from.tummy.to.back, 0.75)
sit <- quantile(cordata$Sits.without.support, 0.75)
walk <- quantile(cordata$Walks.without.support, 0.75)
climb <- quantile(cordata$Climbs.up.stairs.with.help, 0.75)


# categorizing variables
cordata <- cordata %>% mutate(Breastfeeding.group = ifelse(EBF.Time >=6 & 
                                                             Breastfeeding.time >= 24, "group 1","group 2"))

cordata$Breastfeeding.group <- as.factor(cordata$Breastfeeding.group)
cordata$Breastfeeding.group <- relevel(cordata$Breastfeeding.group, ref = "group 2")


cordata <- cordata %>% mutate(roll2 = ifelse(Rolls.from.tummy.to.back <= roll, 
                                             1, 0))
cordata <- cordata %>% mutate(sit2 = ifelse(Sits.without.support <= sit, 1, 0))
cordata <- cordata %>% mutate(walk2 = ifelse(Walks.without.support <= walk,
                                             1, 0))
cordata <- cordata %>% mutate(climb2 = ifelse(Climbs.up.stairs.with.help <= climb,
                                              1, 0))

# Defining domain

cordata <- cordata %>% mutate(Gross.motor = ifelse(roll2 == 1 & sit2 == 1 & 
                                                     walk2 == 1 & climb2 == 1,
                                                   1, 0))


# Fitting the model
full.gross <- glm(Gross.motor~Breastfeeding.group+Age.at.delivery+
                    Maternal.educational.level+Wealth.status+
                    Maternal.employment+ Child.birth.weight +
                    Child.gender + Birth.order,
                  data = cordata, family = binomial)

summary(full.gross)

x <- exp(cbind(OR = coef(full.gross), confint(full.gross, level = 0.9)))

round(x, 3)

