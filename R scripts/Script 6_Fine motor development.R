library(tidyverse)
library(car)

# Loading dataset
cordata <- read_csv("Datasets/Data for fine motor development.csv")

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
grasp <- quantile(cordata$Reaches.and.grasps.objects, 0.75)
move <- quantile(cordata$Moves.things, 0.75)
pick <- quantile(cordata$Picks.things, 0.75)
stack <- quantile(cordata$Stacks.blocks, 0.75)

# categorizing variables
cordata <- cordata %>% mutate(Breastfeeding.group = ifelse(EBF.Time >=6 & 
                       Breastfeeding.time >= 24, "group 1","group 2"))

cordata$Breastfeeding.group <- as.factor(cordata$Breastfeeding.group)
cordata$Breastfeeding.group <- relevel(cordata$Breastfeeding.group, ref = "group 2")

table(cordata$Breastfeeding.group)
round(prop.table(table(cordata$Breastfeeding.group))*100, digits = 1)

cordata <- cordata %>% mutate(grasp2 = ifelse(Reaches.and.grasps.objects <= grasp,
                                              1, 0))
cordata <- cordata %>% mutate(move2 = ifelse(Moves.things <= move, 1, 0))
cordata <- cordata %>% mutate(pick2 = ifelse(Picks.things <= pick, 1, 0))
cordata <- cordata %>% mutate(stack2 = ifelse(Stacks.blocks <= stack,
                                              1, 0))

# Defining domain

cordata <- cordata %>% mutate(Fine.motor = ifelse(grasp2 == 1 & move2 == 1 & 
                                                    pick2 == 1 & stack2 == 1, 
                                                  1, 0))


# Fitting the model
full.fine <- glm(Fine.motor~Breastfeeding.group+Age.at.delivery+
                   Maternal.educational.level+Wealth.status+
                   Maternal.employment+ Child.birth.weight+
                   Child.gender + Birth.order,
                 data = cordata, family = binomial)

summary(full.fine)

x <- exp(cbind(OR = coef(full.fine), confint(full.fine, level = 0.90)))

round(x, 3)

