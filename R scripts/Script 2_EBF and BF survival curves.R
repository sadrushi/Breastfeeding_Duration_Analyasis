# Loading libraries
library(tidyverse)
library(survival)
library(survminer)
library(ggsci)
library(DescTools)
library(readxl)

# Importing dataset
data <- read_excel("Datasets/Breastfeeding_Duration_Data.xlsx")

# Creating Survival Object
attach(data)
surv.obj.1 <- Surv(EBF.time, Censoring.status.1)
surv.obj.2 <- Surv(Breastfeeding.time, Censoring.status.2)

# KP estimates
EBF.fit <- survfit(surv.obj.1 ~ 1, data = data)
BF.fit <- survfit(surv.obj.2 ~ 1, data = data)

# Statistics
print(EBF.fit, print.rmean = TRUE)
surv_summary(EBF.fit, data = data)

print(BF.fit, print.rmean = TRUE)
surv_summary(BF.fit, data = data)

# Median survival time
surv_median(EBF.fit)
surv_median(BF.fit)

# Survival Curve
ggsurvplot(EBF.fit, data = data, 
           conf.int = TRUE,
           conf.int.fill = "#F0E442",
           conf.int.alpha = 0.4,
           break.time.by = 1,
           legend = "none", 
           ggtheme = theme_bw(),
           xlab = "Duration of exclusive breastfeeding (months)",
           palette = "#0072B2",
           font.x = c(13, "black"),
           font.y = c(13, "black"),
           font.tickslab = c(11, "grey10")) 


ggsurvplot(BF.fit, data = data, 
           conf.int = TRUE,
           conf.int.fill = "#F0E442",
           conf.int.alpha = 0.4,
           break.time.by = 6,
           legend = "none", 
           ggtheme = theme_bw(),
           xlab = "Duration of extended breastfeeding (months)",
           palette = "#0072B2",
           censor.shape = 124,
           censor.size = 4,
           font.x = c(13, "black"),
           font.y = c(13, "black"),
           font.tickslab = c(11, "grey10"))

