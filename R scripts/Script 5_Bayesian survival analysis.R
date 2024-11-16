library(tidyverse)
library(survival)
library(survminer)
library(rstanarm)
library(rstan)
library(bayesplot)
library(cowplot)
library(pec)
library(loo)
library(broom.mixed)
library(gridExtra)
library(readxl)
library(writexl)

# Loading datasets
Baysurv.data <- read_excel("Datasets/Data for BF duration analysis.xlsx")
Priors.data <- Final_prior_data <- read_excel("Datasets/Final prior data.xlsx")


# Factorizing variables
Baysurv.data <- Baysurv.data %>%
  mutate(Age.at.delivery = factor(Age.at.delivery, levels = c("18-30","31-45")),
         Maternal.educational.level = factor(Maternal.educational.level, levels = c("Up to ordinary level","Up to advanced level","Higher education")),
         Maternal.employment = factor(Maternal.employment, levels = c("Unemployed","Employed")),
         Delivery.mode = factor(Delivery.mode, levels = c("Normal","Cesarean")),
         Wealth.status = factor(Wealth.status, levels = c("Lower","Middle","Higher")),
         Child.gender = factor(Child.gender, levels = c("Female","Male")),
         Birth.order = factor(Birth.order, levels = c("First","Second","Third or above")),
         Age.at.first.feeding = factor(Age.at.first.feeding, levels = c("< 6 months",">= 6 months")),
         Child.birth.weight = factor(Child.birth.weight, levels = c("Low","Normal")),
         Parity = factor(Parity, levels = c("1-2 children","3-4 children")))

# Defining reference levels
Baysurv.data$Age.at.delivery <- relevel(Baysurv.data$Age.at.delivery, ref = "31-45")
Baysurv.data$Maternal.educational.level <- relevel(Baysurv.data$Maternal.educational.level,ref = "Higher education")
Baysurv.data$Maternal.employment <- relevel(Baysurv.data$Maternal.employment,ref = "Employed")
Baysurv.data$Delivery.mode <- relevel(Baysurv.data$Delivery.mode,ref = "Cesarean")
Baysurv.data$Wealth.status <- relevel(Baysurv.data$Wealth.status,ref = "Higher")
Baysurv.data$Child.gender <- relevel(Baysurv.data$Child.gender,ref = "Female")
Baysurv.data$Birth.order <- relevel(Baysurv.data$Birth.order,ref = "First")
Baysurv.data$Age.at.first.feeding <- relevel(Baysurv.data$Age.at.first.feeding,ref = ">= 6 months")
Baysurv.data$Child.birth.weight <- relevel(Baysurv.data$Child.birth.weight, ref = "Low")
Baysurv.data$Parity <- relevel(Baysurv.data$Parity, ref = "1-2 children")

################################################################################
###################### Classical Cox PH Model #########################
################################################################################

# Fitting classical Cox PH model
attach(Baysurv.data)
surv.obj <- Surv(BF.time, Censoring.status)

Clcox.fit <- coxph(surv.obj~Age.at.delivery+Maternal.educational.level+
                     Maternal.employment+Wealth.status+Delivery.mode+
                     Child.gender+Birth.order+Child.birth.weight+
                     Age.at.first.feeding+Parity)

# Summary of fitted model
summary(Clcox.fit)

# Testing proportional Hazards assumption
cox.zph(Clcox.fit)

################################################################################
###################### Bayesian Cox PH Model #########################
################################################################################

## Finding appropriate baseline hazard

# Exponential model
mod_exp <- stan_surv(
  formula = Surv(BF.time, Censoring.status) ~ 1,
  data = Baysurv.data,
  basehaz = "exp",
  prior_intercept = normal(0,1),
  chains = 4,
  iter = 2*5000,
  seed = 345)

# Weibull model
mod_weibull <- stan_surv(
  formula = Surv(BF.time, Censoring.status) ~ 1,
  data = Baysurv.data,
  basehaz = "weibull",
  prior_intercept = normal(0,1),
  chains = 4,
  iter = 2*5000,
  seed = 345)

# Gompertz model
mod_gompertz <- stan_surv(
  formula = Surv(BF.time, Censoring.status) ~ 1,
  data = Baysurv.data,
  basehaz = "gompertz",
  prior_intercept = normal(0,1),
  chains = 4,
  iter = 2*5000,
  seed = 345)

# plot of the estimated baseline hazard function
fits_stan <- list("Exponential" = mod_exp,
                  "Weibull" = mod_weibull,
                  "Gompertz" = mod_gompertz)

plot.base <- map(fits_stan, ~ plot(., plotfun = "basehaz", 
                                   limits = "none") )

bayesplot_grid(
  plots = plot.base,
  titles = names(fits_stan),
  grid_args = list(ncol = 3))

# Compare models with different baseline with survival curve
plots <- map(fits_stan, ~ ps_check(.))

bayesplot_grid(
  plots = plots,
  titles = names(fits_stan),
  grid_args = list(ncol = 3))

# Compare models with different baseline
loo_1 <- loo(mod_exp)
loo_2 <- loo(mod_weibull)
loo_3 <- loo(mod_gompertz)

loo_1$estimates
loo_2$estimates
loo_3$estimates

loo_compare(loo_1, loo_2, loo_3)

# ------------------------------------------------------------------------------

## Setting priors

# Setting priors from the Random-effects model

informative_means.R <- Priors.data$beta.Random
informative_sds.R <- Priors.data$sd.beta.Random

informative_priors.R <- normal(location = informative_means.R, scale = informative_sds.R)

# ------------------------------------------------------------------------------

## Fitting Bayesian survival model

# Gompertz model

mod_gompertz_final <- stan_surv(
  formula = Surv(BF.time, Censoring.status) ~ Age.at.delivery+Maternal.educational.level+
    Maternal.employment+Wealth.status+Delivery.mode+
    Child.gender+Birth.order+Child.birth.weight+
    Age.at.first.feeding+Parity,
  data = Baysurv.data,
  basehaz = "gompertz",
  prior = informative_priors.R,
  prior_intercept = normal(0,1),
  chains = 4,
  cores = 2,
  iter = 2*5000,
  seed = 345)

# Prior summary
prior_summary(mod_gompertz_final)

# Summary of fitted model: with estimates
BaySurv.outcome.estimate <- tidy(mod_gompertz_final, conf.int = TRUE, conf.level = 0.95)

# Summary of fitted model: with hazard ratios
BaySurv.outcome.HR <- tidy(mod_gompertz_final, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95)

#BaySurv.outcome <- merge(BaySurv.outcome.estimate, BaySurv.outcome.HR, by = "term")
#write_xlsx(BaySurv.outcome, "Datasets/Bayesian model outcome final.xlsx")

print(mod_gompertz_final, digits = 3)


# ------------------------------------------------------------------------------

# MCMC diagnostics

# trace plots

mcmc_trace(mod_gompertz_final, size = 0.1, 
           pars = c("(Intercept)", 
                    "Age.at.delivery18-30")) + 
  ggplot2::scale_color_discrete(type = c("#3B0F70FF", "#E7298A", "#E6AB02", 
                                         "#1B9E77")) +
  ggplot2::theme(axis.text = element_text(size = 12, color = "black"),
                 legend.title = element_text(size = 13.5),
                 strip.text = element_text(size = 14, color = "black"))


mcmc_trace(mod_gompertz_final, size = 0.1, 
           pars = c("Maternal.educational.levelUp to ordinary level", 
                    "Maternal.educational.levelUp to advanced level")) +  
  ggplot2::scale_color_discrete(type = c("#3B0F70FF", "#E7298A", "#E6AB02", 
                                         "#1B9E77")) +
  ggplot2::theme(axis.text = element_text(size = 12, color = "black"),
                 legend.title = element_text(size = 13.5),
                 strip.text = element_text(size = 14, color = "black"))

mcmc_trace(mod_gompertz_final, size = 0.1, 
           pars = c("Maternal.employmentUnemployed",
                    "Wealth.statusLower")) +  
  ggplot2::scale_color_discrete(type = c("#3B0F70FF", "#E7298A", 
                                         "#E6AB02", "#1B9E77")) +
  ggplot2::theme(axis.text = element_text(size = 12, color = "black"),
                 legend.title = element_text(size = 13.5),
                 strip.text = element_text(size = 14, color = "black"))



mcmc_trace(mod_gompertz_final, size = 0.1, pars = c("Wealth.statusMiddle", 
                                              "Delivery.modeNormal")) + 
  ggplot2::scale_color_discrete(type = c("#3B0F70FF", "#E7298A", "#E6AB02", 
                                         "#1B9E77")) +
  ggplot2::theme(axis.text = element_text(size = 12, color = "black"),
                 legend.title = element_text(size = 13.5),
                 strip.text = element_text(size = 14, color = "black"))


mcmc_trace(mod_gompertz_final, size = 0.1, pars = c("Child.genderMale", 
                                              "Birth.orderSecond")) + 
  ggplot2::scale_color_discrete(type = c("#3B0F70FF", "#E7298A", "#E6AB02", 
                                         "#1B9E77")) +
  ggplot2::theme(axis.text = element_text(size = 12, color = "black"),
                 legend.title = element_text(size = 13.5),
                 strip.text = element_text(size = 14, color = "black"))


mcmc_trace(mod_gompertz_final, size = 0.1, pars = c("Birth.orderThird or above", 
                                              "Child.birth.weightNormal")) +
  ggplot2::scale_color_discrete(type = c("#3B0F70FF", "#E7298A", "#E6AB02", 
                                         "#1B9E77")) +
  ggplot2::theme(axis.text = element_text(size = 12, color = "black"),
                 legend.title = element_text(size = 13.5),
                 strip.text = element_text(size = 14, color = "black"))


mcmc_trace(mod_gompertz_final, size = 0.1, pars = c("Age.at.first.feeding< 6 months", 
                                              "Parity3-4 children")) +
  ggplot2::scale_color_discrete(type = c("#3B0F70FF", "#E7298A", "#E6AB02", 
                                         "#1B9E77")) +
  ggplot2::theme(axis.text = element_text(size = 12, color = "black"),
                 legend.title = element_text(size = 13.5),
                 strip.text = element_text(size = 14, color = "black"))



# density plots

mcmc_dens_overlay(mod_gompertz_final,  
                  pars = c("(Intercept)", 
                           "Age.at.delivery18-30")) + 
  ggplot2::scale_color_discrete(type = c("#3B0F70FF", "#E7298A", "#E6AB02", 
                                         "#1B9E77")) +
  ggplot2::theme(axis.text = element_text(size = 12, color = "black"),
                 legend.title = element_text(size = 13.5),
                 strip.text = element_text(size = 14, color = "black"))

mcmc_dens_overlay(mod_gompertz_final,  
                  pars = c("Maternal.educational.levelUp to ordinary level", 
                           "Maternal.educational.levelUp to advanced level")) + 
  ggplot2::scale_color_discrete(type = c("#3B0F70FF", "#E7298A", "#E6AB02", 
                                         "#1B9E77")) +
  ggplot2::theme(axis.text = element_text(size = 12, color = "black"),
                 legend.title = element_text(size = 13.5),
                 strip.text = element_text(size = 14, color = "black"))


mcmc_dens_overlay(mod_gompertz_final, 
                  pars = c("Maternal.employmentUnemployed",
                           "Wealth.statusLower")) +  
  ggplot2::scale_color_discrete(type = c("#3B0F70FF", "#E7298A", "#E6AB02", 
                                         "#1B9E77")) +
  ggplot2::theme(axis.text = element_text(size = 12, color = "black"),
                 legend.title = element_text(size = 13.5),
                 strip.text = element_text(size = 14, color = "black"))

mcmc_dens_overlay(mod_gompertz_final, 
                  pars = c("Wealth.statusMiddle", 
                           "Delivery.modeNormal")) +  
  ggplot2::scale_color_discrete(type = c("#3B0F70FF", "#E7298A", "#E6AB02", 
                                         "#1B9E77")) +
  ggplot2::theme(axis.text = element_text(size = 12, color = "black"),
                 legend.title = element_text(size = 13.5),
                 strip.text = element_text(size = 14, color = "black"))

mcmc_dens_overlay(mod_gompertz_final, 
                  pars = c("Child.genderMale", 
                           "Birth.orderSecond")) +  
  ggplot2::scale_color_discrete(type = c("#3B0F70FF", "#E7298A", "#E6AB02", 
                                         "#1B9E77")) +
  ggplot2::theme(axis.text = element_text(size = 12, color = "black"),
                 legend.title = element_text(size = 13.5),
                 strip.text = element_text(size = 14, color = "black"))


mcmc_dens_overlay(mod_gompertz_final, 
                  pars = c("Birth.orderThird or above", 
                           "Child.birth.weightNormal"
                  )) +   ggplot2::scale_color_discrete(type = c("#3B0F70FF", "#E7298A",
                                                                "#E6AB02", "#1B9E77")) +
  ggplot2::theme(axis.text = element_text(size = 12, color = "black"),
                 legend.title = element_text(size = 13.5),
                 strip.text = element_text(size = 14, color = "black"))



mcmc_dens_overlay(mod_gompertz_final, pars = c("Age.at.first.feeding< 6 months", 
                                         "Parity3-4 children")) + 
  ggplot2::scale_color_discrete(type = c("#3B0F70FF", "#E7298A", "#E6AB02", 
                                         "#1B9E77")) +
  ggplot2::theme(axis.text = element_text(size = 12, color = "black"),
                 legend.title = element_text(size = 13.5),
                 strip.text = element_text(size = 14, color = "black"))


# R-hat
color_scheme_set(c("#800080", "#400040", "#bf7fbf", "#bf7fbf",  "#800080",
                   "#400040"))

r <- rhat(mod_gompertz_final)
mcmc_rhat_hist(r) +
  ggplot2::theme(axis.text = element_text(size = 12, color = "black"),
                 legend.title = element_text(size = 12),
                 strip.text = element_text(size = 12, color = "black"))


# Effective sample size
color_scheme_set(c("#800080", "#400040", "#bf7fbf", "#bf7fbf",  "#800080",
                   "#400040"))

n <- neff_ratio(mod_gompertz_final)
mcmc_neff(n, size = 2) + yaxis_text(hjust = 1) + 
  ggplot2::theme(axis.text = element_text(size = 12, color = "black"),
                 legend.title = element_text(size = 12))

