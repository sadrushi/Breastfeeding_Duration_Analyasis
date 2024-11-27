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
         Child.birth.weight = factor(Child.birth.weight, levels = c("Low","Normal")))

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


################################################################################
###################### Classical Cox PH Model #########################
################################################################################

# Fitting classical Cox PH model
attach(Baysurv.data)
surv.obj <- Surv(BF.time, Censoring.status)

Clcox.fit <- coxph(surv.obj~Age.at.delivery+Maternal.educational.level+
                     Maternal.employment+Wealth.status+Delivery.mode+
                     Child.gender+Birth.order+Child.birth.weight+
                     Age.at.first.feeding)

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
informative_sds.R <- sqrt(Priors.data$variance.beta.Random+1)

round(informative_sds.R, 4)

informative_priors.R <- normal(location = informative_means.R, scale = informative_sds.R)

# ------------------------------------------------------------------------------

## Fitting Bayesian survival model

# Gompertz model

# Fit the model with prior predictive simulation
mod_gompertz1.prior <- stan_surv(
  formula = Surv(BF.time, Censoring.status) ~ Age.at.delivery+Maternal.educational.level+
    Maternal.employment+Wealth.status+Delivery.mode+
    Child.gender+Birth.order+Child.birth.weight+
    Age.at.first.feeding,
  data = Baysurv.data,
  basehaz = "gompertz",
  prior_PD = TRUE,
  prior = informative_priors.R,
  prior_intercept = normal(0,1),
  chains = 4,
  cores = 2,
  iter = 2*5000,
  seed = 345)

# Summarize prior results
summary(mod_gompertz1.prior)
tidy(mod_gompertz1.prior, conf.int = TRUE, conf.level = 0.95)

# Extract prior predictive samples
prior_samples <- as.data.frame(mod_gompertz1.prior)

# Fit the Cox proportional hazards model
mod_gompertz1 <- stan_surv(
  formula = Surv(BF.time, Censoring.status) ~ Age.at.delivery+Maternal.educational.level+
    Maternal.employment+Wealth.status+Delivery.mode+
    Child.gender+Birth.order+Child.birth.weight+
    Age.at.first.feeding,
  data = Baysurv.data,
  basehaz = "gompertz",
  prior = informative_priors.R,
  prior_intercept = normal(0,1),
  chains = 4,
  cores = 2,
  iter = 2*5000,
  seed = 345)

# Summarize posterior results
summary(mod_gompertz1)
tidy(mod_gompertz1, conf.int = TRUE, conf.level = 0.95)

# Summary of fitted model: with hazard ratios
tidy(mod_gompertz1, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95)

print(mod_gompertz1, digits = 3)

# Extract posterior samples
posterior_samples <- as.data.frame(mod_gompertz1)

# ------------------------------------------------------------------------------

# Comparison between prior and posterior

ggplot() +
  geom_density(data = prior_samples, aes(x = `Age.at.delivery18-30`, color = "Prior")) +
  geom_density(data = posterior_samples, aes(x = `Age.at.delivery18-30`, color = "Posterior")) +
  labs(title = "Comparison of Prior and Posterior for Maternal Age Effect",
       x = "Hazard Ratio",
       y = "Density",
       color = "Distribution") +
  theme_minimal()

ggplot() +
  geom_density(data = prior_samples, aes(x = `Maternal.educational.levelUp to ordinary level`, color = "Prior")) +
  geom_density(data = posterior_samples, aes(x = `Maternal.educational.levelUp to ordinary level`, color = "Posterior")) +
  labs(title = "Comparison of Prior and Posterior for Maternal Edu OL Effect",
       x = "Hazard Ratio",
       y = "Density",
       color = "Distribution") +
  theme_minimal()

ggplot() +
  geom_density(data = prior_samples, aes(x = `Maternal.educational.levelUp to advanced level`, color = "Prior")) +
  geom_density(data = posterior_samples, aes(x = `Maternal.educational.levelUp to advanced level`, color = "Posterior")) +
  labs(title = "Comparison of Prior and Posterior for Maternal Edu AL Effect",
       x = "Hazard Ratio",
       y = "Density",
       color = "Distribution") +
  theme_minimal()

ggplot() +
  geom_density(data = prior_samples, aes(x = `Maternal.employmentUnemployed`, color = "Prior")) +
  geom_density(data = posterior_samples, aes(x = `Maternal.employmentUnemployed`, color = "Posterior")) +
  labs(title = "Comparison of Prior and Posterior for Maternal Emp Unemp Effect",
       x = "Hazard Ratio",
       y = "Density",
       color = "Distribution") +
  theme_minimal()

ggplot() +
  geom_density(data = prior_samples, aes(x = `Wealth.statusLower`, color = "Prior")) +
  geom_density(data = posterior_samples, aes(x = `Wealth.statusLower`, color = "Posterior")) +
  labs(title = "Comparison of Prior and Posterior for Wealth lower Effect",
       x = "Hazard Ratio",
       y = "Density",
       color = "Distribution") +
  theme_minimal()

ggplot() +
  geom_density(data = prior_samples, aes(x = `Wealth.statusMiddle`, color = "Prior")) +
  geom_density(data = posterior_samples, aes(x = `Wealth.statusMiddle`, color = "Posterior")) +
  labs(title = "Comparison of Prior and Posterior for Wealth middle Effect",
       x = "Hazard Ratio",
       y = "Density",
       color = "Distribution") +
  theme_minimal()

ggplot() +
  geom_density(data = prior_samples, aes(x = `Delivery.modeNormal`, color = "Prior")) +
  geom_density(data = posterior_samples, aes(x = `Delivery.modeNormal`, color = "Posterior")) +
  labs(title = "Comparison of Prior and Posterior for Delivery normal Effect",
       x = "Hazard Ratio",
       y = "Density",
       color = "Distribution") +
  theme_minimal()

ggplot() +
  geom_density(data = prior_samples, aes(x = `Child.genderMale`, color = "Prior")) +
  geom_density(data = posterior_samples, aes(x = `Child.genderMale`, color = "Posterior")) +
  labs(title = "Comparison of Prior and Posterior for Child's gender male Effect",
       x = "Hazard Ratio",
       y = "Density",
       color = "Distribution") +
  theme_minimal()

ggplot() +
  geom_density(data = prior_samples, aes(x = `Birth.orderSecond`, color = "Prior")) +
  geom_density(data = posterior_samples, aes(x = `Birth.orderSecond`, color = "Posterior")) +
  labs(title = "Comparison of Prior and Posterior for Birth order second Effect",
       x = "Hazard Ratio",
       y = "Density",
       color = "Distribution") +
  theme_minimal()

ggplot() +
  geom_density(data = prior_samples, aes(x = `Birth.orderThird or above`, color = "Prior")) +
  geom_density(data = posterior_samples, aes(x = `Birth.orderThird or above`, color = "Posterior")) +
  labs(title = "Comparison of Prior and Posterior for Birth order >=3 Effect",
       x = "Hazard Ratio",
       y = "Density",
       color = "Distribution") +
  theme_minimal()

ggplot() +
  geom_density(data = prior_samples, aes(x = `Child.birth.weightNormal`, color = "Prior")) +
  geom_density(data = posterior_samples, aes(x = `Child.birth.weightNormal`, color = "Posterior")) +
  labs(title = "Comparison of Prior and Posterior for Birth weight normal Effect",
       x = "Hazard Ratio",
       y = "Density",
       color = "Distribution") +
  theme_minimal()

ggplot() +
  geom_density(data = prior_samples, aes(x = `Age.at.first.feeding< 6 months`, color = "Prior")) +
  geom_density(data = posterior_samples, aes(x = `Age.at.first.feeding< 6 months`, color = "Posterior")) +
  labs(title = "Comparison of Prior and Posterior for First feed <6 months Effect",
       x = "Hazard Ratio",
       y = "Density",
       color = "Distribution") +
  theme_minimal()
# ------------------------------------------------------------------------------

posterior_vs_prior(mod_gompertz1)

mcmc_intervals(mod_gompertz1) + vline_at(0)

plot_grid(
  bayesplot_grid(mcmc_intervals(mod_gompertz1.prior),
                 mcmc_intervals(mod_gompertz1),
                 titles = c("Prior", "Posterior"),
                 grid_args = list(nrow = 2)),
  
  bayesplot_grid(mcmc_hist(mod_gompertz1.prior),
                 mcmc_hist(mod_gompertz1),
                 titles = c("Prior", "Posterior"),
                 grid_args = list(nrow = 2)),
  
  ncol = 2
)

# ------------------------------------------------------------------------------

# MCMC diagnostics

# trace plots

mcmc_trace(mod_gompertz1, size = 0.1, 
           pars = c("(Intercept)", 
                    "Age.at.delivery18-30")) + 
  ggplot2::scale_color_discrete(type = c("#3B0F70FF", "#E7298A", "#E6AB02", 
                                         "#1B9E77")) +
  ggplot2::theme(axis.text = element_text(size = 12, color = "black"),
                 legend.title = element_text(size = 13.5),
                 strip.text = element_text(size = 14, color = "black"))


mcmc_trace(mod_gompertz1, size = 0.1, 
           pars = c("Maternal.educational.levelUp to ordinary level", 
                    "Maternal.educational.levelUp to advanced level")) +  
  ggplot2::scale_color_discrete(type = c("#3B0F70FF", "#E7298A", "#E6AB02", 
                                         "#1B9E77")) +
  ggplot2::theme(axis.text = element_text(size = 12, color = "black"),
                 legend.title = element_text(size = 13.5),
                 strip.text = element_text(size = 14, color = "black"))

mcmc_trace(mod_gompertz1, size = 0.1, 
           pars = c("Maternal.employmentUnemployed",
                    "Wealth.statusLower")) +  
  ggplot2::scale_color_discrete(type = c("#3B0F70FF", "#E7298A", 
                                         "#E6AB02", "#1B9E77")) +
  ggplot2::theme(axis.text = element_text(size = 12, color = "black"),
                 legend.title = element_text(size = 13.5),
                 strip.text = element_text(size = 14, color = "black"))



mcmc_trace(mod_gompertz1, size = 0.1, pars = c("Wealth.statusMiddle", 
                                              "Delivery.modeNormal")) + 
  ggplot2::scale_color_discrete(type = c("#3B0F70FF", "#E7298A", "#E6AB02", 
                                         "#1B9E77")) +
  ggplot2::theme(axis.text = element_text(size = 12, color = "black"),
                 legend.title = element_text(size = 13.5),
                 strip.text = element_text(size = 14, color = "black"))


mcmc_trace(mod_gompertz1, size = 0.1, pars = c("Child.genderMale", 
                                              "Birth.orderSecond")) + 
  ggplot2::scale_color_discrete(type = c("#3B0F70FF", "#E7298A", "#E6AB02", 
                                         "#1B9E77")) +
  ggplot2::theme(axis.text = element_text(size = 12, color = "black"),
                 legend.title = element_text(size = 13.5),
                 strip.text = element_text(size = 14, color = "black"))


mcmc_trace(mod_gompertz1, size = 0.1, pars = c("Birth.orderThird or above", 
                                              "Child.birth.weightNormal")) +
  ggplot2::scale_color_discrete(type = c("#3B0F70FF", "#E7298A", "#E6AB02", 
                                         "#1B9E77")) +
  ggplot2::theme(axis.text = element_text(size = 12, color = "black"),
                 legend.title = element_text(size = 13.5),
                 strip.text = element_text(size = 14, color = "black"))


mcmc_trace(mod_gompertz1, size = 0.1, pars = c("Age.at.first.feeding< 6 months")) +
  ggplot2::scale_color_discrete(type = c("#3B0F70FF", "#E7298A", "#E6AB02", 
                                         "#1B9E77")) +
  ggplot2::theme(axis.text = element_text(size = 12, color = "black"),
                 legend.title = element_text(size = 13.5),
                 strip.text = element_text(size = 14, color = "black"))



# density plots

mcmc_dens_overlay(mod_gompertz1,  
                  pars = c("(Intercept)", 
                           "Age.at.delivery18-30")) + 
  ggplot2::scale_color_discrete(type = c("#3B0F70FF", "#E7298A", "#E6AB02", 
                                         "#1B9E77")) +
  ggplot2::theme(axis.text = element_text(size = 12, color = "black"),
                 legend.title = element_text(size = 13.5),
                 strip.text = element_text(size = 14, color = "black"))

mcmc_dens_overlay(mod_gompertz1,  
                  pars = c("Maternal.educational.levelUp to ordinary level", 
                           "Maternal.educational.levelUp to advanced level")) + 
  ggplot2::scale_color_discrete(type = c("#3B0F70FF", "#E7298A", "#E6AB02", 
                                         "#1B9E77")) +
  ggplot2::theme(axis.text = element_text(size = 12, color = "black"),
                 legend.title = element_text(size = 13.5),
                 strip.text = element_text(size = 14, color = "black"))


mcmc_dens_overlay(mod_gompertz1, 
                  pars = c("Maternal.employmentUnemployed",
                           "Wealth.statusLower")) +  
  ggplot2::scale_color_discrete(type = c("#3B0F70FF", "#E7298A", "#E6AB02", 
                                         "#1B9E77")) +
  ggplot2::theme(axis.text = element_text(size = 12, color = "black"),
                 legend.title = element_text(size = 13.5),
                 strip.text = element_text(size = 14, color = "black"))

mcmc_dens_overlay(mod_gompertz1, 
                  pars = c("Wealth.statusMiddle", 
                           "Delivery.modeNormal")) +  
  ggplot2::scale_color_discrete(type = c("#3B0F70FF", "#E7298A", "#E6AB02", 
                                         "#1B9E77")) +
  ggplot2::theme(axis.text = element_text(size = 12, color = "black"),
                 legend.title = element_text(size = 13.5),
                 strip.text = element_text(size = 14, color = "black"))

mcmc_dens_overlay(mod_gompertz1, 
                  pars = c("Child.genderMale", 
                           "Birth.orderSecond")) +  
  ggplot2::scale_color_discrete(type = c("#3B0F70FF", "#E7298A", "#E6AB02", 
                                         "#1B9E77")) +
  ggplot2::theme(axis.text = element_text(size = 12, color = "black"),
                 legend.title = element_text(size = 13.5),
                 strip.text = element_text(size = 14, color = "black"))


mcmc_dens_overlay(mod_gompertz1, 
                  pars = c("Birth.orderThird or above", 
                           "Child.birth.weightNormal"
                  )) +   ggplot2::scale_color_discrete(type = c("#3B0F70FF", "#E7298A",
                                                                "#E6AB02", "#1B9E77")) +
  ggplot2::theme(axis.text = element_text(size = 12, color = "black"),
                 legend.title = element_text(size = 13.5),
                 strip.text = element_text(size = 14, color = "black"))



mcmc_dens_overlay(mod_gompertz1, pars = c("Age.at.first.feeding< 6 months")) + 
  ggplot2::scale_color_discrete(type = c("#3B0F70FF", "#E7298A", "#E6AB02", 
                                         "#1B9E77")) +
  ggplot2::theme(axis.text = element_text(size = 12, color = "black"),
                 legend.title = element_text(size = 13.5),
                 strip.text = element_text(size = 14, color = "black"))


# R-hat
color_scheme_set(c("#0072B2", "#400040", "#bf7fbf", "#bf7fbf",  "#0072B2",
                   "#400040"))

r <- rhat(mod_gompertz1)
mcmc_rhat_hist(r) +
  ggplot2::theme(axis.text = element_text(size = 12, color = "black"),
                 legend.title = element_text(size = 12),
                 strip.text = element_text(size = 12, color = "black"))


# Effective sample size
color_scheme_set(c("#0072B2", "#400040", "#bf7fbf", "#bf7fbf",  "#0072B2",
                   "#400040"))

n <- neff_ratio(mod_gompertz1)
mcmc_neff(n, size = 2) + yaxis_text(hjust = 1) + 
  ggplot2::theme(axis.text = element_text(size = 12, color = "black"),
                 legend.title = element_text(size = 12))

