library(readxl)
library(writexl)
library(tidyverse)

Priors <- read_excel("Datasets/Prior data 1.xlsx")

colnames(Priors)

Priors <- Priors %>% mutate(
  beta.Random = round(log(HR.Random),3),
  variance.beta.Random = round(((log(UCI.Random) - log(LCI.Random))/(2*1.96))^2,4),
  sd.beta.Random = round(((log(UCI.Random) - log(LCI.Random))/(2*1.96)),4))

write_xlsx(Priors, "Datasets/Final prior data.xlsx")
