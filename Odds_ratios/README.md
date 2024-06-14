library(tidyverse)
library("geepack")

setwd("~/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/NG_transmission/Odds_ratios")

Metadata = read.table("../Metadata/odds.csv", header = TRUE, sep= ",")

fit <- geeglm(formula = Status ~ Sex + AgeGroupSum + countif + PEN + TET +
                CTRIX + CIPRO + AZITH + Site_summary
               data = Metadata, 
               id = cluster, 
               family = binomial, 
               corstr = "independence")
fit
summary(fit)

broom::tidy(x = fit, exp=T, conf.int = TRUE)