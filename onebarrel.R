# We now consider the case with 3 test fires and one from the crime scene (suspect bullet)
# Therefore,
# For known matches, we have only 18 comparisons
# For known non-matches, 
# we have 15*3+30*3=135 comparisons(same bullet different lands + different bullet different lands)
library(tidyverse)
library(betareg)
library(bulletxtrctr)
# 1. we are going to first fit the mixture beta distribution for FAU 330===========================

# See the data
fau330_km <- land_score_km_generator2("./summer/data/")[[1]] 

fau330_km$bullet2
fau330_km$bullet1
fau330_km %>% filter(bullet2 != "A", bullet1 != "A")

fau330_knm <- land_score_knm_generator2("./summer/data/")[[1]]

fau330_knm$bullet2
fau330_knm$bullet1

fau330_knm %>% filter(bullet2 != "A", bullet1 != "A")

# fit the model with bullet A removed
fau330_km_no_A <- fau330_km %>% filter(bullet2 != "A", bullet1 != "A")
fau330_knm_no_A <- fau330_knm %>% filter(bullet2 != "A", bullet1 != "A")
## the knm 
fau330_knm_beta_no_A <- betamix(ccf~1|1, data = fau330_knm_no_A, k = 2)

fau330_knm_beta_no_A
fau330_knm_beta_no_A$flexmix@prior
fau330_knm_beta_no_A$flexmix@components

## the km
fau330_km_beta_no_A <- betamix(ccf~1|1, data = fau330_km_no_A, k = 2)

fau330_km_beta_no_A
fau330_km_beta_no_A$flexmix@prior
fau330_km_beta_no_A$flexmix@components

## See the plot
beta_params(fau330_km_beta_no_A)
beta_params(fau330_knm_beta_no_A)


#