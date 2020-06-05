# We now consider the case with 3 test fires and one from the crime scene (suspect bullet)
# Therefore,
# For known matches, we have only 18 comparisons
# For known non-matches, 
# we have 15*3+30*3=135 comparisons(same bullet different lands + different bullet different lands)
library(tidyverse)
library(betareg)
library(bulletxtrctr)
# 1. we are going to first fit the mixture beta distribution for FAU 330===========================

# See the data (data not posted on github, use the ccf directly instead)
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

fau330_km_no_A_simulated <- simulate_bitamix(fau330_km_beta_no_A)

ggplot(data = fau330_km_no_A_simulated) +
  geom_line(aes(x = ccf, y = y), color = "green") +
  geom_line(aes(x = ccf, y = y1), color = "blue") +
  geom_line(aes(x = ccf, y = y2), color = "red")

fau330_knm_no_A_simulated <- simulate_bitamix(fau330_knm_beta_no_A)

ggplot(data = fau330_knm_no_A_simulated) +
  geom_line(aes(x = ccf, y = y), color = "green") +
  geom_line(aes(x = ccf, y = y1), color = "blue") +
  geom_line(aes(x = ccf, y = y2), color = "red")

# 3. See the ccfs for the leave out bullet with those known(used) bullets==========================
fau330_km_A <- fau330_km %>% filter(bullet1 == "A"|bullet2 == "A")
fau330_knm_A <- fau330_knm %>% filter(bullet1 == "A"|bullet2 == "A")

# We first estimate the distributions from the A involved comparison. 
# Then we could compare the distributions estimated to the previous ones(Should this use some of the EM algorithm result? Yes, but it seems to recommand to use simulation not theoretical results)
# Or should we use some other measurement to compare how close these distributions are, e.g. simulation based interval/test
# But note that we can known in advance in practice the km or knm.
fau330_km_bata_A <- betamix(ccf~1|1, data = fau330_km_A, k = 2)
fau330_knm_bata_A <- betamix(ccf~1|1, data = fau330_knm_A, k = 2)


# 4. Evaluate how well this is(use a simulation based test procedure)=============================================

# For test fires(known bullets), km

# a function to run the process(do simulations) with better adaption to the other cases

summary_simulated <- function(model, nrun, n_each_run){
  # don't use beta_params here
  param1 <-  cv_summary(list(model))[[1]]
  
  summary_test_km <- lapply(1:nrun, FUN = function(i) {
    simulated <- rbetamix(n = n_each_run, p1 = param1$prior[[1]], a1 = param1$params[[1]]$alpha, b1 = param1$params[[1]]$beta,
                          a2 = param1$params[[2]]$alpha, b2 = param1$params[[2]]$alpha)
    simulated <- data.frame(ccf = simulated)
    
    model <- betamix(ccf~1|1, data = simulated, k = 2)
    
    summarized <- cv_summary(list(model))[[1]]
    
    print(i)
    summarized
  })
  
  summary_test_km
}

# for 99 simulated samples, (not 100, because 5 and 95 are exact the quatile we want to choose)
summary_test_km <- summary_simulated(fau330_km_beta_no_A, nrun = 99, n_each_run = 18)

# For test fires(known bullets), knm
summary_test_knm <- summary_simulated(fau330_knm_beta_no_A, nrun = 99, n_each_run = 135)

# For suspect fire(known bullets), km
summary_test_km <- summary_simulated(fau330_km_bata_A, nrun = 99, n_each_run = 18)

# For suspect fire(known bullets), knm
summary_test_knm <- summary_simulated(fau330_knm_beta_A, nrun = 99, n_each_run = 105)