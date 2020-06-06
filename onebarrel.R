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
# But note that we can't know in advance in practice the km or knm.
fau330_km_beta_A <- betamix(ccf~1|1, data = fau330_km_A, k = 2)
fau330_knm_beta_A <- betamix(ccf~1|1, data = fau330_knm_A, k = 2)


# 4. Evaluate how well this is(use a simulation based test procedure)=============================================

# For test fires(known bullets), km

# a function to run the process(do simulations) with better adaption to the other cases

summary_simulated <- function(model, nrun, n_each_run){
  # don't use beta_params here
  param1 <-  cv_summary(list(model))[[1]]
  
  work_flow <- function(i) {
    simulated <- rbetamix(n = n_each_run, p1 = param1$prior[[1]], a1 = param1$params[[1]]$alpha, b1 = param1$params[[1]]$beta,
                          a2 = param1$params[[2]]$alpha, b2 = param1$params[[2]]$alpha)
    simulated <- data.frame(ccf = simulated)
    
    model <- betamix(ccf~1|1, data = simulated, k = 2)
    
    summarized <- cv_summary(list(model))[[1]]
    
    print(i)
    summarized
  }
  
  summary_test_km <- lapply(1:nrun, FUN = function(i) {try(work_flow(i))})
  
  summary_test_km
}
# This sometimes fail because of the EM algorithm failed to find a starting value for certain data sets
# (sometimes it reports only one component? Sometimes no converges?(if no converges, there is an error/warning(two cases), and will be shown in the following process))

set.seed(20200606)
# for 99 simulated samples, (not 100, because 5 and 95 are exact the quatile we want to choose)
summary_test_km <- summary_simulated(fau330_km_beta_no_A, nrun = 99, n_each_run = 18)

# For test fires(known bullets), knm
summary_test_knm <- summary_simulated(fau330_knm_beta_no_A, nrun = 99, n_each_run = 135)

# For suspect fire(known bullets), km
summary_suspect_km <- summary_simulated(fau330_km_beta_A, nrun = 99, n_each_run = 18)

# For suspect fire(known bullets), knm
summary_suspect_knm <- summary_simulated(fau330_knm_beta_A, nrun = 99, n_each_run = 105)


# Since there is one sample in summary_suspect_km, failed to get starting values and reach a estimation, we might add another one to replace it
# This failed because no proper starting values, but report as error "no convergence", when there is good starting values but no convergence, report a warning
summary_suspect_km_66_fail <- summary_suspect_km[[66]]

summary_suspect_km <- summary_simulated(fau330_km_beta_A, nrun = 1, n_each_run = 18)

# Now it's ready to investigate

# To me, intuitively
# for beta distribution, the absolute values(variance related) of alpla and beta are not good indicator of stablity, but their ratios are more important

# To extract: alpha, beta, mu, phi, prior for both components (note, the direct result for params of the model are not mu, phi, but with a link function)

extract_params <- function(summary_something) {
  
  temp1 <- lapply(summary_something, FUN = function(d){
    p1 <- d$prior[[1]]
    c1_alpha <- d$params[[1]]$alpha
    c1_beta <- d$params[[1]]$beta
    c1_mu <- plogis(d$mu_phi_params[[1]]$mean)
    c1_phi <- exp(d$mu_phi_params[[1]]$precision)
    
    c2_alpha <- d$params[[2]]$alpha
    c2_beta <- d$params[[2]]$beta
    c2_mu <- plogis(d$mu_phi_params[[2]]$mean)
    c2_phi <- exp(d$mu_phi_params[[2]]$precision)
    
    output <- c(p1, c1_alpha, c1_beta, c1_mu, c2_alpha, c2_beta, c2_mu, c2_phi)
    output
  })
  
  temp2 <- lapply(1:8, FUN = function(i){
    
    theparamvector <- lapply(temp1, FUN = function(d2){
      output <- d2[i]
    }) %>% unlist()
    
    theparamvector
  }) 
  
  names(temp2) <-  c("p1", "c1_alpha", "c1_beta", "c1_mu", "c2_alpha", "c2_beta", "c2_mu", "c2_phi")
  
  temp3 <- as.data.frame(temp2)
  temp3
}

summary_test_km_data <- extract_params(summary_test_km)
summary_test_knm_data <- extract_params(summary_test_knm)
summary_suspect_km_data <- extract_params(summary_suspect_km)
summary_suspect_knm_data <- extract_params(summary_suspect_knm)

# The 95% confidence interval for aplha (expect to be extremely wide)

# The 95% confidence interval for beta (expect to be extremely wide)

# The 95% confidence interval for mu (expect to be extremely wide)

# The 95% confidence interval for phi (expect to be extremely wide)

# 5. Do similar estimations to more single barrel cases (within or out of this FAU330 barrel)==========

