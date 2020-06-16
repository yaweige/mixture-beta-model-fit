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

summary_suspect_km[[66]] <- summary_simulated(fau330_km_beta_A, nrun = 1, n_each_run = 18)

# saveRDS(summary_test_km, file = "summary_test_km.rds")
# saveRDS(summary_test_knm, file = "summary_test_knm.rds")
# saveRDS(summary_suspect_km, file = "summary_suspect_km.rds")
# saveRDS(summary_suspect_knm, file = "summary_suspect_knm.rds")
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

# The 95% Bootstrap condifence interval
# lower: 2*estimated-bootstrap[(M+1)*(1-alpha/2)]
# upper: 2*estimated-bootstrap[(M+1)*(alpha/2)]


bootstrap_intervals <- function(model, summary_x_x_data, alpha = 0.1){
  param1 <-  cv_summary(list(model))[[1]]
  p1 <- param1$prior[[1]]
  c1_alpha <- param1$params[[1]]$alpha
  c1_beta <- param1$params[[1]]$beta
  c1_mu <- plogis(param1$mu_phi_params[[1]]$mean)
  c1_phi <- exp(param1$mu_phi_params[[1]]$precision)
  
  c2_alpha <- param1$params[[2]]$alpha
  c2_beta <- param1$params[[2]]$beta
  c2_mu <- plogis(param1$mu_phi_params[[2]]$mean)
  c2_phi <- exp(param1$mu_phi_params[[2]]$precision)
  
  estimated <- c(p1, c1_alpha, c1_beta, c1_mu, c2_alpha, c2_beta, c2_mu, c2_phi)
  
  quantile_position_large <- floor((nrow(summary_x_x_data)+1)*(1-alpha/2))
  quantile_position_small <- floor((nrow(summary_x_x_data)+1)*(alpha/2))
  
  bootstrap_sample_small <- lapply(summary_x_x_data, FUN = function(x){
    output <- x[order(x)]
    output[quantile_position_small]
  }) %>% unlist()
  
  bootstrap_sample_large <- lapply(summary_x_x_data, FUN = function(x){
    output <- x[order(x)]
    output[quantile_position_large]
  }) %>% unlist() 
  
  lower <- 2*estimated - bootstrap_sample_large
  upper <- 2*estimated - bootstrap_sample_small
  
  output <- matrix(c(lower, upper), byrow = T, nrow = 2)
  colnames(output) <- c("p1", "c1_alpha", "c1_beta", "c1_mu", "c2_alpha", "c2_beta", "c2_mu", "c2_phi")
  rownames(output) <- paste0((1-alpha)*100, c("%lower", "%upper"))
  output
}

bootstrap_intervals_diff <- function(model1, model2, sampled_data1, sampled_data2, alpha = 0.1){
  param1 <-  cv_summary(list(model1))[[1]]
  m1_p1 <- param1$prior[[1]]
  m1_c1_alpha <- param1$params[[1]]$alpha
  m1_c1_beta <- param1$params[[1]]$beta
  m1_c1_mu <- plogis(param1$mu_phi_params[[1]]$mean)
  m1_c1_phi <- exp(param1$mu_phi_params[[1]]$precision)
  
  m1_c2_alpha <- param1$params[[2]]$alpha
  m1_c2_beta <- param1$params[[2]]$beta
  m1_c2_mu <- plogis(param1$mu_phi_params[[2]]$mean)
  m1_c2_phi <- exp(param1$mu_phi_params[[2]]$precision)
  
  m1_estimated <- c(m1_p1, m1_c1_alpha, m1_c1_beta, m1_c1_mu, m1_c2_alpha, m1_c2_beta, m1_c2_mu, m1_c2_phi)
  
  param2 <-  cv_summary(list(model2))[[1]]
  m2_p1 <- param2$prior[[1]]
  m2_c1_alpha <- param2$params[[1]]$alpha
  m2_c1_beta <- param2$params[[1]]$beta
  m2_c1_mu <- plogis(param2$mu_phi_params[[1]]$mean)
  m2_c1_phi <- exp(param2$mu_phi_params[[1]]$precision)
  
  m2_c2_alpha <- param2$params[[2]]$alpha
  m2_c2_beta <- param2$params[[2]]$beta
  m2_c2_mu <- plogis(param2$mu_phi_params[[2]]$mean)
  m2_c2_phi <- exp(param2$mu_phi_params[[2]]$precision)
  
  m2_estimated <- c(m2_p1, m2_c1_alpha, m2_c1_beta, m2_c1_mu, m2_c2_alpha, m2_c2_beta, m2_c2_mu, m2_c2_phi)
  
  difference <- m1_estimated-m2_estimated
  
  if (nrow(sampled_data1) >= nrow(sampled_data2)) {
    nrow_used <- nrow(sampled_data2)} else {
      nrow_used <- nrow(sampled_data1)}
  
  if (ncol(sampled_data1) != ncol(sampled_data2)) stop("data has different number of columns")
  
  sampled_diff <- lapply(1:ncol(sampled_data1), FUN = function(i) {
    output <- sampled_data1[1:nrow_used, i] - sampled_data2[1:nrow_used, i]
    output
  }) %>% as.data.frame()
  
  quantile_position_large <- floor((nrow_used+1)*(1-alpha/2))
  quantile_position_small <- floor((nrow_used+1)*(alpha/2))
  
  bootstrap_sample_small <- lapply(sampled_diff, FUN = function(x){
    output <- x[order(x)]
    output[quantile_position_small]
  }) %>% unlist()
  
  bootstrap_sample_large <- lapply(sampled_diff, FUN = function(x){
    output <- x[order(x)]
    output[quantile_position_large]
  }) %>% unlist() 
  
  
  
  lower <- 2*difference - bootstrap_sample_large
  upper <- 2*difference - bootstrap_sample_small
  
  output <- matrix(c(lower, upper), byrow = T, nrow = 2)
  colnames(output) <- c("p1_diff", "c1_alpha_diff", "c1_beta_diff", "c1_mu_diff", "c2_alpha_diff", "c2_beta_diff", "c2_mu_diff", "c2_phi_diff")
  rownames(output) <- paste0((1-alpha)*100, c("%lower", "%upper"))
  output
}


# The four models we estimated, bootstrap CI

aa <- bootstrap_intervals(fau330_km_beta_no_A, summary_test_km_data)
bb <- bootstrap_intervals(fau330_knm_beta_no_A, summary_test_knm_data)
cc <- bootstrap_intervals(fau330_km_beta_A, summary_suspect_km_data)
dd <- bootstrap_intervals(fau330_knm_beta_A, summary_suspect_knm_data)

bootstrap90_CI <- list(No_A_km = aa, NO_A_knm = bb, A_km = cc, A_knm = dd)

# saveRDS(bootstrap90_CI, file = "bootstrap90_CI.rds")

# How to test if two parameters are same or not, should we do a bootstrap interval for that? How?
# Inference conclusions:
# 1. Some of the intervals cover beyond the parameter space. (Is there anything wrong? This could be the case if the parameter varies a lot)
# 2. Should have a closer look of the bootstrap samples, e.g. boxplot
# 3. How to interpret the interval and test results if it goes beyond the parameter space?(should we cut it or how to report this, or shoud we have larger simulated samples?)
# 4. Obviously, the mean values are well seperated between KM and KNM
# 5. p1 are also well seperated except A_km

# we should have a look at the intervals of differences
# How to construct the difference? If we do every possible pairwise difference, is that close to form a U statistic?
# Or we do for the corresponding pairs in the simulated list (Do this now)

# bootstrap intervals for difference
ee <- bootstrap_intervals_diff(fau330_knm_beta_no_A, fau330_knm_beta_A, 
                               summary_test_knm_data, summary_suspect_knm_data,
                               alpha = 0.1)
ff <- bootstrap_intervals_diff(fau330_km_beta_no_A, fau330_km_beta_A, 
                               summary_test_km_data, summary_suspect_km_data,
                               alpha = 0.1)

bootstrap90_CI_diff <- list(km_noA_minus_A = ff, knm_noA_minus_A = ee)
#saveRDS(bootstrap90_CI_diff, file = "bootstrap90_CI_diff.rds")
# 5. Redo the above with lager data set run on the server========================================

larger_summary_test_km <- readRDS("./upload/summary_test_km.rds")
larger_summary_test_knm <- readRDS("./upload/summary_test_knm.rds")
larger_summary_suspect_km <- readRDS("./upload/summary_suspect_km.rds")
larger_summary_suspect_knm <- readRDS("./upload/summary_suspect_knm.rds")


# Detect the failed cases
larger_summary_test_km_failed <- lapply(larger_summary_test_km, FUN = is.list) %>% unlist()
sum(larger_summary_test_km_failed)

larger_summary_test_knm_failed <- lapply(larger_summary_test_knm, FUN = is.list) %>% unlist()
sum(larger_summary_test_knm_failed)

larger_summary_suspect_km_failed <- lapply(larger_summary_suspect_km, FUN = is.list) %>% unlist()
sum(larger_summary_test_km_failed)

larger_summary_suspect_knm_failed <- lapply(larger_summary_suspect_knm, FUN = is.list) %>% unlist()
sum(larger_summary_test_knm_failed)

larger_summary_test_km_data <- extract_params(larger_summary_test_km[larger_summary_test_km_failed])
larger_summary_test_knm_data <- extract_params(larger_summary_test_knm[larger_summary_test_knm_failed])
larger_summary_suspect_km_data <- extract_params(larger_summary_suspect_km[larger_summary_suspect_km_failed])
larger_summary_suspect_knm_data <- extract_params(larger_summary_suspect_knm[larger_summary_suspect_knm_failed])
gg <- bootstrap_intervals(fau330_km_beta_no_A, larger_summary_test_km_data)
hh <- bootstrap_intervals(fau330_knm_beta_no_A, larger_summary_test_knm_data)
ii <- bootstrap_intervals(fau330_km_beta_A, larger_summary_suspect_km_data)
jj <- bootstrap_intervals(fau330_knm_beta_A, larger_summary_suspect_knm_data)

bootstrap90_CI2 <- list(No_A_km = gg, NO_A_knm = hh, A_km = ii, A_knm = jj)
bootstrap90_CI2
#saveRDS(bootstrap90_CI2, file = "bootstrap90_CI2.rds")

kk <- bootstrap_intervals_diff(fau330_knm_beta_no_A, fau330_knm_beta_A, 
                               larger_summary_test_knm_data, larger_summary_suspect_knm_data,
                               alpha = 0.1)
ll <- bootstrap_intervals_diff(fau330_km_beta_no_A, fau330_km_beta_A, 
                               larger_summary_test_km_data, larger_summary_suspect_km_data,
                               alpha = 0.1)

bootstrap90_CI_diff2 <- list(km_noA_minus_A = kk, knm_noA_minus_A = ll)
bootstrap90_CI_diff2 
#saveRDS(bootstrap90_CI_diff2, file = "bootstrap90_CI_diff2.rds")
# 6. Do similar estimations to more single barrel cases (within or out of this FAU330 barrel)==========


