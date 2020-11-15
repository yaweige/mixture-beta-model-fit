library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(gridExtra)

theme_set(theme_bw())

# README==================================================
# This scripts are under active changes, and some latest parts may not included
# Ignore most of the comments since they are not directly related to anything
# README END====================================================

# we have 442 sets of bullets in our data

# want to know:
# 1. how big a sample size is needed to reach a decent distribution estimation (stable)?
# 2. how the error rates change accordingly?

# 1 barrel, 2 barrel, 3 barrel, 4 barrel, 5 barrel (80 random draws for each)
# 4 bullets, 2*4 bullets, 3*4bullets... 
# So this is not to double the data, we need some assumptions (especially for KM), 
# what's our original assumption when estimating the full data distribution
# We need another measure of needed known source evidence? (We just use the number of sets)
# calculate the sample sd of each parameter, visual difference
# find the thresholds, and quantify the variation
# find the errors rates, quantify the variation

# conclude the sensitivity with respect to sample size, make suggestion to necessary sample size needed

length(ccf_km_full)
length(ccf_knm_full)
dim(ccf_km_full_data)
dim(ccf_knm_full_data)
ccf_km_full[[1]]
ccf_knm_full[[1]]

# parameter transformation (alpha, beta) -> (mu, phi), for k = 2
param_trans <- function(param) {
  p1 <- param[1]
  mu1 <- param[2] / (param[2] + param[3])
  phi1 <- param[2] + param[3]
  mu2 <- param[4] / (param[4] + param[5])
  phi2 <- param[4] + param[5]
  
  output <- c(p1, mu1, phi1, mu2, phi2)
  output
}

# to simulate and draw estimated curves
# ab_param = c(p1, a1, b1, a2, b2)
simulate_betamix_c2 <- function(ab_param) {
  x <- 1:99/100
  p1 <- ab_param[1]
  a1 <- ab_param[2]
  b1 <- ab_param[3]
  a2 <- ab_param[4]
  b2 <- ab_param[5]
  
  y1 <- dbeta(x, shape1 = a1, shape2 = b1)
  y2 <- dbeta(x, shape1 = a2, shape2 = b2)
  y <- p1*y1 + (1-p1)*y2
  output <- data.frame(x = x, y1 = y1, y2 = y2, y = y)
  output
}

# The KM run (with hessian matrix not estimated, use the same starting values as in full data)==========

# 1, 2, 3, 4, 5 barrels===================================

# output a list, 
# [[1]] model, [[2]] parameter data, [[3]] simulated data, 
# [[4]] boxplot matrix, [[5]] curves, [[6]] averaged and 1-sd-band curve

set.seed(202002021)
index_km_n1 <- sample(1:442, 80, replace = FALSE)
index_km_n2 <- sample(1:442, 160, replace = FALSE)
index_km_n3 <- sample(1:442, 240, replace = FALSE)
index_km_n4 <- sample(1:442, 320, replace = FALSE)
index_km_n5 <- sample(1:442, 400, replace = FALSE)

# index: samples used, length of index = neach*nsample
# neach: number of barrels used for each estimation
# nsample: total number of estimations
changing_size_function <- function(data, index, neach, nsample) {
  
  if (neach*nsample > length(data)) stop("used samples exceeded the size of the data")
  
  ccf_km_c2_n2 <- lapply(1:nsample, FUN = function(i) {
    b <- index %>% split(rep(1:nsample, each = neach)) %>% `[[`(i)
    
    selected_data <- vector()
    selected_data <- lapply(1:neach, FUN = function(x) {
      output <- c(selected_data, data[[b[x]]]$ccf)
    }) %>% unlist()
    
    model <- mine_betamix(data = selected_data, 
                          par = c(0.3, 5, 5, 5, 2),
                          k = 2,
                          control = list(maxit = 3000),
                          hessian = FALSE)
    print(i)
    model
  })
  
  # find the converged estimation: all converged
  num_of_unconverged <- map(ccf_km_c2_n2, "convergence") %>% unlist() %>% sum()
  print("num of converged out of all samples: ")
  print(nsample-num_of_unconverged)
  
  if(num_of_unconverged != 0) {
    print("index of un-converged sample, data[[x]] to get the data:")
    print(index[which(map(ccf_km_c2_n2, "convergence") %>% unlist() != 0)])
    print("un-converged sample in the output, output[[x]] to get the estimation:")
    print(which(map(ccf_km_c2_n2, "convergence") %>% unlist() != 0))
    print("further analysis not done, only estimation returned")
    return(ccf_km_c2_n2)
  }
  
  convergence_n2 <- map(ccf_km_c2_n2, "convergence") %>% unlist()
  ccf_km_c2_n2_converged <- subset(ccf_km_c2_n2, convergence_n2 == 0)
  map(ccf_km_c2_n2_converged, "convergence") %>% unlist() %>% sum()
  
  # output a data.frame with all those needed
  ccf_km_c2_n2_data <- lapply(1:5, FUN = function(i) {
    column <- ccf_km_c2_n2_converged %>% map("par") %>% map_dbl(.f = function(x) x[i])
    column
  })  %>% as.data.frame() %>% set_names(c("p1", "a1", "b1", "a2", "b2"))
  
  temp_data <- lapply(1:5, FUN = function(i) {
    column <- ccf_km_c2_n2_converged %>% map("par") %>% map(param_trans) %>% map_dbl(.f = function(x) x[i])
    column
  })  %>% as.data.frame() %>% set_names(c("p1", "mu1", "phi1", "mu2", "phi2"))
  
  ccf_km_c2_n2_data <- bind_cols(ccf_km_c2_n2_data, temp_data %>% select(-p1))
  ccf_km_c2_n2_data <- ccf_km_c2_n2_data %>%
    mutate(mu = p1*mu1 + (1-p1)*mu2)
  
  ccf_km_c2_n2_data_long <- ccf_km_c2_n2_data %>% gather()
  p1 <- ccf_km_c2_n2_data_long %>% 
    ggplot(aes(y =value)) + 
    geom_boxplot() + 
    facet_wrap(~key, scales = "free")
  
  # to plot
  n2_simulated_data <- lapply(1:nrow(ccf_km_c2_n2_data), FUN = function(i) {
    output <- simulate_betamix_c2(ccf_km_c2_n2_data[i,] %>% as.numeric())
    output
  }) %>% bind_rows(.id = "sampleID")
  
  p2 <- n2_simulated_data %>%
    ggplot(aes(x = x, y = y, group = rep(1:nsample, each = 99))) + 
    geom_line()
  
  # averaged curve (average of the curves (densities), and a 2*(1.96?)*sd band)
  n2_simulated_avg_data <- n2_simulated_data %>%
    group_by(avg_id = rep(1:99, times = nsample)) %>%
    summarise(x = mean(x),
              y_avg = mean(y),
              y_sd = sd(y))
  
  p3 <- n2_simulated_avg_data %>%
    mutate(y_low = y_avg - y_sd,
           y_up = y_avg + y_sd) %>%
    ggplot() + 
    geom_line(aes(x = x, y = y_avg)) + 
    geom_line(aes(x = x, y = y_low), color = "blue", linetype = "dashed") + 
    geom_line(aes(x = x, y = y_up), color = "blue", linetype = "dashed")
  
  output <- list(model = ccf_km_c2_n2,
                 param = ccf_km_c2_n2_data,
                 simulation = n2_simulated_data,
                 boxplot = p1,
                 curves = p2,
                 averaged = p3)
  output
}

# See all the one barrel case, some may fail to converge: only one barrel failed to converge
test_one_barrel <- changing_size_function(data = ccf_km_full,
                                          index = 1:442,
                                          neach = 1,
                                          nsample = 442)
test_one_barrel[[8]]

ccf_km_c2_n1_results <- changing_size_function(data = ccf_km_full,
                                               index = index_km_n1,
                                               neach = 1,
                                               nsample = 80)
ccf_km_c2_n1_results$boxplot
ccf_km_c2_n1_results$curves
ccf_km_c2_n1_results$averaged

ccf_km_c2_n2_results <- changing_size_function(data = ccf_km_full,
                                               index = index_km_n2,
                                               neach = 2,
                                               nsample = 80)
ccf_km_c2_n2_results$boxplot
ccf_km_c2_n2_results$curves
ccf_km_c2_n2_results$averaged

ccf_km_c2_n3_results <- changing_size_function(data = ccf_km_full,
                                               index = index_km_n3,
                                               neach = 3,
                                               nsample = 80)
ccf_km_c2_n3_results$boxplot
ccf_km_c2_n3_results$curves
ccf_km_c2_n3_results$averaged

ccf_km_c2_n4_results <- changing_size_function(data = ccf_km_full,
                                               index = index_km_n4,
                                               neach = 4,
                                               nsample = 80)
ccf_km_c2_n4_results$boxplot
ccf_km_c2_n4_results$curves
ccf_km_c2_n4_results$averaged

ccf_km_c2_n5_results <- changing_size_function(data = ccf_km_full,
                                               index = index_km_n5,
                                               neach = 5,
                                               nsample = 80)
ccf_km_c2_n5_results$boxplot
ccf_km_c2_n5_results$curves
ccf_km_c2_n5_results$averaged

# more analysis on the results=========================================
# Based on our assumption, the estimated distributions are realizations of the one with full data (or the underlying one)
# inspired by the plots of Q1, should we define a criteria for convergece? or is there any? in terms of samples used
# Question 1. Are the two components same? Is beta distribution a more reasonable choice?

# Draw (mu1, mu2) colored by p1
mu1mu2p1 <- bind_rows(ccf_km_c2_n1_results$param, 
                      ccf_km_c2_n2_results$param, 
                      ccf_km_c2_n3_results$param, 
                      ccf_km_c2_n4_results$param, 
                      ccf_km_c2_n5_results$param, .id = "neach")

mu1mu2p1 %>%
  ggplot(aes(x = mu1, y = mu2, color = p1)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  facet_wrap(~neach) + 
  scale_color_gradient2(mid = "gray", midpoint = 0.5)
# observations from the above plot
# 1. clearly, there are two groups of points, and a few points on the top and left boundary
# 2. call the upper group as group 1, the lower group as group 2. The group 2 has clearly evidence supporting only beta
# 3. mu1 and mu2 has a clear linear relationship
# 4. For group 1, blue on top, gray in the middle, red in the bottom
# 5. red points are more than blue points
# 6. The above observations, are clearer with more data (except the two group division)
# How about other dark red blue and dark blue in group 1
# Can we conclude that for one barrel case, the majority of weight (p1, p2), goes to c1 or c2, or p1 around 0.5, c1 and c2 similar?
# It seems even for one barrel cases, there are sometimes two modes
mu1mu2p1 %>%
  ggplot(aes(x = phi1, y = phi2, color = p1)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  facet_wrap(~neach) + 
  xlim(0,100) +
  ylim(0,100) +
  scale_color_gradient2(mid = "gray", midpoint = 0.5)


# Question 2. How to interpret the curves plotted, obviously, they don't look the same?


# Question 3. Why 5 barrel estimation is better than 1 barrel estimation? In what sense, how to measure it?

# Analysis of the results=============================================================

# plot1: variation of mean of overall, c1 and c2 for KM and KNM

all5_params <- bind_rows(ccf_km_c2_n1_results$param, 
                         ccf_km_c2_n2_results$param, 
                         ccf_km_c2_n3_results$param, 
                         ccf_km_c2_n4_results$param, 
                         ccf_km_c2_n5_results$param,
                         .id = "size")

all5_params <- all5_params %>%
  mutate(var1 = a1*b1/(a1+b1)^2/(a1+b1+1),
         var2 = a2*b2/(a2+b2)^2/(a2+b2+1),
         var = p1^2*var1 + (1-p1)^2*var2)

p1 <- all5_params %>%
  ggplot(aes(y = mu1, x = size)) + 
  geom_boxplot() + 
  facet_grid(.~"mu1")

p2 <- all5_params %>%
  ggplot(aes(y = mu2, x = size)) + 
  geom_boxplot() + 
  facet_grid(.~"mu2")

p3 <- all5_params %>%
  ggplot(aes(y = mu, x = size)) + 
  geom_boxplot() + 
  facet_grid(.~"mu")

grid.arrange(p1, p2, p3, ncol=1)

# plot2: variation of variance of overall, c1 and c2 for KM and KNM

p4 <- all5_params %>%
  ggplot(aes(y = var1, x = size)) + 
  geom_boxplot() + 
  facet_grid(.~"var1")

p5 <- all5_params %>%
  ggplot(aes(y = var2, x = size)) + 
  geom_boxplot() + 
  facet_grid(.~"var2")

p6 <- all5_params %>%
  ggplot(aes(y = var, x = size)) + 
  geom_boxplot() + 
  facet_grid(.~"var")

grid.arrange(p4, p5, p6, ncol=1)

# plot3: variation of p1 for KM and KNM

all5_params %>%
  ggplot(aes(y = p1, x = size)) + 
  geom_boxplot() + 
  facet_grid(.~"p1")
