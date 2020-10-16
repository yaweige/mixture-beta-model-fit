library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

theme_set(theme_bw())
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

# The KM run (with hessian matrix not estimated, use the same starting values as in full data)==========
# 1 barrel=====================================
# select the random sample
set.seed(202002021)
index_km_n1 <- sample(1:442, 80, replace = FALSE)

ccf_km_c2_n1 <- lapply(index_km_n1, FUN = function(i) {
  selected_data <- ccf_km_full[[i]]
  model <- mine_betamix(data = selected_data$ccf, 
                        par = c(0.3, 5, 5, 5, 2),
                        k = 2,
                        control = list(maxit = 3000),
                        hessian = FALSE)
  print(i)
  model
})

# find the converged estimation: all converged
map(ccf_km_c2_n1, "convergence") %>% unlist() %>% sum()
convergence_n1 <- map(ccf_km_c2_n1, "convergence") %>% unlist()
ccf_km_c2_n1_converged <- subset(ccf_km_c2_n1, convergence_n1 == 0)
map(ccf_km_c2_n1_converged, "convergence") %>% unlist() %>% sum()

# get the estimation, should order the components, if we want to compare a single parameter
ccf_km_c2_n1_con_param <- ccf_km_c2_n1_converged %>% map("par") %>% map(param_trans)
hist(map_dbl(ccf_km_c2_n1_con_param, .f = function(x) x[1]))
plot(density(map_dbl(ccf_km_c2_n1_con_param, .f = function(x) x[1])))
boxplot(map_dbl(ccf_km_c2_n1_con_param, .f = function(x) x[1]))
boxplot(map_dbl(ccf_km_c2_n1_con_param, .f = function(x) x[2]))
boxplot(map_dbl(ccf_km_c2_n1_con_param, .f = function(x) x[3]))
boxplot(map_dbl(ccf_km_c2_n1_con_param, .f = function(x) x[4]))
boxplot(map_dbl(ccf_km_c2_n1_con_param, .f = function(x) x[5]))

# output a data.frame with all those needed
ccf_km_c2_n1_data <- lapply(1:5, FUN = function(i) {
  column <- ccf_km_c2_n1_converged %>% map("par") %>% map_dbl(.f = function(x) x[i])
  column
  })  %>% as.data.frame() %>% set_names(c("p1", "a1", "b1", "a2", "b2"))

temp_data <- lapply(1:5, FUN = function(i) {
  column <- ccf_km_c2_n1_converged %>% map("par") %>% map(param_trans) %>% map_dbl(.f = function(x) x[i])
  column
})  %>% as.data.frame() %>% set_names(c("p1", "mu1", "phi1", "mu2", "phi2"))

ccf_km_c2_n1_data <- bind_cols(ccf_km_c2_n1_data, temp_data %>% select(-p1))
ccf_km_c2_n1_data <- ccf_km_c2_n1_data %>%
  mutate(mu = p1*mu1 + (1-p1)*mu2)

ccf_km_c2_n1_data_long <- ccf_km_c2_n1_data %>% gather()
ccf_km_c2_n1_data_long %>% 
  ggplot(aes(y =value)) + 
  geom_boxplot() + 
  facet_wrap(~key, scales = "free")

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

n1_simulated_data <- lapply(1:nrow(ccf_km_c2_n1_data), FUN = function(i) {
  output <- simulate_betamix_c2(ccf_km_c2_n1_data[i,] %>% as.numeric())
  output
}) %>% bind_rows(.id = "sampleID")

n1_simulated_data %>%
  ggplot(aes(x = x, y = y, group = rep(1:80, each = 99))) + 
  geom_line()

# averaged curve (average of the curves (densities), and a 2*(1.96?)*sd band)
n1_simulated_avg_data <- n1_simulated_data %>%
  group_by(avg_id = rep(1:99, times = 80)) %>%
  summarise(x = mean(x),
            y_avg = mean(y),
            y_sd = sd(y))

n1_simulated_avg_data %>%
  mutate(y_low = y_avg - y_sd,
         y_up = y_avg + y_sd) %>%
  ggplot() + 
  geom_line(aes(x = x, y = y_avg)) + 
  geom_line(aes(x = x, y = y_low), color = "blue", linetype = "dashed") + 
  geom_line(aes(x = x, y = y_up), color = "blue", linetype = "dashed")

temp <- n1_simulated_avg_data %>%
  mutate(y_low = y_avg - y_sd,
         y_up = y_avg + y_sd)

n1_simulated_band <- data.frame(rep(1:99/100, times = 2), c(temp$y_low, rev(temp$y_up)))


p1 <- n1_simulated_avg_data %>%
  mutate(y_low = y_avg - y_sd,
         y_up = y_avg + y_sd) %>%
  ggplot() + 
  geom_line(aes(x = x, y = y_avg)) + 
  geom_line(aes(x = x, y = y_low), color = "blue", linetype = "dashed") + 
  geom_line(aes(x = x, y = y_up), color = "blue", linetype = "dashed") + 
  geom_ribbon(aes(x = x, ymin = y_low, ymax = y_up), alpha = 0.3) + 
  geom_point(aes(x = mu, y = 1), data = ccf_km_c2_n1_data, alpha = 0)

ggExtra::ggMarginal(
  p = p1,
  type = 'boxplot',
  margins = 'x',
  size = 5,
  colour = 'black',
  fill = 'gray'
)

# 2 barrels






