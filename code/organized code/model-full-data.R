# Make use of the RF data to fit the full data model, and generate some figures and tables
# Make sure to use your path to read and save

library(dplyr)
library(ggplot2)
library(betareg)
library(kableExtra)
theme_set(theme_bw())
# read in and prepare the data =========
# use path like: "/media/Raven/LAPD-processing/data/rfscore/rf_km_full.rds"
rf_km_full <- readRDS("./organized data/rfscore/rf_km_full.rds")
rf_knm_full <- readRDS("./organized data/rfscore/rf_knm_full.rds")

# data processing
rf_km_full_data <- bind_rows(rf_km_full)
rf_knm_full_data <- bind_rows(rf_knm_full)

# note, there are 0 and 1s in the rfscore which can't be used to fit beta distribution
# We do adjustment, set 1 as 0.997 (round up), 0 as 0.006 (round down) 
# (the quantiles can provide some support for doing this)
range(rf_km_full_data$rfscore)
range(rf_knm_full_data$rfscore)
boxplot(rf_km_full_data$rfscore)
boxplot(rf_knm_full_data$rfscore)
quantile(rf_km_full_data$rfscore, probs = c(0.95, 0.94, 0.93, 0.92, 0.91))
quantile(rf_knm_full_data$rfscore, probs = c(0.005, 0.01))

rf_km_full_data$rfscore <- ifelse(rf_km_full_data$rfscore == 1, yes = 0.997, no = rf_km_full_data$rfscore)
rf_km_full_data$rfscore <- ifelse(rf_km_full_data$rfscore == 0, yes = 0.006, no = rf_km_full_data$rfscore)
rf_knm_full_data$rfscore <- ifelse(rf_knm_full_data$rfscore == 1, yes = 0.997, no = rf_knm_full_data$rfscore)
rf_knm_full_data$rfscore <- ifelse(rf_knm_full_data$rfscore == 0, yes = 0.006, no = rf_knm_full_data$rfscore)

# draw some histograms
ggplot(rf_km_full_data) + 
  geom_histogram(aes(x = rfscore, y = ..ndensity..), bins = 20, fill = "blue4", color = "gray99") + 
  ylab("relative frequency")

ggsave("./organized code/figures/rf_km_full_2components.png")

ggplot(rf_knm_full_data) + 
  geom_histogram(aes(x = rfscore, y = ..ndensity..), bins = 20, fill = "blue4", color = "gray99") + 
  ylab("relative frequency")

ggsave("./organized code/figures/rf_knm_full_2components.png")

ggplot(bind_rows(rf_km_full_data, rf_knm_full_data, .id = "Class")) + 
  geom_histogram(aes(x = rfscore, y = ..density.., fill = Class), bins = 20, color = "gray99",
                 position = "identity", alpha = 0.5) + 
  ylab("density") +
  scale_fill_discrete(breaks = c("1", "2"), label = c("KM", "KNM"))

ggsave("./organized code/figures/rf_full_2components_combined.png")

# fit the model (two components for both, 5 mins, 20 mins)===========================
rf_km_full_model <- betamix(rfscore~1|1, data = rf_km_full_data, k = 2)
rf_knm_full_model <- betamix(rfscore~1|1, data = rf_knm_full_data, k = 2)

#saveRDS(rf_km_full_model, file = "./organized data/models/rf_km_full_model.rds")
#saveRDS(rf_knm_full_model, file = "./organized data/models/rf_knm_full_model.rds")

# bootstrap CIs for this estimation is too burdensome and may not be realistic(couple weeks for 1000sample(one core))

# it assigns observations to each component, and all observations are assigned to one group
# But this is decided by the "posteriors" of the components, it actually just compare p1 and p2, however, we combine them
rf_km_full_model$flexmix
rf_knm_full_model$flexmix

rf_km_full_model$flexmix@posterior$scaled
rf_km_full_model$flexmix@logLik
rf_knm_full_model$flexmix@logLik

# Draw the estimated distribution (there are helper functions used here)======================
rf_km_full_simulated <- simulate_bitamix(rf_km_full_model)
rf_knm_full_simulated <- simulate_bitamix(rf_knm_full_model)
rf_full_simulated_combined <- bind_rows(rf_km_full_simulated, rf_knm_full_simulated, .id = "Class")
colnames(rf_full_simulated_combined)[2] <- "rfscore"

rf_full_simulated_combined %>%
  ggplot(aes(color = Class, fill = Class)) + 
  geom_line(aes(x = rfscore, y = y)) + 
  geom_histogram(aes(x = rfscore, y = ..density..), bins = 20, color = "gray99",
                 position = "identity", alpha = 0.5, 
                 data = bind_rows(rf_km_full_data, rf_knm_full_data, .id = "Class")) +
  scale_color_discrete(breaks = c("1", "2"), label = c("KM", "KNM")) + 
  scale_fill_discrete(breaks = c("1", "2"), label = c("KM", "KNM")) +
  ylab("density")

ggsave("./organized code/figures/rf_full_2components_estimated.png")

# Inference: asymptotic (EM) inference, bootstrap, likelihood ratio for sufficiency?=================
# Details shown in each section if available
# wonder if direct computation of information matrix is possible since gamma function involved (seem to be available)
# But should we use the same way for both large and small sample sizes? 
# (bootstrap burdensome for large data, wonder if information matrix is also good for relatively small data?)

# See if beta sufficient by likelihood ratio==========================
rf_km_full_beta_model <- betamix(rfscore~1|1, data = rf_km_full_data, k = 1)
rf_knm_full_beta_model <- betamix(rfscore~1|1, data = rf_knm_full_data, k = 1)

rf_km_full_beta_model$flexmix@logLik
rf_knm_full_beta_model$flexmix@logLik

# KM not significant
-2*(rf_km_full_beta_model$flexmix@logLik-rf_km_full_model$flexmix@logLik)
# KNM significant
-2*(rf_knm_full_beta_model$flexmix@logLik-rf_knm_full_model$flexmix@logLik)

# draw the plot in this single beta case
rf_km_full_beta_simulated <- simulate_bitamix(rf_km_full_beta_model)
rf_knm_full_beta_simulated <- simulate_bitamix(rf_knm_full_beta_model)
rf_full_beta_simulated_combined <- bind_rows(rf_km_full_beta_simulated, rf_knm_full_beta_simulated, .id = "Class")
colnames(rf_full_beta_simulated_combined)[2] <- "rfscore"

rf_full_beta_simulated_combined %>%
  ggplot(aes(color = Class, fill = Class)) + 
  geom_line(aes(x = rfscore, y = y)) + 
  geom_histogram(aes(x = rfscore, y = ..density..), bins = 20, color = "gray99",
                 position = "identity", alpha = 0.5, 
                 data = bind_rows(rf_km_full_data, rf_knm_full_data, .id = "Class")) +
  scale_color_discrete(breaks = c("1", "2"), label = c("KM", "KNM")) + 
  scale_fill_discrete(breaks = c("1", "2"), label = c("KM", "KNM")) +
  ylab("density")

ggsave("./organized code/figures/rf_full_beta_estimated.png")

# Make table tables for the estimation and tests (SE and CI for each estimate to be produced)=======
# How to test if betamix is sufficient w.r.t the saturated model(what is it?) 
beta_params(rf_km_full_model)

table_rf_full_model <- extract_params(cv_summary(list(rf_km_full_model, rf_knm_full_model))) 
rownames(table_rf_full_model ) <- c("rf_km_full_model", "rf_knm_full_model")

table_rf_full_model %>%
  round(digits = 2) %>% 
  kable() %>% 
  kable_styling(bootstrap_options = c("striped", "hover"))

# some tables directly created in RMD/latex with the estimated values here
















