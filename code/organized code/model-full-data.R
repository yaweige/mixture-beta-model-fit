# Make use of the RF data to fit the full data model, and generate some figures and tables
# Make sure to use your path to read and save

library(dplyr)
library(betareg)
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

ggsave("./organized code/figures/rf_knm_full_2components_combined.png")

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

# Draw the estimated distribution (there are helper functions used here)======================





