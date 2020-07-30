library(dplyr)
library(betareg)
library(ggplot2)
# Run the following to get the full data
# ccf_km_full <- land_score_km_generator2("./data/organized data/ccf/ccf_km_full.rds")
# ccf_knm_full <- land_score_knm_generator2("./data/organized data/ccf/ccf_knm_full.rds")

ccf_km_full <- readRDS("./data/organized data/ccf/ccf_km_full.rds")
ccf_knm_full <- readRDS("./data/organized data/ccf/ccf_knm_full.rds")

ccf_km_full_data <- bind_rows(ccf_km_full)
ccf_knm_full_data <- bind_rows(ccf_knm_full)

range(ccf_km_full_data$ccf)
range(ccf_knm_full_data$ccf)

ccf_knm_full_data %>% filter(ccf == 1)
ccf_knm_full_data %>% filter(ccf >= 0.9)

boxplot(ccf_km_full_data$ccf)
boxplot(ccf_knm_full_data$ccf)
