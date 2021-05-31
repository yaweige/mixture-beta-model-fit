library(dplyr)
# 36 each for km, 240 each for knm
# within barrel comparisons===================================================
ccf_km_full <- readRDS("./data/ccf/ccf_km_2nd.rds")
ccf_knm_full <- readRDS("./data/ccf/ccf_knm_2nd.rds")

ccf_km_full_data <- bind_rows(ccf_km_full)
ccf_knm_full_data <- bind_rows(ccf_knm_full)

ccf_km_full_data$source <- "KM"
ccf_km_full_data$usage <- "part1"

ccf_knm_full_data$source <- "KNM"
ccf_knm_full_data$usage <- "part1"

ccf_full_data <- bind_rows(ccf_km_full_data,
                           ccf_knm_full_data)

# non-matching bullet comparisons for test====================================
ccf_nonmatched_km <- readRDS("./data/ccf/ccf_nonmatched_km.rds")
ccf_nonmatched_knm <- readRDS("./data/ccf/ccf_nonmatched_knm.rds")

ccf_nonmatched_km_data <- ccf_nonmatched_km %>% bind_rows()
ccf_nonmatched_km_data$source <- "KM"
ccf_nonmatched_km_data$usage <- "part2"

ccf_nonmatched_knm_data <- ccf_nonmatched_knm %>% bind_rows()
ccf_nonmatched_knm_data$source <- "KNM"
ccf_nonmatched_knm_data$usage <- "part2"

ccf_nonmatched_data <- bind_rows(ccf_nonmatched_km_data,
                                 ccf_nonmatched_knm_data)

# save as csv======================

ccf_data <- bind_rows(ccf_full_data, ccf_nonmatched_data)
write.csv(ccf_data, file = "./data/ccf/ccf_data.csv", row.names = F)




