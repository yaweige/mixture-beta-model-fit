head(km_ccf)
head(knm_ccf)
dim(km_ccf)
dim(knm_ccf)


# Use the data we have currently, to fit the full data model, compare this with the single barrel model we just evaluated

# This takes much longer time than the single barrel fit
full_km_model <- betamix(ccf~1|1, data = km_ccf, k = 2)
full_knm_model <- betamix(ccf~1|1, data = knm_ccf, k = 2)

saveRDS(full_km_model, "full_km_model.rds")
saveRDS(full_knm_model, "full_knm_model.rds")

system.time(betamix(ccf~1|1, data = km_ccf, k = 2))
# about 72 seconds, 1000 samples takes 20 hours

# Do a bootstrap inference similar to the single barrel case (This run on the server with 3000 samples each)
# summary_full_km <- summary_simulated(full_km_model, nrun = 99, n_each_run = 2232)
# summary_full_knm <- summary_simulated(full_knm_model, nrun = 99, n_each_run = 3276)

# read in the data from server (this ran in two different times, just combine them)

summary_full_km1 <- readRDS("./upload/summary_full_km.rds")
summary_full_km2 <- readRDS("./upload/summary_full_km2.rds")

summary_full_km <- c(summary_full_km1, summary_full_km2)
rm(summary_full_km1)
rm(summary_full_km2)

summary_full_knm1 <- readRDS("./upload/summary_full_knm.rds")
summary_full_knm2 <- readRDS("./upload/summary_full_knm2.rds")

summary_full_knm <- c(summary_full_knm1, summary_full_knm2)
rm(summary_full_knm1)
rm(summary_full_knm2)

# Detect the failed cases
summary_full_km_nonfailed <- lapply(summary_full_km, FUN = is.list) %>% unlist()
sum(summary_full_km_nonfailed)

summary_full_knm_nonfailed <- lapply(summary_full_knm, FUN = is.list) %>% unlist()
sum(summary_full_knm_nonfailed)
# All succeeded

summary_full_km_data <- extract_params(summary_full_km)
summary_full_knm_data <- extract_params(summary_full_knm)

mm <- bootstrap_intervals(full_km_model, summary_x_x_data = summary_full_km_data)
nn <- bootstrap_intervals(full_knm_model, summary_x_x_data = summary_full_knm_data)

bootstrap_full_90_CI <- list(full_km = mm, full_knm = nn)
#saveRDS(bootstrap_full_90_CI, "bootstrap_full_90_CI.rds")

# the difference

oo <- bootstrap_intervals_diff(model1 = fau330_km_beta_no_A, model2 = full_km_model, 
                               sampled_data1 = larger_summary_test_km_data, sampled_data2 = summary_full_km_data)
pp <- bootstrap_intervals_diff(model1 = fau330_knm_beta_no_A, model2 = full_knm_model, 
                               sampled_data1 = larger_summary_test_knm_data, sampled_data2 = summary_full_knm_data)

bootstrap_full_90_CI_diff <- list(km_330_minus_full = oo, knm_330_minus_full = pp)
#saveRDS(bootstrap_full_90_CI_diff, "bootstrap_full_90_CI_diff.rds")
# I suspect the order adjustment is not correct in cv_summary (checked, this is correct)
# fau330_km_beta_no_A is kind of special, that it has larger p1(reversed), it has second componets with even larger mean value
# we should generate this to more single barrel case.

rr <- bootstrap_intervals_diff(model1 = fau330_km_beta_A, model2 = full_km_model, 
                               sampled_data1 = larger_summary_suspect_km_data, sampled_data2 = summary_full_km_data)
ss <- bootstrap_intervals_diff(model1 = fau330_knm_beta_A, model2 = full_knm_model, 
                               sampled_data1 = larger_summary_suspect_knm_data, sampled_data2 = summary_full_knm_data)

bootstrap_full_90_CI_diff_withA <- list(km_330withA_minus_full = rr, knm_330WithA_minus_full = ss)
#saveRDS(bootstrap_full_90_CI_diff_withA, "bootstrap_full_90_CI_diff_withA.rds")
# summarize, and make a plan for next step

# summarize the models we estimated










