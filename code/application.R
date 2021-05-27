library(purrr)
library(dplyr)
library(stringr)
library(ggplot2)

theme_set(theme_bw())
# README========================================================
# Ignore most of the comments since they are the thoughts when I made up this, but not directly related to anything
# Some randomly chosen index, sequences, are only meaningful in my initial run. For later run, read in the data.
# The path of data should all be changed to /organized data folder if they are not
# END README====================================================

# We distinguish two cases:
# 1. (database case) use the full data estimation (may leave some for test) to get threshold and likelihood ratio
# 2. (single case) use 3 bullets to get likelihood ratio, threshold
# Expected results: based on the data base estimation, we can get better accuracy (overall).
# This could justify our use of background distribution (database estimation) 
# to get likelihood ratio, threshold and finally accuracy
# And also guide our other analysis

# We need a bullet level likelihood ratio, to compare bullets
# We only need a land level likelihood ratio, to compare lands

## How do we define ground truth
## How do we decide the sample size, size of KM and size of KNM
## Do we have to distinguish matched bullets and non-matched bullets for land comparisons

# 1. We can do cross validation and rotate the samples used to train and test
# 2. Or we can divide the data to training set and testing set at once.
# We want to avoid reusing samples to fit and to compare (but how to achieve this in full data case?)
# we follow 2. for now

# Since we currently have 442 sets of barrels, we use one bullet from each set as a test bullet
# there are 442 test bullets. How do we produce a test set with both matched and non-matched lands?

# Two or three test fires case==========================
# 1. the distributions are estimated only with the test fired bullets (the same one hundred barrels as the other case?)
# 2. prepare KM bullets, and KNM bullets (new data we need to generate)
# 3. compare the lands, and select by SAM, the 18 land comparisons (in KM, they are matches, in KNM, they are those assumed matches (suspected ones))
# 4. find the threshold (case by case), and do classification
# 5. report accuracy
# 6. ROC curve (since the threshold is selected case by case, we may not able to get an ROC here)

# for the data we have now (only KM)================
# we select 342 barrels as training sets for database
set.seed(20201021)
# database_index will be used in the next section, to train the population distribution
database_index <- sample(1:442, size = 342)
database_index <- readRDS("./2020fall/related-data/database_index.rds")
#saveRDS(database_index, "./2020fall/related-data/database_index.rds")
# this set is used in both sections, to train case by case distribution, and get questioned comparisons
test_index <- setdiff(1:442, database_index)
# questioned bullets in the test set, each set has one questioned bullet, used in both cases
questioned_bullet_index <- sample(LETTERS[1:4], 100, replace = T)
questioned_bullet_index <- readRDS("./2020fall/related-data/questioned_bullet_index.rds")
#saveRDS(questioned_bullet_index, "./2020fall/related-data/questioned_bullet_index.rds")

test_km_sets <- subset(ccf_km_full, subset = 1:442 %in% test_index)
test_km_sets[[1]]

test_knm_sets <- subset(ccf_knm_full, subset = 1:442 %in% test_index)
test_knm_sets[[1]]

# test and questioned sets, split each set to tested and questioned bullets
test_questioned_sets <- Map(function(x, y, z) {
  
  test_km <- x %>% filter(bullet1 != z, bullet2 != z)
  test_knm <- y %>% filter(bullet1 != z, bullet2 != z)
  questioned_km <- x %>% filter(bullet1 == z|bullet2 == z)
  
  output <- list(test_km = test_km, test_knm = test_knm,
                 questioned_km = questioned_km)
  
}, test_km_sets, test_knm_sets, questioned_bullet_index)

length(test_questioned_sets)
test_questioned_sets[[1]]

# Used Map instead of pmap:
# a pmap potential bug, as the code reflected
# and a reported issue about function parameter and list names: https://github.com/tidyverse/purrr/issues/203
# could easily lead to mistakes
# example <- function(a, b, c) {
#   output = list(a, b, c)
#   output
# }
# pmap(.l = list(a = 1:3, b = 3:1, c = letters[1:3]), .f = example)
# example <- function(x, y, z) {
#   output = list(x, y, z)
#   output
# }
# pmap(.l = list(a = 1:3, b = 3:1, c = letters[1:3]), .f = example)

# fit the distributions for each set (three test fires)

case_by_case_models <- lapply(test_questioned_sets, FUN = function(x) {
  # the starting values are those from the full data models
  km_model <- mine_betamix(data = x$test_km$ccf, 
                           par = c(0.419, 7.37, 8.45, 14.03, 4.45),
                           k = 2,
                           control = list(maxit = 3000),
                           hessian = FALSE)
  knm_model <- mine_betamix(data = x$test_knm$ccf, 
                            par = c(0.666, 11.23, 20.33, 7.06, 7.22),
                            k = 2,
                            control = list(maxit = 3000),
                            hessian = FALSE)
  output <- list(km_model = km_model, knm_model = knm_model)
  output
})

# KM convergence results, 1 indicates not converged
case_by_case_models %>% map_dbl(.f = function(x) x$km_model$convergence)
# KNM convergence results
case_by_case_models %>% map_dbl(.f = function(x) x$knm_model$convergence)
# The ones we can used: (both knm and km are converged)
converged_index <- case_by_case_models %>% map_dbl(.f = function(x) x$km_model$convergence)==0 &
  case_by_case_models %>% map_dbl(.f = function(x) x$knm_model$convergence) == 0

sum(converged_index)
which(!converged_index)

case_by_case_models[[17]]

# The function: find_cross() is defined in helper-functions.r to find the cross of the two densities
case_by_case_models <- case_by_case_models %>% lapply(FUN = function(x) {
  x$threshold_equal_pdf <- find_cross(x$km_model, x$knm_model,start = 0.2,end = 0.9, step = 0.0001)
  x
})
case_by_case_models[[1]]

# do classification

case_by_case_classifcation <- Map(function(x, y) {
  # True: identification; False: elimination
  class <- x$questioned_km$ccf > y$threshold_equal_pdf
  class
}, test_questioned_sets, case_by_case_models)

case_by_case_classifcation_converged <- subset(case_by_case_classifcation, converged_index)
case_by_case_classifcation_converged %>% unlist() %>% sum
case_by_case_classifcation_converged %>% unlist() %>% length

# have a look at full data model classification, this is much easier

Map(function(x, y) {
  # True: identification; False: eliminations
  class <- x$questioned_km$ccf > 0.5304
  class
}, test_questioned_sets, case_by_case_models) %>% subset(converged_index) %>% unlist() %>% sum() 

# for KNM bullets (new data generated from the server)====================================
# The new comp and extracted ccf are from the server

ccf_nonmatched_km <- readRDS("./2020fall/related-data/ccf_nonmatched_km.rds")
ccf_nonmatched_knm <- readRDS("./2020fall/related-data/ccf_nonmatched_knm.rds")

nonmatched_test_questioned_sets <- Map(function(x, y) {
  
  test_km <- x %>% filter(bullet1 != "Q", bullet2 != "Q")
  test_knm <- y %>% filter(bullet1 != "Q", bullet2 != "Q")
  questioned_km <- x %>% filter(bullet1 == "Q"|bullet2 == "Q")
  
  output <- list(test_km = test_km, test_knm = test_knm,
                 questioned_km = questioned_km)
  
}, ccf_nonmatched_km, ccf_nonmatched_knm)

nonmatched_test_questioned_sets[[1]]

nonmatched_case_by_case_models <- lapply(nonmatched_test_questioned_sets, FUN = function(x) {
  # the starting values are those from the full data models
  km_model <- mine_betamix(data = x$test_km$ccf, 
                           par = c(0.419, 7.37, 8.45, 14.03, 4.45),
                           k = 2,
                           control = list(maxit = 3000),
                           hessian = FALSE)
  knm_model <- mine_betamix(data = x$test_knm$ccf, 
                            par = c(0.666, 11.23, 20.33, 7.06, 7.22),
                            k = 2,
                            control = list(maxit = 3000),
                            hessian = FALSE)
  output <- list(km_model = km_model, knm_model = knm_model)
  output
})


# KM convergence results, 1 indicates not converged
nonmatched_case_by_case_models %>% map_dbl(.f = function(x) x$km_model$convergence)
# KNM convergence results
nonmatched_case_by_case_models %>% map_dbl(.f = function(x) x$knm_model$convergence)
# The ones we can used: (both knm and km are converged)
nonmatched_converged_index <- nonmatched_case_by_case_models %>% map_dbl(.f = function(x) x$km_model$convergence)==0 &
  nonmatched_case_by_case_models %>% map_dbl(.f = function(x) x$knm_model$convergence) == 0

sum(nonmatched_converged_index)
which(!nonmatched_converged_index)

# The function: find_cross() is defined in helper-functions.r to find the threshold
nonmatched_case_by_case_models <- nonmatched_case_by_case_models %>% lapply(FUN = function(x) {
  x$threshold_equal_pdf <- find_cross(x$km_model, x$knm_model,start = 0.1,end = 0.9, step = 0.0001)
  x
})
nonmatched_case_by_case_models[[1]]

# do classification

nonmatched_case_by_case_classifcation <- Map(function(x, y) {
  # True: identification; False: elimination
  class <- x$questioned_km$ccf > y$threshold_equal_pdf
  class
}, nonmatched_test_questioned_sets, nonmatched_case_by_case_models)

nonmatched_case_by_case_classifcation_converged <- subset(nonmatched_case_by_case_classifcation, nonmatched_converged_index)
nonmatched_case_by_case_classifcation_converged %>% unlist() %>% sum
nonmatched_case_by_case_classifcation_converged %>% unlist() %>% length

# have a look at full data model classification

Map(function(x, y) {
  # True: identification; False: eliminations
  class <- x$questioned_km$ccf > 0.5304
  class
}, nonmatched_test_questioned_sets, nonmatched_case_by_case_models) %>% subset(nonmatched_converged_index) %>% unlist() %>% sum() 

1044/1656 - 309/1044

# Database distribution case (partly included in the previous sections, we just re-fitted the model here as a formal section)===================================
# 1. the distributions are estimated with 342 barrels, the other one hundred barrels as testing set
# 2. prepare KM and KNM similar to the previous
# 3. compare the lands, and select by SAM, the 18 land comparisons (in KM, they are matches, in KNM, they are those assumed matches (suspected ones))
# 4. find the threshold, and do classification
# 5. report accuracy
# 6. ROC curve (this then doesn't depend on threshold selection, which is also an important aspect affecting accuracy)

# matched bullets
database_matchedB_km <- mine_betamix(data = ccf_km_full %>% 
                                       subset(subset = !1:442 %in% test_index) %>% 
                                       bind_rows() %>%
                                       as.data.frame() %>%
                                       pull(ccf),
                                     par = c(0.419, 7.36, 8.45, 14, 4.45),
                                     k = 2,
                                     hessian = FALSE,
                                     control = list(maxit = 1000))
database_matchedB_km$converge 
database_matchedB_km

database_matchedB_knm <- mine_betamix(data = ccf_knm_full %>% 
                                        subset(subset = !1:442 %in% test_index) %>% 
                                        bind_rows() %>%
                                        as.data.frame() %>%
                                        filter(ccf<=0.95) %>%
                                        pull(ccf),
                                      par = c(0.666, 11.23, 20.33, 7.065, 7.22),
                                      k = 2,
                                      hessian = FALSE,
                                      control = list(maxit = 1000))


database_matchedB_knm$converge 
database_matchedB_knm

find_cross(m1 = database_matchedB_km, m2 = database_matchedB_knm)
# 0.53
database_matchedB_km_result <- lapply(test_questioned_sets, FUN = function(x) x$questioned_km$ccf > 0.53)

database_matchedB_km_result %>% subset(converged_index) %>% unlist() %>% sum() 


# non-matched bullets
database_nonmatchedB_result <- lapply(nonmatched_test_questioned_sets, FUN = function(x) x$questioned_km$ccf > 0.53)

# Analysis of the results=========================
# Matched bullets, individual cases data:
# case_by_case_classifcation

case_by_case_classifcation_converged <- subset(case_by_case_classifcation, converged_index)
case_by_case_classifcation_converged %>% unlist() %>% sum
case_by_case_classifcation_converged %>% unlist() %>% length

case_by_case_classifcation_converged[[1]] %>% sum /18

data.frame(trueRate = map_dbl(case_by_case_classifcation_converged, .f = function(x) sum(x)/18)) %>%
  ggplot(aes(y = trueRate)) + 
  geom_boxplot()

# database
database_matchedB_km_result_converged <- database_matchedB_km_result %>% subset(converged_index)

data.frame(trueRate = map_dbl(database_matchedB_km_result_converged, .f = function(x) sum(x)/18)) %>%
  ggplot(aes(y = trueRate)) + 
  geom_boxplot()

# in one plot
temp_match_km_cases <- map_dbl(case_by_case_classifcation_converged, .f = function(x) sum(x)/18)
temp_match_km_database <- map_dbl(database_matchedB_km_result_converged, .f = function(x) sum(x)/18)
matched_results <- data.frame(trueRate = c(temp_match_km_cases, temp_match_km_database),
                              class = c(rep("cases", times = length(temp_match_km_cases)),
                                        rep("database", times = length(temp_match_km_database))))
matched_results %>%
  ggplot(aes(y = trueRate, x= class)) + 
  geom_boxplot()


# Non-matched bullets
# individual cases
# nonmatched_case_by_case_classifcation
nonmatched_case_by_case_classifcation_converged <- subset(nonmatched_case_by_case_classifcation, nonmatched_converged_index)
nonmatched_case_by_case_classifcation_converged %>% unlist() %>% sum
nonmatched_case_by_case_classifcation_converged %>% unlist() %>% length

# in one plot

temp_nonmatch_cases <- map_dbl(nonmatched_case_by_case_classifcation_converged, .f = function(x) 1-sum(x)/18)
temp_nonmatch_database <- map_dbl(database_nonmatchedB_result, .f = function(x) 1-sum(x)/18)
matched_results <- data.frame(trueRate = c(temp_nonmatch_cases, temp_nonmatch_database),
                              class = c(rep("cases", times = length(temp_nonmatch_cases)),
                                        rep("database", times = length(temp_nonmatch_database))))
matched_results %>%
  ggplot(aes(y = trueRate, x= class)) + 
  geom_boxplot()

# 















