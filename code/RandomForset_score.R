# Make use of the full comp data, generate and compare the Random Forest scores to mimic the real case, see the performance

# We're considerring compare the score based likelihood ratio for 
# 1) no distributional assumption? What exactly is this
# 2) beta distribution
# 3) beta mix distribution

library(randomForest)
library(bulletxtrctr)
library(pROC)

# To generate the RF scores we will use=======================================

# Use land_score_km_generator3_RF and land_score_knm_generator3_RF (in helper functions)

# Run the following to get the full data
# rf_km_full <- land_score_km_generator3_RF("~/comp_full_manual_grooves_notMine/")
# rf_knm_full <- land_score_knm_generator3_RF("~/comp_full_manual_grooves_notMine/")

rf_km_full[[1]]


# Leave one out procedure(For now I would not reuse any barrel, i.e. only leave bullet A out)==============================
# Since RF scores are different from ccf, we will draw some plots for the RF scores latter

# Leave A out
onebarrel_km_example <- rf_km_full[[1]]
onebarrel_example_km_noA <- onebarrel_km_example %>% filter(bullet2 != "A", bullet1 != "A")


onebarrel_knm_example <- rf_knm_full[[1]]
onebarrel_example_knm_noA <- onebarrel_knm_example %>% filter(bullet2 != "A", bullet1 != "A")


# All data with A
onebarrel_example_km_A <- onebarrel_km_example %>% filter(bullet1 == "A"|bullet2 == "A")
onebarrel_example_km_A$groundtruth <- "match"

# exclude comparison within A for KNM (for KM, this has been excluded already)
onebarrel_example_knm_A <- onebarrel_knm_example %>% 
  filter(bullet1 == "A"|bullet2 == "A", !bullet1 == "A"&bullet2 == "A")
  
onebarrel_example_knm_A$groundtruth <- "nonmatch"

onebarrel_example_A <- bind_rows(onebarrel_example_km_A, onebarrel_example_knm_A)

# Beta Fit============================================================================
km_noA_beta <-  betamix(rfscore~1|1, data = onebarrel_example_km_noA, k = 1)
knm_noA_beta <- betamix(rfscore~1|1, data = onebarrel_example_knm_noA, k = 1)


# Beta Mix Fit=========================================================================
km_noA_betamix <-  betamix(rfscore~1|1, data = onebarrel_example_km_noA, k = 2)
knm_noA_betamix <-  betamix(rfscore~1|1, data = onebarrel_example_knm_noA, k = 2)

# To get a cutoff value from the Known comparisons, i.e. noA=================================
# A cutoff for the SLR, a cutoff for the distribution itself(with more than one possible nominal type 1 errors).
# How did everyone make decisions with SLR
# Shall we control the cutoff selection criteria, such as: always equate specificity and sensitivity
# These are impportant questions to be answered

# In land level, in bullet level?

# we might do
# 1. get the observed SLR under different models, get the equal error cutoff, evaluate the performance on the test bullet A
# 2. direct compare the score, get equal error cutoff, evaluate the performance on the test bullet A?
# Any other way to make use of the distributions? Are these actually same?


# example, re-run this latter
noA_mm <- loglr_generator_general(c(onebarrel_example_km_noA$rfscore, onebarrel_example_knm_noA$rfscore), 
                                  km_model = list(km_noA_betamix), knm_model = list(knm_noA_betamix))
noA_mm_roc <- roc(response = c(rep(1, 18), rep(0, 135)), predictor = noA_mm)
plot(noA_mm_roc)


temp <- coords(noA_mm_roc, x = "all", transpose = F) 
temp
# Note! This is not the best we can choose from the ROC, need more work, 
# since sometimes one of specificity, sensitivity constant, the other keeps improving (since discrete)
temp[which.min(abs(temp$specificity - temp$sensitivity)),]
# threshold specificity sensitivity
# -1.543217   0.8888889   0.8888889

eer_extractor <- function(data, response, km_model, knm_model) {
  noA_mm <- loglr_generator_general(data, 
                                    km_model = list(km_model), knm_model = list(knm_model))
  noA_mm_roc <- roc(response = response, predictor = noA_mm)
  temp <- coords(noA_mm_roc, x = "all", transpose = F) 
  output <- temp[which.min(abs(temp$specificity - temp$sensitivity)),]
  c(output[1,], AUC = noA_mm_roc$auc)
}



# There are three possible combinations of the estimated models we will explore
# (hard to imagine that betamix outperforms beta...)
# 1. KM: betamix, KNM: betamix
# 2. KM: beta,    KNM: betamix
# 3. KM: betamix, KNM: beta
# 4. KM: beta,    KNM: beta
# And direct comparison (how to construct straight up SLR without distributional assumption? Some data driven, bootstrap thing?)

# 1. KM: betamix, KNM: betamix========================================

# rfscore = 0 & 1 replaced by 0.005, 0.995
mm_cutoff <- eer_extractor(c(onebarrel_example_km_noA$rfscore, onebarrel_example_knm_noA$rfscore), 
              response = c(rep(1, 18), rep(0, 135)),
              km_model = km_noA_betamix, knm_model = knm_noA_betamix)
mm_cutoff 

A_mm <- loglr_generator_general(onebarrel_example_A$rfscore, km_model = list(km_noA_betamix), knm_model = list(knm_noA_betamix))

# 2. KM: beta, KNM: betamix========================================
bm_cutoff <- eer_extractor(c(onebarrel_example_km_noA$rfscore, onebarrel_example_knm_noA$rfscore), 
                           response = c(rep(1, 18), rep(0, 135)),
                           km_model = km_noA_beta, knm_model = knm_noA_betamix)
bm_cutoff

A_bm <- loglr_generator_general(onebarrel_example_A$rfscore, km_model = list(km_noA_beta), knm_model = list(knm_noA_betamix))

# 3. KM: betamix, KNM: beta========================================
mb_cutoff <- eer_extractor(c(onebarrel_example_km_noA$rfscore, onebarrel_example_knm_noA$rfscore), 
                           response = c(rep(1, 18), rep(0, 135)),
                           km_model = km_noA_betamix, knm_model = knm_noA_beta)
mb_cutoff

A_mb <- loglr_generator_general(onebarrel_example_A$rfscore, km_model = list(km_noA_betamix), knm_model = list(knm_noA_beta))

# 4. KM: beta, KNM: beta========================================
bb_cutoff <- eer_extractor(c(onebarrel_example_km_noA$rfscore, onebarrel_example_knm_noA$rfscore), 
                           response = c(rep(1, 18), rep(0, 135)),
                           km_model = km_noA_beta, knm_model = knm_noA_beta)
bb_cutoff

A_bb <- loglr_generator_general(onebarrel_example_A$rfscore, km_model = list(km_noA_beta), knm_model = list(knm_noA_beta))

# 5. Direct cutoff from the RF scores

noA_direct_roc <- roc(response = c(rep(1, 18), rep(0, 135)), 
                      predictor = c(onebarrel_example_km_noA$rfscore, onebarrel_example_knm_noA$rfscore))
plot(noA_direct_roc)

noA_direct_roc$auc

temp <- coords(noA_direct_roc, x = "all", transpose = F) 
temp <- temp[which.min(abs(temp$specificity - temp$sensitivity)),]
direct_cutoff <- c(temp[1,], AUC = noA_direct_roc$auc)
direct_cutoff


# Do predictions==========================================================================

mm_pred <- ifelse(A_mm>=mm_cutoff$threshold, yes = 1, no = 0)
bm_pred <- ifelse(A_bm>=bm_cutoff$threshold, yes = 1, no = 0)
mb_pred <- ifelse(A_mb>=mb_cutoff$threshold, yes = 1, no = 0)
bb_pred <- ifelse(A_bb>=bb_cutoff$threshold, yes = 1, no = 0)
direct_pred <- ifelse(A_direct>=direct_cutoff$threshold, yes = 1, no = 0)
A_true <- ifelse(onebarrel_example_A$groundtruth=="match", yes = 1, no = 0)

# Wrap a function to do ALL above for ALL barrels available, and report ALL results needed================================

work_flow_rf_roc_leaveAout <- function(onebarrel_km_example, onebarrel_knm_example) {
  # Leave A out
  onebarrel_example_km_noA <- onebarrel_km_example %>% filter(bullet2 != "A", bullet1 != "A")
  
  onebarrel_example_knm_noA <- onebarrel_knm_example %>% filter(bullet2 != "A", bullet1 != "A")
  
  
  # All data with A
  onebarrel_example_km_A <- onebarrel_km_example %>% filter(bullet1 == "A"|bullet2 == "A")
  onebarrel_example_km_A$groundtruth <- "match"
  
  # exclude comparison within A for KNM (for KM, this has been excluded already)
  onebarrel_example_knm_A <- onebarrel_knm_example %>% 
    filter(bullet1 == "A"|bullet2 == "A", !bullet1 == "A"&bullet2 == "A")
  
  onebarrel_example_knm_A$groundtruth <- "nonmatch"
  
  onebarrel_example_A <- bind_rows(onebarrel_example_km_A, onebarrel_example_knm_A)
  
  # Beta Fit============================================================================
  km_noA_beta <-  betamix(rfscore~1|1, data = onebarrel_example_km_noA, k = 1)
  knm_noA_beta <- betamix(rfscore~1|1, data = onebarrel_example_knm_noA, k = 1)
  
  
  # Beta Mix Fit=========================================================================
  km_noA_betamix <-  betamix(rfscore~1|1, data = onebarrel_example_km_noA, k = 2)
  knm_noA_betamix <-  betamix(rfscore~1|1, data = onebarrel_example_knm_noA, k = 2)
  
  # To get a cutoff value from the Known comparisons, i.e. noA=================================
  # A cutoff for the SLR, a cutoff for the distribution itself(with more than one possible nominal type 1 errors).
  # How did everyone make decisions with SLR
  # Shall we control the cutoff selection criteria, such as: always equate specificity and sensitivity
  # These are impportant questions to be answered
  
  # In land level, in bullet level?
  
  # we might do
  # 1. get the observed SLR under different models, get the equal error cutoff, evaluate the performance on the test bullet A
  # 2. direct compare the score, get equal error cutoff, evaluate the performance on the test bullet A?
  # Any other way to make use of the distributions? Are these actually same?
  
  
  # example, re-run this latter
  # noA_mm <- loglr_generator_general(c(onebarrel_example_km_noA$rfscore, onebarrel_example_knm_noA$rfscore), 
  #                                   km_model = list(km_noA_betamix), knm_model = list(knm_noA_betamix))
  # noA_mm_roc <- roc(response = c(rep(1, 18), rep(0, 135)), predictor = noA_mm)
  # 
  # temp <- coords(noA_mm_roc, x = "all", transpose = F) 
  # temp %>% filter(abs(specificity - sensitivity) <= 0.001)
  # threshold specificity sensitivity
  # -1.543217   0.8888889   0.8888889
  
  eer_extractor <- function(data, response, km_model, knm_model) {
    noA_mm <- loglr_generator_general(data, 
                                      km_model = list(km_model), knm_model = list(knm_model))
    noA_mm_roc <- roc(response = response, predictor = noA_mm)
    temp <- coords(noA_mm_roc, x = "all", transpose = F) 
    output <- temp[which.min(abs(temp$specificity - temp$sensitivity)),]
    c(output[1,], AUC = noA_mm_roc$auc)
  }
  
  
  
  # There are three possible combinations of the estimated models we will explore
  # (hard to imagine that betamix outperforms beta...)
  # 1. KM: betamix, KNM: betamix
  # 2. KM: beta,    KNM: betamix
  # 3. KM: betamix, KNM: beta
  # 4. KM: beta,    KNM: beta
  # And direct comparison (how to construct straight up SLR without distributional assumption? Some data driven, bootstrap thing?)
  
  # 1. KM: betamix, KNM: betamix========================================
  
  # rfscore = 0 & 1 replaced by 0.005, 0.995
  mm_cutoff <- eer_extractor(c(onebarrel_example_km_noA$rfscore, onebarrel_example_knm_noA$rfscore), 
                             response = c(rep(1, 18), rep(0, 135)),
                             km_model = km_noA_betamix, knm_model = knm_noA_betamix)
  #mm_cutoff 
  
  A_mm <- loglr_generator_general(onebarrel_example_A$rfscore, km_model = list(km_noA_betamix), knm_model = list(knm_noA_betamix))
  
  # 2. KM: beta, KNM: betamix========================================
  bm_cutoff <- eer_extractor(c(onebarrel_example_km_noA$rfscore, onebarrel_example_knm_noA$rfscore), 
                             response = c(rep(1, 18), rep(0, 135)),
                             km_model = km_noA_beta, knm_model = knm_noA_betamix)
  #bm_cutoff
  
  A_bm <- loglr_generator_general(onebarrel_example_A$rfscore, km_model = list(km_noA_beta), knm_model = list(knm_noA_betamix))
  
  # 3. KM: betamix, KNM: beta========================================
  mb_cutoff <- eer_extractor(c(onebarrel_example_km_noA$rfscore, onebarrel_example_knm_noA$rfscore), 
                             response = c(rep(1, 18), rep(0, 135)),
                             km_model = km_noA_betamix, knm_model = knm_noA_beta)
  #mb_cutoff
  
  A_mb <- loglr_generator_general(onebarrel_example_A$rfscore, km_model = list(km_noA_betamix), knm_model = list(knm_noA_beta))
  
  # 4. KM: beta, KNM: beta========================================
  bb_cutoff <- eer_extractor(c(onebarrel_example_km_noA$rfscore, onebarrel_example_knm_noA$rfscore), 
                             response = c(rep(1, 18), rep(0, 135)),
                             km_model = km_noA_beta, knm_model = knm_noA_beta)
  #bb_cutoff
  
  A_bb <- loglr_generator_general(onebarrel_example_A$rfscore, km_model = list(km_noA_beta), knm_model = list(knm_noA_beta))
  
  # 5. Direct cutoff from the RF scores
  
  noA_direct_roc <- roc(response = c(rep(1, 18), rep(0, 135)), 
                        predictor = c(onebarrel_example_km_noA$rfscore, onebarrel_example_knm_noA$rfscore))

  
  temp <- coords(noA_direct_roc, x = "all", transpose = F) 
  temp <- temp[which.min(abs(temp$specificity - temp$sensitivity)),]
  direct_cutoff <- c(temp[1,], AUC = noA_direct_roc$auc)
  #direct_cutoff
  
  
  # Do predictions==========================================================================
  
  mm_pred <- ifelse(A_mm>=mm_cutoff$threshold, yes = 1, no = 0)
  bm_pred <- ifelse(A_bm>=bm_cutoff$threshold, yes = 1, no = 0)
  mb_pred <- ifelse(A_mb>=mb_cutoff$threshold, yes = 1, no = 0)
  bb_pred <- ifelse(A_bb>=bb_cutoff$threshold, yes = 1, no = 0)
  direct_pred <- ifelse(onebarrel_example_A$rfscore>=direct_cutoff$threshold, yes = 1, no = 0)
  A_true <- ifelse(onebarrel_example_A$groundtruth=="match", yes = 1, no = 0)
  
  output <- list(mm_cutoff = mm_cutoff,
                 bm_cutoff = bm_cutoff,
                 mb_cutoff = mb_cutoff,
                 bb_cutoff = bb_cutoff,
                 direct_cutoff = direct_cutoff,
                 mm_pred = mm_pred,
                 bm_pred = bm_pred,
                 mb_pred = mb_pred,
                 bb_pred = bb_pred,
                 direct_pred = direct_pred,
                 A_true = A_true)
  
  output
  
}
work_flow_rf_roc_leaveAout(onebarrel_km_example = rf_km_full[[1]],
                           onebarrel_knm_example = rf_knm_full[[1]])

rf_roc_leaveAout <- lapply(1:441, FUN = function(i) {
  print(i)
  output <- try(work_flow_rf_roc_leaveAout(onebarrel_km_example = rf_km_full[[i]],
                             onebarrel_knm_example = rf_knm_full[[i]]))
  output
  })
#saveRDS(rf_roc_leaveAout, file = "~/some_saved_data/rf_roc_leaveAout.rds")

# Work to produce performance results==========================================================

# output <- list(mm_cutoff = mm_cutoff,
#                bm_cutoff = bm_cutoff,
#                mb_cutoff = mb_cutoff,
#                bb_cutoff = bb_cutoff,
#                direct_cutoff = direct_cutoff,
#                mm_pred = mm_pred,
#                bm_pred = bm_pred,
#                mb_pred = mb_pred,
#                bb_pred = bb_pred,
#                direct_pred = direct_pred,
#                A_true = A_true)

succeed_fail_cases <- lapply(rf_roc_leaveAout, FUN = is.list) %>% unlist()

rf_roc_leaveAout_succeed <- rf_roc_leaveAout[succeed_fail_cases]

rf_roc_leaveAout_succeed[[2]]

lapply(rf_roc_leaveAout_succeed, FUN = function(x) {
  is.na(x$mm_pred[1])|is.na(x$mb_pred[1])|is.na(x$bm_pred[1])|is.na(x$bb_pred[1])}) %>% 
  unlist() %>% sum()

example <- bind_rows(rf_roc_leaveAout_succeed)

example <- lapply(rf_roc_leaveAout_succeed, FUN = function(x) bind_rows(x[1:5]))

example <- lapply(rf_roc_leaveAout_succeed, FUN = function(x) {
 bind_rows(x[1])
}) %>% bind_rows()



rf_roc_extractor1 <- function(data) {
  prefix <- c("mm", "bm", "mb", "bb", "direct")
  output2 <- lapply(1:5, FUN = function(i){
    output1 <- lapply(data, FUN = function(x) {
      bind_rows(x[i])
    }) %>% bind_rows()
    
    colnames(output1) <- paste(prefix[i], colnames(output1), sep = "_")
    output1
  }) %>% bind_cols()
 output2 
}

rf_roc_result1 <- rf_roc_extractor1(rf_roc_leaveAout_succeed)
#saveRDS(rf_roc_result1, "~/some_saved_data/rf_roc_result1.rds")
rf_roc_result1 %>% colMeans() %>% t %>% as.data.frame() %>% select(ends_with("AUC"))%>% knitr::kable()
rf_roc_result1 %>% colMeans() %>% t %>% as.data.frame() %>% select(ends_with("sensitivity"))%>% knitr::kable()

rf_roc_result1_more <- rf_roc_result1 %>% colMeans() %>% matrix(nrow = 5, byrow = T)
rf_roc_result1_more
colnames(rf_roc_result1_more) <- c("Threshold", "Specificity", "Sensitivity", "AUC")
rownames(rf_roc_result1_more) <- c("mm", "bm", "mb", "bb", "direct")

rf_roc_extractor2 <- function(data) {
  output1 <- lapply(6:10, FUN = function(i) {
    sum_table <- matrix(rep(0,4),nrow = 2, ncol = 2)
    example <- lapply(rf_roc_leaveAout_succeed, FUN = function(x) {
      sum_table <<- sum_table + table(x[c(i, 11)])
    })
    sum_table
  })

  output2 <- lapply(output1, FUN = function(x) {
    test_sensitivity <- x[2, 2]/(x[2, 2] + x[1, 2])
    test_specificity <- x[1, 1]/(x[1, 1] + x[2, 1])
    
    output3 <- list(classtable = x, 
                    test_sensitivity = test_sensitivity, 
                    test_specificity = test_specificity)
  })
  
  names(output2) <- c("mm", "bm", "mb", "bb", "direct")
  
  output2
}

rf_roc_result2 <- rf_roc_extractor2(rf_roc_leaveAout_succeed)
#saveRDS(rf_roc_result2, file = "~/some_saved_data/rf_roc_result2.rds")

rf_roc_result2 %>% lapply(FUN = function(x) {x$test_sensitivity})
rf_roc_result2 %>% lapply(FUN = function(x) {x$test_specificity})

# How about variation (training and testing)
par(mfrow = c(2, 3))
boxplot(rf_roc_result1$mm_sensitivity, main = "mm_train_sensitivity")
boxplot(rf_roc_result1$bm_sensitivity, main = "bm_train_sensitivity")
boxplot(rf_roc_result1$mb_sensitivity, main = "mb_train_sensitivity")
boxplot(rf_roc_result1$bb_sensitivity, main = "bb_train_sensitivity")
boxplot(rf_roc_result1$direct_sensitivity, main = "direct_train_sensitivity")


par(mfrow = c(2, 3))
boxplot(rf_roc_result1$mm_AUC, main = "mm_train_AUC")
boxplot(rf_roc_result1$bm_AUC, main = "bm_train_AUC")
boxplot(rf_roc_result1$mb_AUC, main = "mb_train_AUC")
boxplot(rf_roc_result1$bb_AUC, main = "bb_train_AUC")
boxplot(rf_roc_result1$direct_AUC, main = "direct_train_AUC")


# extract test sensitivity and specificity 
rf_roc_extractor3 <- function(data) {
  output1 <- lapply(6:10, FUN = function(i) {
    output2 <- lapply(rf_roc_leaveAout_succeed, FUN = function(x) {
      output3 <- table(x[c(i, 11)])
      test_sensitivity <- output3[2, 2]/(output3[2, 2] + output3[1, 2])
      test_specificity <- output3[1, 1]/(output3[1, 1] + output3[2, 1])
      output4 <- list(test_sensitivity = test_sensitivity, 
                      test_specificity = test_specificity)
      output4
      
    })
    
    output2
  })
  
  names(output1) <- c("mm", "bm", "mb", "bb", "direct")
  
  output1
}

rf_roc_result3 <- rf_roc_extractor3(rf_roc_leaveAout_succeed)
rf_roc_result3_more <- lapply(rf_roc_result3, FUN = bind_rows)
rf_roc_result3_more2 <- lapply(rf_roc_result3_more, FUN = colMeans)

par(mfrow = c(2, 3))
boxplot(rf_roc_result3_more$mm$test_sensitivity, main = "mm_test_sensitivity")
boxplot(rf_roc_result3_more$bm$test_sensitivity, main = "bm_test_sensitivity")
boxplot(rf_roc_result3_more$mb$test_sensitivity, main = "mb_test_sensitivity")
boxplot(rf_roc_result3_more$bb$test_sensitivity, main = "bb_test_sensitivity")
boxplot(rf_roc_result3_more$direct$test_sensitivity, main = "direct_test_sensitivity")

par(mfrow = c(2, 3))
boxplot(rf_roc_result3_more$mm$test_specificity, main = "mm_test_specificity")
boxplot(rf_roc_result3_more$bm$test_specificity, main = "bm_test_specificity")
boxplot(rf_roc_result3_more$mb$test_specificity, main = "mb_test_specificity")
boxplot(rf_roc_result3_more$bb$test_specificity, main = "bb_test_specificity")
boxplot(rf_roc_result3_more$direct$test_specificity, main = "direct_test_specificity")



























