library(bulletxtrctr)
library(dplyr)

# To generate the ccf scores we will use=======================================

# Use land_score_km_generator2 and land_score_knm_generator2 (in helper functions)

# Run the following to get the full data
ccf_km_full <- land_score_km_generator2("~/comp_full_manual_grooves_notMine_updated/")
ccf_knm_full <- land_score_knm_generator2("~/comp_full_manual_grooves_notMine_updated/")
#
saveRDS(ccf_km_full, file = "~/some_saved_data/ccf_km_full.rds")
saveRDS(ccf_knm_full, file = "~/some_saved_data/ccf_knm_full.rds")

# we will make use of the following two functions, which are included in the helper function file============
# If there are any other implicit dependencies, see the help function file (will be posted, should be good for these code)

# this works within a barrel to extract the km land-land ccf comparisons
# except the first argument "dir", other arguments are not well defined

# read in km
land_score_km_generator2 <- function(dir, FAUno = NULL, firstno = NULL, remove = NULL) {
  files <- paste0(dir, list.files(dir))
  
  if(!is.null(firstno)){
    files <- files[1:firstno]
  }
  
  if(!is.null(remove)){
    files <- files[setdiff(1:length(files), remove)]
  }
  
  output <- lapply(files, function(x) {
    comp_manual_fau <- readRDS(x)
    features_fau <- comp_manual_fau %>% tidyr::unnest(legacy_features)
    features_slc <- features_fau %>% select(bullet1, bullet2, landA, landB, ccf, land1, land2)
    
    reference <- list(n1 = data.frame(landA = 1:6, landB = 1:6),
                      n2 = data.frame(landA = 1:6, landB = rep(1:6, times = 2)[2:7]),
                      n3 = data.frame(landA = 1:6, landB = rep(1:6, times = 2)[3:8]),
                      n4 = data.frame(landA = 1:6, landB = rep(1:6, times = 2)[4:9]),
                      n5 = data.frame(landA = 1:6, landB = rep(1:6, times = 2)[5:10]),
                      n6 = data.frame(landA = 1:6, landB = rep(1:6, times = 2)[6:11])) 
    
    features_sum <- features_slc %>% 
      group_by(bullet1, bullet2) %>%
      mutate(max_diagonal_score = which.max(compute_average_scores(as.numeric(landA), as.numeric(landB), ccf))) %>%
      ungroup() %>%
      group_by(bullet1, bullet2, landA) %>%
      filter(as.numeric(landB) == reference[[max_diagonal_score[1]]][as.numeric(landA[1]),2]) %>%
      ungroup()
    
    
    # reduced duplicated ones, but it doesn't affect the distribution
    output <- features_sum %>% 
      filter(as.numeric(factor(bullet1)) > as.numeric(factor(bullet2))) %>%
      select(ccf, bullet1, bullet2, land1, land2)
    
    output
  })
  
  output
}

#onebarrel_km_temp <- land_score_km_generator2(dir = "D:/LAPD-comp/")

# read in knm

# this works within a barrel to extract the knm land-land comparasions
land_score_knm_generator2 <- function(dir, FAUno = NULL, firstno = NULL, remove = NULL) {
  files <- paste0(dir, list.files(dir))
  
  if(!is.null(firstno)){
    files <- files[1:firstno]
  }
  
  if(!is.null(remove)){
    files <- files[setdiff(1:length(files), remove)]
  }
  
  output <- lapply(files, function(x) {
    comp_manual_fau <- readRDS(x)
    features_fau <- comp_manual_fau %>% tidyr::unnest(legacy_features)
    features_slc <- features_fau %>% select(bullet1, bullet2, landA, landB, ccf, land1, land2)
    
    reference <- list(n1 = data.frame(landA = 1:6, landB = 1:6),
                      n2 = data.frame(landA = 1:6, landB = rep(1:6, times = 2)[2:7]),
                      n3 = data.frame(landA = 1:6, landB = rep(1:6, times = 2)[3:8]),
                      n4 = data.frame(landA = 1:6, landB = rep(1:6, times = 2)[4:9]),
                      n5 = data.frame(landA = 1:6, landB = rep(1:6, times = 2)[5:10]),
                      n6 = data.frame(landA = 1:6, landB = rep(1:6, times = 2)[6:11])) 
    
    features_sum <- features_slc %>% 
      group_by(bullet1, bullet2) %>%
      mutate(max_diagonal_score = which.max(compute_average_scores(as.numeric(landA), as.numeric(landB), ccf))) %>%
      ungroup() %>%
      group_by(bullet1, bullet2, landA) %>%
      filter(as.numeric(landB) != reference[[max_diagonal_score[1]]][as.numeric(landA[1]),2]) %>%
      ungroup()
    
    
    # reduced duplicated ones, but it doesn't affect the distribution
    output <- features_sum %>% 
      filter(as.numeric(factor(bullet1)) >= as.numeric(factor(bullet2))) %>%
      filter(as.numeric(factor(bullet1)) > as.numeric(factor(bullet2))|as.numeric(factor(landA)) > as.numeric(factor(landB))) %>%
      select(ccf, bullet1, bullet2, land1, land2)
    
    output
  })
  
  output
}
