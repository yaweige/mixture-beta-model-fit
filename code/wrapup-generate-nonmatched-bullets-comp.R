library(tidyverse)
library(bulletxtrctr)
# input data==================================================================
# if not stored somewhere, just randomly generate the input: database_index <- sample(1:626, 100)
# randomly selected locally
excluded_index <- c(116, 611, 186, 529, 605, 225, 602,
                    246, 247, 250, 571, 608,
                    477, 478, 479, 480)

database_index <- readRDS("./data/database_index.rds")
database_index <- setdiff(database_index, excluded_index)

all_used_index <- setdiff(1:626, excluded_index)
fau_test_index <- setdiff(all_used_index, database_index)

# questioned bullets (ABCD) randomly selected locally
questioned_bullet_index <- readRDS("./data/questioned_bullet_index.rds")
# make sure the length matches. Since the some fau_test_index are eliminated, use first length(fau_test_index) bullet index

# randomization (used to get pairs ok KNM bullets pairs)=================================================
get_rearranged_index <- function(fau_test_index) {
  x <- 0
  while(x<1000) {
    x <- x + 1
    if (x==1000) {print("failed to find a good one"); break}
    rearranged_index <- sample(1:length(fau_test_index), length(fau_test_index))
    if (any(rearranged_index - 1:length(fau_test_index) == 0)) next
    # the following condition seem to be able to cover the two cases we want to avoids
    if (!any(rearranged_index - (1:length(fau_test_index))[order(rearranged_index)] == 0)) {print("find it!"); break}
  }
  
  if(any(rearranged_index - 1:length(fau_test_index) == 0)) stop("not good")
  if(any(rearranged_index - (1:length(fau_test_index))[order(rearranged_index)] == 0)) stop("not good2")
  
  rearranged_index
}

rearranged_index <- get_rearranged_index(fau_test_index)

# functions to produce test_questioned_pairs_data===============
# from now, we only use the test_question_pairs_data

get_test_questionded_pairs_data <- function(fau_test_index, 
                                            questioned_bullet_index,
                                            rearranged_index){
  
  if (length(fau_test_index) != length(questioned_bullet_index)) stop("input not same length")
  if (length(fau_test_index) != length(rearranged_index)) stop("input not same length2")
  
  questioned_bullet <- paste(fau_test_index, questioned_bullet_index)
  
  test_3_bullet <- lapply(questioned_bullet_index, FUN = function(x) {
    output <- setdiff(LETTERS[1:4], x)
    output
  }) 
  
  test_3_bullet_compact <-  test_3_bullet %>% lapply(FUN = function(x) {paste0(x, collapse = "-")}) %>% unlist()
  test_quesetioned_pairs <- paste(questioned_bullet,  
                                  paste(fau_test_index[rearranged_index], 
                                        test_3_bullet_compact[rearranged_index]), 
                                  sep = " & ")
  
  test_quesetioned_pairs_data <- data.frame(FAU1 = test_quesetioned_pairs %>% str_extract("[0-9]+(?= )"),
                                            bullet1 = test_quesetioned_pairs %>% str_extract("[A-D]+(?= &)"),
                                            FAU2 = test_quesetioned_pairs %>% str_extract("(?<=& )[0-9]+"),
                                            bullet2_1 = test_quesetioned_pairs %>% str_extract("(?<= )[A-D](?=-)"),
                                            bullet2_2 = test_quesetioned_pairs %>% str_extract("(?<=-)[A-D](?=-)"),
                                            bullet2_3 = test_quesetioned_pairs %>% str_extract("[A-D]$"))
  test_quesetioned_pairs_data
}

test_quesetioned_pairs_data <- get_test_questionded_pairs_data(fau_test_index, 
                                                               questioned_bullet_index[1:length(fau_test_index)],
                                                               rearranged_index)

write.csv(test_quesetioned_pairs_data, "./data/test_quesetioned_pairs_data.csv", row.names = FALSE)

# groove data preparation=============================================================

grooves_manual_lapd1 <- read_csv("./data/grooves-manual-lapd1.csv")
# need to make use of the crosscuts manually selected for some problematic lands. And also groove locations
grooves_manual_lapd_problematic_lands <- read_csv("./data/grooves_manual_lapd_problematic_lands.csv")

prepare_grooves <- function(grooves_manual_lapd1, 
                            grooves_manual_lapd_problematic_lands) {
  # some manipulation of the ID data
  grooves_manual_lapd1 <- grooves_manual_lapd1 %>% 
    mutate(FAU = as.numeric(str_extract(scan_id, pattern = "(?<=FAU)\\d+")),
           bullet = str_extract(scan_id, pattern = "(?<=B)[ABCD]"),
           land = as.numeric(str_extract(scan_id, pattern = "(?<=L)[123456]"))) %>%
    arrange(FAU, bullet, land)
  
  grooves_manual_lapd_problematic_lands <- grooves_manual_lapd_problematic_lands %>% 
    mutate(FAU = as.numeric(str_extract(scan_id, pattern = "(?<=FAU)\\d+")),
           bullet = str_extract(scan_id, pattern = "(?<=B)[ABCD]"),
           land = as.numeric(str_extract(scan_id, pattern = "(?<=L)[123456]"))) %>%
    arrange(FAU, bullet, land)
  
  na_grooves_manual_lapd1 <- grooves_manual_lapd1 %>% filter(is.na(crosscut))
  
  filled_na_data <- na_grooves_manual_lapd1 %>% select(-crosscut, -groove_left_manual, -groove_right_manual)  %>% 
    left_join(grooves_manual_lapd_problematic_lands %>% 
                select(FAU, bullet, land, crosscut, groove_left_manual, groove_right_manual))
  
  # 16 out 17 are filled in, the last one (FAU 416, bullet A, land 1) is manually checked and filled in
  # FAU 416, bullet A, land 1: (crosscut = 120, left groove = 337.6, right groove = 2153)
  if(filled_na_data[14,]$FAU != 416) stop("unexpected error, fau not 416")
  if(filled_na_data[14,]$bullet != "A") stop("unexpected error1")
  if(filled_na_data[14,]$land != 1) stop("unexpected error2")
  filled_na_data[14,]$crosscut <- 120
  filled_na_data[14,]$groove_left_manual <- 337.6
  filled_na_data[14,]$groove_right_manual <- 2153
  
  filled_na_data <- filled_na_data %>% select(study, source, scan_id, crosscut, manual_code,
                                              manual_rep, groove_left_manual, groove_right_manual,
                                              FAU, bullet, land)
  
  grooves_manual_lapd1[is.na(grooves_manual_lapd1$crosscut),] <- filled_na_data
  
  grooves_manual_lapd1 <- grooves_manual_lapd1 %>% 
    filter(!FAU %in% fau_not_used) %>%
    ungroup() %>%
    arrange(FAU, bullet, land)
  
  grooves_manual_lapd1
}

grooves_manual_lapd1 <- prepare_grooves(grooves_manual_lapd1, 
                                        grooves_manual_lapd_problematic_lands)

# the function to produce the comp using the test_quesetioned_pairs_data and grooves==========================
# replace the groove locations in the grooves_template to make sure we only changed the IDs not the format
# grooves_template <- readRDS("~/lapd_manual_grooves/lapd-grooves-FAU-1.rda")
# this should be the form wanted
grooves_template <- list(list(groove = c(left = 100, right = 1900)),
                         list(groove = c(left = 100, right = 1900)),
                         list(groove = c(left = 100, right = 1900)),
                         list(groove = c(left = 100, right = 1900)),
                         list(groove = c(left = 100, right = 1900)),
                         list(groove = c(left = 100, right = 1900)),
                         list(groove = c(left = 100, right = 1900)),
                         list(groove = c(left = 100, right = 1900)),
                         list(groove = c(left = 100, right = 1900)),
                         list(groove = c(left = 100, right = 1900)),
                         list(groove = c(left = 100, right = 1900)),
                         list(groove = c(left = 100, right = 1900)),
                         list(groove = c(left = 100, right = 1900)),
                         list(groove = c(left = 100, right = 1900)),
                         list(groove = c(left = 100, right = 1900)),
                         list(groove = c(left = 100, right = 1900)),
                         list(groove = c(left = 100, right = 1900)),
                         list(groove = c(left = 100, right = 1900)),
                         list(groove = c(left = 100, right = 1900)),
                         list(groove = c(left = 100, right = 1900)),
                         list(groove = c(left = 100, right = 1900)),
                         list(groove = c(left = 100, right = 1900)),
                         list(groove = c(left = 100, right = 1900)),
                         list(groove = c(left = 100, right = 1900)))

work_flow_manual_grooves3 <- function(testFAU, 
                                      testFAUrmBullet,
                                      questionedFAU, 
                                      questionedBullet,
                                      savefolder) {
  print(questionedFAU)
  path <- paste("/media/Raven/LAPD/FAU", testFAU)
  bullets <- read_dir(path)
  
  bullets <- bullets %>% 
    mutate(newbullet = toupper(bullet),
           newland = str_extract(land, pattern = "(?<=Land )[123456]")) %>%
    arrange(newbullet, newland) %>%
    select(-newbullet, -newland)
  
  # browser()
  # Add in: questioned bullet==============================================
  pathQ <- paste0("/media/Raven/LAPD/FAU ", questionedFAU, "/Bullet ", questionedBullet)
  questioned <- read_bullet(pathQ)
  
  questioned <- questioned %>%
    mutate(bullet = paste("Bullet",questionedBullet),
           land = str_extract(source, pattern = "(?<=Bullet [A-D]/).+"))
  
  # questioned <- questioned %>%
  #   mutate(bullet = "Q",
  #          land = str_extract(source, pattern = "(?<=Land )[1-6]"))
  
  row_replace <- bullets$bullet %in% paste("Bullet", testFAUrmBullet)
  bullets[row_replace, "bullet"] <- questioned$bullet
  bullets[row_replace, "land"] <- questioned$land
  
  # the commented doesn't work
  # bullets[row_replace, "x3p"] <- questioned$x3p
  row_replace_start <- min(which(row_replace))
  row_replace_end <- max(which(row_replace))
  last_index <- 0
  if (row_replace_end != 24) {
    last_index <- (row_replace_end + 1):24
  }
  bullets$x3p <- c(bullets$x3p[0:(row_replace_start - 1)], 
                   questioned$x3p, 
                   bullets$x3p[last_index])
  
  bullets[row_replace, "path"] <- paste0("/media/Raven/LAPD/FAU ", questionedFAU)
  bullets[row_replace, "files"] <- paste0("Bullet ", questionedBullet, "/", questioned$land)
  
  # end====================================================================
  # new cc section=========================================================
  
  # find the corresponding rows in the manual groove data
  my_cc_groove_data_test <- grooves_manual_lapd1 %>% filter(FAU == testFAU)
  my_cc_groove_data_questioned <- grooves_manual_lapd1 %>% filter(FAU == questionedFAU,
                                                                  bullet == questionedBullet)
  
  # cross section for all lands
  bullets <- bullets %>% mutate(crosscut = my_cc_groove_data_test$crosscut)
  bullets[row_replace, "crosscut"] <- my_cc_groove_data_questioned$crosscut
  
  # end new cc section=====================================================
  
  bullets <- bullets %>% mutate(ccdata = purrr:::map2(.x = x3p, .y = crosscut, .f = x3p_crosscut))
  
  # read in manual grooves
  # path_manual_grooves <- paste0("~/lapd_manual_grooves/lapd-grooves-FAU-", testFAU, ".rda")
  # bullets <- bullets %>% mutate(grooves = readRDS(path_manual_grooves))
  
  # new for the particular function to make use of manual groove IDs
  manual_groove_data <- grooves_manual_lapd1 %>% filter(FAU == testFAU)
  manual_groove_data[row_replace, ] <- grooves_manual_lapd1 %>% filter(FAU == questionedFAU, bullet == questionedBullet)
  
  
  if(nrow(bullets) != 24) {stop("not 24 lands? nrow of bullets is not 24")}
  if(nrow(manual_groove_data) != 24) {stop("not 24 lands? nrow of manual groove IDs is not 24")}
  
  grooves <- lapply(1:nrow(bullets), FUN = function(i){
    grooves_template[[i]]$groove[1] <- as.numeric(manual_groove_data[i, "groove_left_manual"])
    grooves_template[[i]]$groove[2] <- as.numeric(manual_groove_data[i, "groove_right_manual"])
    
    return(grooves_template[[i]])
  })
  
  bullets <- bullets %>% mutate(grooves = grooves)
  
  # get signiture
  bullets <- bullets %>% mutate(sigs = purrr::map2(.x = ccdata, .y = grooves, 
                                                   .f = function(x, y) {
                                                     cc_get_signature(ccdata = x, grooves = y, span1 = 0.75, span2 = 0.03)
                                                   }))
  
  # compare four lands
  bullets$bulletland <- paste0(bullets$bullet, "-", bullets$land)
  #bullets$bulletland[row_replace] <- paste0(bullets$bullet, " - ", "Land ", bullets$land)[row_replace]
  lands <- unique(bullets$bulletland)
  
  comparisons <- data.frame(expand.grid(land1 = lands, land2 = lands),
                            stringsAsFactors = FALSE)
  
  comparisons <- comparisons %>% mutate(aligned = purrr::map2(.x = land1, .y = land2, .f = function(xx, yy) {
    land1 <- bullets$sigs[bullets$bulletland == xx][[1]]
    land2 <- bullets$sigs[bullets$bulletland == yy][[1]]
    land1$bullet <- "first-land"
    land2$bullet <- "second-land"
    
    sig_align(land1$sig, land2$sig)
  }))
  
  comparisons <- comparisons %>% mutate(
    striae = aligned %>% purrr::map(.f = sig_cms_max, span = 75),
    bullet1 = str_extract(land1, "Bullet ([A-Z]) - Land ([1-6])") %>% str_replace("Bullet ([A-Z]) - Land ([1-6])", "\\1"), 
    bullet2 = str_extract(land2, "Bullet ([A-Z]) - Land ([1-6])") %>% str_replace("Bullet ([A-Z]) - Land ([1-6])", "\\1"), 
    landA = str_extract(land1, "Bullet ([A-Z]) - Land ([1-6])") %>% str_replace("Bullet ([A-Z]) - Land ([1-6])", "\\2"), 
    landB = str_extract(land2, "Bullet ([A-Z]) - Land ([1-6])") %>% str_replace("Bullet ([A-Z]) - Land ([1-6])", "\\2"))
  
  
  comparisons <- comparisons %>% mutate(features = purrr::map2(.x = aligned, .y = striae, 
                                                               .f = extract_features_all, resolution = 0.645), 
                                        legacy_features = purrr::map(striae, 
                                                                     extract_features_all_legacy, resolution = 0.645))
  
  p1 <- paste(testFAU, paste0(setdiff(LETTERS[1:4], testFAUrmBullet), collapse = ""), sep = "-")
  # this seq should be sep!!
  p2 <- paste(questionedFAU, questionedBullet, seq = "-")
  
  #saveRDS(comparisons, file=sprintf("./data/KNM-bullets-comp/comp-KNM-FAU-%s-vs-FAU-%s.rds", p2, p1))
  savefiles <- paste0(savefolder, "comp-KNM-FAU-%s-vs-FAU-%s.rds")
  saveRDS(comparisons, file=sprintf(savefiles, p2, p1))
  #print("test!")
}


# run the process==========================================================================================

# Have a closer look at FAU 3, 4, 5, 43 which report errors saying "crosscut data must have >0 rows"
# This is wired since the groove data we have rely on those cross cut, but works well.

# FAU 85 doesn't exist? it did exist in July...
test_quesetioned_pairs_data <- read.csv("./data/test_quesetioned_pairs_data.csv")

test_quesetioned_pairs_data <- test_quesetioned_pairs_data %>%
  mutate(temp = paste(bullet2_1, bullet2_2, bullet2_3)) %>%
  mutate(testFAUrmBullet = lapply(temp, FUN = function(x) {
    output <- setdiff(LETTERS[1:4], str_extract_all(x, pattern = "[A-D]") %>% unlist())
    output
    }) %>% unlist()) %>%
  select(-temp)

Map(f = function(a, b, c, d) {try(work_flow_manual_grooves3(testFAU = a,
                                                            testFAUrmBullet = b,
                                                            questionedFAU = c,
                                                            questionedBullet = d,
                                                            savefolder = "./data/KNM-bullets-comp/"))},
    test_quesetioned_pairs_data$FAU2,
    test_quesetioned_pairs_data$testFAUrmBullet,
    test_quesetioned_pairs_data$FAU1,
    test_quesetioned_pairs_data$bullet1)

# Map(f = function(a, b, c, d) {try(work_flow_manual_grooves3(testFAU = a, testFAUrmBullet = b, questionedFAU = c, questionedBullet = d))},
#     fau_test_index_rearranged,
#     questioned_bullet_rearranged,
#     fau_test_index,
#     questioned_bullet_index)







































