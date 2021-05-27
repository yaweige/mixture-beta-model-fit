library(dplyr)
library(x3ptools)
library(bulletxtrctr)
library(ggplot2)
library(purrr)
library(stringr)
library(readr)

# Input============================================
# 1. need to properly label the groove data to be able to directly use
grooves_manual_lapd1 <- read_csv("./data/grooves-manual-lapd1.csv")

# Or use the URL and the repo, but the code itself doesn't work for me, not sure where it's wrong
# grooves_manual_lapd1 <- read.csv("https://github.com/CSAFE-ISU/grooves-manual-id/blob/master/grooves-manual-lapd1.csv")

# 2. need to make use of the crosscuts manually selected for some problematic lands. And also groove locations
grooves_manual_lapd_problematic_lands <- read_csv("./data/grooves_manual_lapd_problematic_lands.csv")

# 3. To replace the goove locations in the grooves_template to make sure we only changed the IDs not the format
# Exact which FAU are used here is not important, we just need the format
grooves_template <- readRDS("/media/Raven/LAPD-processing/data/example groove form/lapd-grooves-FAU-1.rda")
#grooves_template <- readRDS("./data/lapd-grooves-FAU-1.rda")


# 4. the LAPD x3p data on the server, used inside the function, we used:
# "/media/Raven/LAPD/"

# Prepare the grooves, crosscut information=====================================
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



work_flow_manual_grooves2 <- function(FAUno, savefolder) {
  print(FAUno)
  path <- paste("/media/Raven/LAPD/FAU", FAUno)
  bullets <- read_dir(path)
  if(nrow(bullets) != 24) stop("lands not 24")
  
  bullets <- bullets %>% 
    mutate(newbullet = toupper( bullet),
           newland = str_extract(land, pattern = "(?<=Land )[123456]")) %>%
    arrange(newbullet, newland) %>%
    select(-newbullet, -newland)
  # find the corresponding rows in the manual groove data
  my_cc_groove_data <- grooves_manual_lapd1 %>% filter(FAU == FAUno)
  
  # cross section for all lands
  # bullets <- bullets %>% mutate(crosscut = x3p %>% purrr::map_dbl(.f = x3p_crosscut_optimize))
  bullets <- bullets %>% mutate(crosscut = my_cc_groove_data$crosscut)
  
  bullets <- bullets %>% mutate(ccdata = purrr:::map2(.x = x3p, .y = crosscut, .f = x3p_crosscut))
  
  left_groove <- my_cc_groove_data$groove_left_manual
  right_groove <- my_cc_groove_data$groove_right_manual
  
  grooves <- lapply(1:nrow(bullets), FUN = function(i){
    grooves_template[[i]]$groove[1] <- as.numeric(left_groove[i])
    grooves_template[[i]]$groove[2] <- as.numeric(right_groove[i])
    
    return(grooves_template[[i]])
  })
  
  bullets <- bullets %>% mutate(grooves = grooves)
  
  # get signiture
  bullets <- bullets %>% mutate(sigs = purrr::map2(.x = ccdata, .y = grooves, 
                                                   .f = function(x, y) {
                                                     cc_get_signature(ccdata = x, grooves = y, span1 = 0.75, span2 = 0.03)
                                                   }))
  
  # compare lands
  bullets$bulletland <- paste0(bullets$bullet, "-", bullets$land)
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
    bullet1 = str_extract(land1, "Bullet ([A-D]) - Land ([1-6])") %>% str_replace("Bullet ([A-D]) - Land ([1-6])", "\\1"), 
    bullet2 = str_extract(land2, "Bullet ([A-D]) - Land ([1-6])") %>% str_replace("Bullet ([A-D]) - Land ([1-6])", "\\1"), 
    landA = str_extract(land1, "Bullet ([A-D]) - Land ([1-6])") %>% str_replace("Bullet ([A-D]) - Land ([1-6])", "\\2"), 
    landB = str_extract(land2, "Bullet ([A-D]) - Land ([1-6])") %>% str_replace("Bullet ([A-D]) - Land ([1-6])", "\\2"))
  
  
  comparisons <- comparisons %>% mutate(features = purrr::map2(.x = aligned, .y = striae, 
                                                               .f = extract_features_all, resolution = 0.645), 
                                        legacy_features = purrr::map(striae, 
                                                                     extract_features_all_legacy, resolution = 0.645)) 
  
  
  features <- comparisons %>% tidyr::unnest(legacy_features)
  
  titles <- paste("Manual-FAU", FAUno)
  
  features %>% ggplot(aes(x = landA, y = landB, fill = ccf)) + 
    geom_tile() + 
    scale_fill_gradient2(low = "grey80", high = "darkorange", midpoint = 0.5) + 
    facet_grid(bullet2 ~ bullet1, labeller = "label_both") + 
    xlab("Land A") + 
    ylab("Land B") + 
    ggtitle(titles) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(aspect.ratio = 1)
  # output
  # saveRDS(comparisons, file=sprintf("./data/comp/comp-full-manual-groove-FAU-%s.rds", FAUno))
  
  savefiles <- paste0(savefolder, "comp-full-manual-groove-FAU-%s.rds")
  saveRDS(comparisons, file=sprintf(savefiles, FAUno))
  # only if you also want to save the pictures
  # ggsave(filename = sprintf("results-full-manual-grooves-FAU-%s.png", FAUno), path = "~/Pictures/LAPD_full_manual_grooves")
}


excluded_index <- c(116, 611, 186, 529, 605, 225, 602,
                    246, 247, 250, 571, 608,
                    477, 478, 479, 480)

lapply(setdiff(1:626, excluded_index), FUN = function(x) {
  try(work_flow_manual_grooves2(x, savefolder = "./data/comp/"))})



















































































