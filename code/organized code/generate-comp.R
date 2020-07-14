# The code is run on the server
# Change the save path before run the function, define the saved names, check all inputs
# Functions are not well defined since it depends on some values outside the function, 
# be sure to run all of this script.
# Some comments are from previous runs, and they are kept

library(dplyr)
library(x3ptools)
library(bulletxtrctr)
library(ggplot2)
library(purrr)
library(stringr)
library(readr)

# This is a modified function from ccf_compare to make use of the .csv full groove IDs which are different from those
# manually generated from Yawei previously

# Input============================================
# 1. need to properly label the groove data to be able to directly use
grooves_manual_lapd1 <- read_csv("/media/Raven/LAPD-processing/data/grooveID/grooves-manual-lapd1.csv")

# Or use the URL and the repo, but the code itself doesn't work for me, not sure where it's wrong
# grooves_manual_lapd1 <- read.csv("https://github.com/CSAFE-ISU/grooves-manual-id/blob/master/grooves-manual-lapd1.csv")

# 2. To replace the goove locations in the grooves_template to make sure we only changed the IDs not the format
# Exact which FAU are used here is not important, we just need the format
grooves_template <- readRDS("/media/Raven/LAPD-processing/data/example groove form/lapd-grooves-FAU-1.rda")

# 3. the LAPD x3p data on the server, used inside the function, we used:
# "/media/Raven/LAPD/"

# Output======================================
# 1. It produces results in the "/media/Raven/LAPD-processing/data/comp", unless you changed it
# 2. It can produce 4(bullets) by 4(bullets) prictures for ccf, but it turned off in default

# Run====================================
# All of them are complete data? No!
sum(complete.cases(grooves_manual_lapd1))

# some manipulation of the ID data
grooves_manual_lapd1 <- grooves_manual_lapd1 %>% 
  mutate(FAU = as.numeric(str_extract(scan_id, pattern = "(?<=FAU)\\d+")),
         bullet = str_extract(scan_id, pattern = "(?<=B)[ABCD]"),
         land = as.numeric(str_extract(scan_id, pattern = "(?<=L)[123456]"))) %>%
  arrange(FAU, bullet, land)

# Note that some FAUs has less/more than 24 lands data (FAU 116, 186, 529, 605) for various reasons
grooves_manual_lapd1 %>% group_by(FAU) %>% tally() %>% filter(n!=24)
example <- grooves_manual_lapd1 %>% filter(FAU == 529)

# example <- lapply(1:24, FUN = function(i){
#   grooves_template[[i]]$groove[1] <- as.numeric(grooves_manual_lapd1[i, "groove_left_manual"])
#   grooves_template[[i]]$groove[2] <- as.numeric(grooves_manual_lapd1[i, "groove_right_manual"])
#   
#   return(grooves_template[[i]])
# })

work_flow_manual_grooves2 <- function(FAUno) {
  print(FAUno)
  path <- paste("/media/Raven/LAPD/FAU", FAUno)
  bullets <- read_dir(path)
  # cross section for all lands
  bullets <- bullets %>% mutate(crosscut = x3p %>% purrr::map_dbl(.f = x3p_crosscut_optimize))
  
  bullets <- bullets %>% mutate(ccdata = purrr:::map2(.x = x3p, .y = crosscut, .f = x3p_crosscut))
  
  # read in manual grooves
  # path_manual_grooves <- paste0("~/lapd_manual_grooves/lapd-grooves-FAU-", FAUno, ".rda")
  # bullets <- bullets %>% mutate(grooves = readRDS(path_manual_grooves))
  
  # new for the particular function to make use of manual groove IDs
  manual_groove_data <- grooves_manual_lapd1 %>% filter(FAU == FAUno)
  
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
  saveRDS(comparisons, file=sprintf("/media/Raven/LAPD-processing/data/comp/comp-full-manual-groove-FAU-%s.rds", FAUno))
  # only if you also want to save the pictures
  # ggsave(filename = sprintf("results-full-manual-grooves-FAU-%s.png", FAUno), path = "~/Pictures/LAPD_full_manual_grooves")
}

# Have a closer look at FAU 3, 4, 5, 43 which report errors saying "crosscut data must have >0 rows"
# This is wired since the groove data we have rely on those cross cut, but works well.

# Provide the FAU number as in put, (116, 186, 529, 605) are known to have bullets different from 4.

lapply(setdiff(1:626, c(116, 186, 529, 605)), FUN = function(x) {try(work_flow_manual_grooves2(x))})

# The following is some code to double check if the missed ones really can't find crosscuts. And Yes.
# newlist <- setdiff(1:626,c(as.numeric(str_extract(list.files("~/comp_full_manual_grooves_notMine_updated/"), pattern = "[0123456789]+")),116, 186, 529, 605))
# lapply(newlist, FUN = function(x) {try(work_flow_manual_grooves2(x))})

