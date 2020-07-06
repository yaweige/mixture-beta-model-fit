library(dplyr)
library(x3ptools)
library(bulletxtrctr)
library(ggplot2)
library(purrr)
library(stringr)
library(readr)

# This is a modified function from ccf_compare to make use of the .csv full groove IDs which are different from those
# manually generated from Yawei previously

# need to properly label the groove data to be able to directly use
grooves_manual_lapd1 <- read_csv("grooves-manual-lapd1.csv")

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

# subset the goove locations in thie grooves_template to make sure we only changed the IDs not the format
grooves_template <- readRDS("~/lapd_manual_grooves/lapd-grooves-FAU-1.rda")
example <- readRDS("~/lapd_manual_grooves/lapd-grooves-FAU-1.rda")

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
  
  saveRDS(comparisons, file=sprintf("~/comp_full_manual_grooves_notMine_updated/comp-full-manual-groove-FAU-%s.rds", FAUno))
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
  
  ggsave(filename = sprintf("results-full-manual-grooves-FAU-%s.png", FAUno), path = "~/Pictures/LAPD_full_manual_grooves")
}

# Have a closer look at FAU 3, 4, 5, 43 which report errors saying "crosscut data must have >0 rows"
# This is wired since the groove data we have rely on those cross cut, but works well.

lapply(setdiff(381:626, c(116, 186, 529, 605)), FUN = function(x) {try(work_flow_manual_grooves2(x))})


###Attach the session infor to when run the above code (2020.06.20)===============================================================


# R version 3.6.3 (2020-02-29)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 18.04.4 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/openblas/libblas.so.3
# LAPACK: /usr/lib/x86_64-linux-gnu/libopenblasp-r0.2.20.so
# 
# Random number generation:
#   RNG:     Mersenne-Twister 
# Normal:  Inversion 
# Sample:  Rounding 
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] betareg_3.1-3       forcats_0.4.0       readr_1.3.1         tidyr_1.0.0         tibble_2.1.3        ggplot2_3.2.1      
# [7] tidyverse_1.3.0     stringr_1.4.0       shiny_1.4.0         purrr_0.3.3         bulletxtrctr_0.2.0  x3ptools_0.0.2.9000
# [13] dplyr_0.8.3         flexmix_2.3-15      lattice_0.20-38    
# 
# loaded via a namespace (and not attached):
#   [1] nlme_3.1-147            fs_1.3.1                xts_0.11-2              lubridate_1.7.4         webshot_0.5.2          
# [6] httr_1.4.1              tools_3.6.3             backports_1.1.5         R6_2.4.1                DBI_1.0.0              
# [11] lazyeval_0.2.2          colorspace_1.4-1        nnet_7.3-12             manipulateWidget_0.10.0 withr_2.1.2            
# [16] readbitmap_0.1.5        tidyselect_0.2.5        smoother_1.1            curl_4.3                compiler_3.6.3         
# [21] rvest_0.3.5             cli_2.0.1               bulletcp_1.0.0          xml2_1.2.2              sandwich_2.5-1         
# [26] scales_1.1.0            lmtest_0.9-37           mvtnorm_1.0-12          digest_0.6.23           tiff_0.1-5             
# [31] jpeg_0.1-8.1            pkgconfig_2.0.3         htmltools_0.4.0         bibtex_0.4.2.2          dbplyr_1.4.2           
# [36] fastmap_1.0.1           htmlwidgets_1.5.1       rlang_0.4.2             readxl_1.3.1            TTR_0.23-6             
# [41] rstudioapi_0.10         generics_0.0.2          zoo_1.8-7               jsonlite_1.6            crosstalk_1.0.0        
# [46] magrittr_1.5            modeltools_0.2-22       Formula_1.2-3           Rcpp_1.0.3              munsell_0.5.0          
# [51] fansi_0.4.1             grooveFinder_0.0.1      lifecycle_0.1.0         stringi_1.4.5           gbRd_0.4-11            
# [56] MASS_7.3-51.6           plyr_1.8.5              promises_1.1.0          crayon_1.3.4            miniUI_0.1.1.1         
# [61] haven_2.2.0             hms_0.5.3               locfit_1.5-9.1          zeallot_0.1.0           knitr_1.27             
# [66] pillar_1.4.3            igraph_1.2.4.2          imager_0.41.2           stats4_3.6.3            reprex_0.3.0           
# [71] glue_1.3.1              bmp_0.3                 modelr_0.1.5            png_0.1-7               vctrs_0.2.1            
# [76] httpuv_1.5.2            Rdpack_0.11-1           cellranger_1.1.0        gtable_0.3.0            assertthat_0.2.1       
# [81] xfun_0.12               mime_0.8                xtable_1.8-4            broom_0.5.2             later_1.0.0            
# [86] rgl_0.100.30 






