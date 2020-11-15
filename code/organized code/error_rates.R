library(tidyverse)

ccf_km_c2$value
ccf_km_c2$par
ccf_km_c2$counts
ccf_km_c2$convergence

# Two-component betamix error rates===================

# find the cutoff where the pdfs crossing

dbetamix2(0.5, ccf_km_c2)
dbetamix2(0.5, ccf_knm_c2)
dbetamix2(0.25, ccf_km_c2)
dbetamix2(0.25, ccf_knm_c2)

find_cross <- function(m1, m2, start = 0.4){
  a <- 1
  b <- 2
  x <- start
  while(abs(a - b)>=0.002) {
    x <- x + 0.001
    a <- dbetamix2(x, m1)
    b <- dbetamix2(x, m2)
    
    if (x > 0.7){
      break
    }
  }
  x
}

cutoff1 <- find_cross(ccf_km_c2, ccf_knm_c2, 0.4)
cutoff1
# 0.529
dbetamix2(0.529, ccf_km_c2)
dbetamix2(0.529, ccf_knm_c2)

# FPR: 0.148
1 - pbetamix2(0.529, ccf_knm_c2)
# FNR: 0.301
pbetamix2(0.529, ccf_km_c2)
# FIR: 0.212
(1 - pbetamix2(0.529, ccf_knm_c2))/(1 - pbetamix2(0.529, ccf_km_c2))
# FER: 0.354
pbetamix2(0.529, ccf_km_c2)/pbetamix2(0.529, ccf_knm_c2)

# The plot for FNR and FPR==========================================
polydata_fnr <- myfit1n2_simulation %>% 
  filter(ccf<=0.53, Class == "1") %>%
  select(ccf, y) %>%
  bind_rows(data.frame(ccf = 0.53, y = 0))
polydata_fpr <- myfit1n2_simulation %>% 
  filter(ccf>=0.53, Class == "2") %>%
  select(ccf, y) %>%
  bind_rows(data.frame(ccf = 0.53, y = 0))
polydata_fnrfpr <- bind_rows(polydata_fnr, polydata_fpr, .id = "Class")

myfit1n2_simulation %>%
  ggplot() + 
  geom_polygon(aes(x = ccf, y = y, fill = Class), data = polydata_fnrfpr, alpha = 0.5) +
  geom_vline(xintercept = 0.529, linetype = 2) +
  geom_line(aes(x = ccf, y = y, color = Class)) +  
  geom_text(aes(x = ccf, y = y, label = label), data = data.frame(ccf = 0.72, y = 3, label = "Threshold = 0.529")) + 
  geom_text(aes(x = ccf, y = y, label = label), data = data.frame(ccf = c(0.375, 0.6), y = c(0.5, 0.45), label = c("FNR", "FPR"))) + 
  scale_color_discrete(breaks = c("1", "2"), label = c("KM", "KNM")) + 
  scale_fill_discrete(breaks = c("1", "2"), label = c("KM", "KNM")) +
  ylab("density")


# The plot for FIR and FER===================================
polydata_part3 <- myfit1n2_simulation %>%
  filter(ccf<=0.53)

polydata_part3[1:53,] <- polydata_part3[1:53,] %>% arrange(desc(ccf))

polydata_part4 <- myfit1n2_simulation %>%
  filter(ccf>=0.53) 

polydata_part4[1:47,] <- polydata_part4[1:47,] %>% arrange(desc(ccf))

polydata_part34 <- bind_rows(polydata_part3 %>% select(ccf, y), polydata_part4 %>% select(ccf, y), .id = "Class")
polydata_part34 <- polydata_part34 %>% mutate(Class = as.character(as.numeric(Class) + 2))

polydata <- bind_rows(polydata_fnrfpr, polydata_part34)

polydata %>% 
  ggplot() + 
  geom_polygon(aes(x = ccf, y = y, group = Class), alpha = 0.1) +
  geom_text(aes(x = ccf, y = y, label = label), 
            data = data.frame(ccf = c(0.375, 0.6, 0.375, 0.75), y = c(0.5, 0.45, 2, 1), label = c("A1", "A2", "A3", "A4"))) + 
  geom_line(aes(x = ccf, y = y, color = Class), data = myfit1n2_simulation) +
  geom_vline(xintercept = 0.529, linetype = 2) +
  geom_text(aes(x = ccf, y = y, label = label), data = data.frame(ccf = 0.72, y = 3, label = "Threshold = 0.529")) + 
  scale_color_discrete(breaks = c("1", "2"), label = c("KM", "KNM")) + 
  ylab("density")
  



















