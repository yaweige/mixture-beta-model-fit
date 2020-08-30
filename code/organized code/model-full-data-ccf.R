library(ggplot2)
theme_set(theme_bw())
# note that, to avoid the (small if any) variation of the model fit each time,
# I stored the models separately, we will base on that
# read in and prepare the data =========
# use path like: "/media/Raven/LAPD-processing/data/ccf/ccf_km_full.rds"
ccf_km_full <- readRDS("./organized data/ccf/ccf_km_full.rds")
ccf_knm_full <- readRDS("./organized data/ccf/ccf_knm_full.rds")

# data processing
ccf_km_full_data <- bind_rows(ccf_km_full)
ccf_knm_full_data <- bind_rows(ccf_knm_full)

# note, there are 1's in the knm, which very likely are misscan (so, there are also those scans in km)
range(ccf_km_full_data$ccf)
range(ccf_knm_full_data$ccf)
boxplot(ccf_km_full_data$ccf)
boxplot(ccf_knm_full_data$ccf)

dim(ccf_km_full_data)
dim(ccf_knm_full_data)

ccf_knm_full_data %>% filter(ccf >= 0.90)
ccf_km_full_data %>% filter(ccf >= 0.95)
# only 29 knms out of 106080 knms are greater than 0.9, remove those

ccf_knm_full_data <- ccf_knm_full_data %>% filter(ccf < 0.90)
ccf_knm_full_data %>% dim()

# function: mine_betamix for fitting the data==============================
#par: starting values provided, different by each k
#par = c(a, b) when k = 1,
#par = c(p1, a1, b1, a2, b2), when k = 2
#par = c(p1, p2, a1, b1, a2, b2, a3, b3), when k =3
#k: number of components 1, 2 or 3.
#data: numeric vector
#control: provided to "optim"
#if the function produce warnings of NaN, that's not a problem, since I don't set hard bounds for parameters
#the "optim" function will handle that well
#suggestion: control mean increases from component 1 to 3, since
#in the sense of exchanging p and components, beta mix with two components have two optimization points,
#which could cost some difficulty in estimation. (but actually not much in the two components cases).
#However, either of these is good enough (are actually identical) as long as we are able to find one.
#This fact together with the unspecified bounds of parameters, will have some more obvious effect in
#three components cases, however, for the same reason, we should not worry too much.
#If we have components more than 3, we should restrict the parameter space explicitly.

mine_betamix <- function(data, par, k, control = NULL, hessian = FALSE){
  beta_likelihood <- function(par, data) {
    a1 <- par[1]
    b1 <- par[2]
    
    beta_1 <- dbeta(data, shape1 = a1, shape2 = b1)
    
    output <- -sum(log(beta_1))
  }
  
  betamix_likelihood <- function(par = c(0.5, 1, 1, 1, 1), data) {
    p1 <- par[1]
    a1 <- par[2]
    b1 <- par[3]
    a2 <- par[4]
    b2 <- par[5]
    
    if (p1<=0 | p1>=1 | a1<= 0 | b1 <= 0| a2<=0| b2<=0){
      output <- Inf
    } else if ((a1/(a1 + b1))>=(a2/(a2 + b2))) {
      output <- Inf
    } else {
      beta_1 <- dbeta(data, shape1 = a1, shape2 = b1)
      beta_2 <- dbeta(data, shape1 = a2, shape2 = b2)
      output <- -sum(log(p1*beta_1 + (1-p1)*beta_2))
    }
    
    output
  }
  
  beta_c3_likelihood <- function(par, data) {
    p1 <- par[1]
    p2 <- par[2]
    a1 <- par[3]
    b1 <- par[4]
    a2 <- par[5]
    b2 <- par[6]
    a3 <- par[7]
    b3 <- par[8]
    
    beta_1 <- dbeta(data, shape1 = a1, shape2 = b1)
    beta_2 <- dbeta(data, shape1 = a2, shape2 = b2)
    beta_3 <- dbeta(data, shape1 = a3, shape2 = b3)
    
    output <- -sum(log(p1*beta_1 + p2*beta_2 + (1 - p1 - p2)*beta_3))
  }
  
  if (k == 1) {
    output <- optim(par = par, fn = beta_likelihood, 
                    data = data,
                    control = control,
                    hessian = hessian)
  } else if (k == 2) {
    output <- optim(par = par, fn = betamix_likelihood, 
                    data = data,
                    control = control,
                    hessian = hessian)
  } else if (k == 3) {
    output <- optim(par = par, fn = beta_c3_likelihood, 
                    data = data,
                    control = control,
                    hessian = hessian)
  }
  
  output$value <- -output$value
  names(output)[which(names(output) == "value")] <- "logLik"
  output
}

# single beta============================
ccf_km_c1 <- mine_betamix(data = ccf_km_full_data$ccf, 
                          par = c(4.16, 2.38),
                          k = 1)
# be sure to check the convergence
ccf_km_c1$convergence

ccf_knm_c1 <- mine_betamix(data = ccf_knm_full_data$ccf, 
                           par = c(6.457, 9.54),
                           k = 1)
ccf_knm_c1$convergence

# two-component beta=========================
ccf_km_c2 <- mine_betamix(data = ccf_km_full_data$ccf, 
                          par = c(0.3, 5, 5, 5, 2),
                          k = 2)

ccf_km_c2$convergence

# starting values borrowed from betamix fit result, however, it doesn't need to be so informative,
# par = c(0.5, 5, 10, 5, 5) works well too
ccf_knm_c2 <- mine_betamix(data = ccf_knm_full_data$ccf, 
                           par = c(0.634, 11.5, 21, 6.8, 7.2),
                           k = 2,
                           control = list(maxit = 1000))
ccf_knm_c2$convergence

# three-component beta========================

ccf_km_c3 <- mine_betamix(data = ccf_km_full_data$ccf,
                          par = c(0.548, 0.2, 7.5, 4.19, 10.69, 17, 26, 6),
                          k = 3,
                          control = list(maxit = 3000))
ccf_km_c3$convergence

ccf_knm_c3 <- mine_betamix(data = ccf_knm_full_data$ccf,
                          par = c(0.2, 0.2, 10, 20, 10, 15, 10, 10),
                          k = 3,
                          control = list(maxit = 1000))
ccf_knm_c3$convergence

# plot and see===================================
# one component
myfit3_simulation <- data.frame(ccf = (1:99)/100, 
                                y = dbeta((1:99)/100, shape1 = 4.1447, shape2 = 2.384))

myfit4_simulation <- data.frame(ccf = (1:99)/100, 
                                y = dbeta((1:99)/100, shape1 = 6.38, shape2 = 9.42))

myfit3n4_simulation <- bind_rows(myfit3_simulation, myfit4_simulation, .id = "Class")

myfit3n4_simulation %>%
  ggplot(aes(color = Class, fill = Class)) + 
  geom_line(aes(x = ccf, y = y)) + 
  geom_histogram(aes(x = ccf, y = ..density..), bins = 20, color = "gray99",
                 position = "identity", alpha = 0.5, 
                 data = bind_rows(ccf_km_full_data, ccf_knm_full_data, .id = "Class")) +
  scale_color_discrete(breaks = c("1", "2"), label = c("KM", "KNM")) + 
  scale_fill_discrete(breaks = c("1", "2"), label = c("KM", "KNM")) +
  ylab("density")

# two components
myfit1_simulation <- data.frame(ccf = (1:99)/100, 
                                y1 = dbeta((1:99)/100, shape1 = 7.366, shape2 = 8.450),
                                y2 = dbeta((1:99)/100, shape1 = 14.026, shape2 = 4.451))

myfit1_simulation <- myfit1_simulation %>% mutate(y = 0.41938*y1 + (1-0.41938)*y2)

myfit2_simulation <- data.frame(ccf = (1:99)/100, 
                                y1 = dbeta((1:99)/100, shape1 = 11.0696, shape2 = 19.8388),
                                y2 = dbeta((1:99)/100, shape1 = 6.5794, shape2 = 6.745))

myfit2_simulation <- myfit2_simulation %>% mutate(y = 0.6737*y1 + (1-0.6737)*y2)

myfit1n2_simulation <- bind_rows(myfit1_simulation, myfit2_simulation, .id = "Class")

myfit1n2_simulation %>%
  ggplot(aes(color = Class, fill = Class)) + 
  geom_line(aes(x = ccf, y = y)) + 
  geom_histogram(aes(x = ccf, y = ..density..), bins = 20, color = "gray99",
                 position = "identity", alpha = 0.5, 
                 data = bind_rows(ccf_km_full_data, ccf_knm_full_data, .id = "Class")) +
  scale_color_discrete(breaks = c("1", "2"), label = c("KM", "KNM")) + 
  scale_fill_discrete(breaks = c("1", "2"), label = c("KM", "KNM")) +
  ylab("density")


# likelihood ratio test===========================================
2*(ccf_km_c2$logLik - ccf_km_c1$logLik)
1-pchisq(2*(ccf_km_c2$logLik - ccf_km_c1$logLik), df = 3)

2*(ccf_km_c3$logLik - ccf_km_c2$logLik)
1-pchisq(2*(ccf_km_c3$logLik - ccf_km_c2$logLik), df = 3)

2*(ccf_knm_c2$logLik - ccf_knm_c1$logLik)
1-pchisq(2*(ccf_knm_c2$logLik - ccf_knm_c1$logLik), df = 3)

2*(ccf_knm_c3$logLik - ccf_knm_c2$logLik)
1-pchisq(2*(ccf_knm_c3$logLik - ccf_knm_c2$logLik), df = 3)
# BIC=============================================
dim(ccf_km_full_data)
# two component beta v.s. single beta
# two component beta v.s. three component beta
# the three component beta just slightly better (about 0.0032, 0.3%), 
# we can still choose two component beta considering estimation cost (big for three-component)
# and the fact of asymptotic result of BIC may not fully achieved
bic1 <- 2*log(15912)- 2*ccf_km_c1$logLik
bic2 <- 5*log(15912)- 2*ccf_km_c2$logLik
bic3 <- 8*log(15912)- 2*ccf_km_c3$logLik
bic1
bic2
bic3

dim(ccf_knm_full_data)
# two component beta v.s. single beta
# two component beta v.s. three component beta (similar reason as above, we choose two component one)
# bic6 greater than bic5 by 0.00017, 0.017%
# also the third component has only (prior) probability of 0.04 while the first two are almost the same
bic4 <- 2*log(106051)- 2*ccf_knm_c1$logLik
bic5 <- 5*log(106051)- 2*ccf_knm_c1$logLik
bic6 <- 8*log(106051)- 2*ccf_knm_c1$logLik
bic4
bic5
bic6

# I don't know the formal proof of BIC, but I think it doesn't account the data size effect of these cases
# since the log function goes to a flat curve in such cases, 
# and doesn't make much sense when the data size increases if the data large enough, 
# compared to the magnitude of Loglik
log(15912)
log(106051)
# increased by 2

# while logLik, increased by a multiplicative factor more than 10!
ccf_km_c2$logLik
ccf_knm_c2$logLik



