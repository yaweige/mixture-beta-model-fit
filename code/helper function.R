## beta mixture params extractor(produce Beta(alpha, beta))
# note: unlike cv_summary, this function didn't reorder the components according to their means
# use cv_summary instead if not just want to see the result
beta_params <- function(model){
  
  output <- lapply(model$flexmix@components, FUN = function(x){
    beta_mean <- plogis(x[[1]]@parameters$mean)
    beta_phi <- exp(x[[1]]@parameters$precision)
    
    alpha <- beta_mean*beta_phi
    beta <- beta_phi - alpha
    
    output <- list(alpha = alpha, beta = beta)
  })
  
  output
}

# A function to produce the cumulative probability of the betamix model given an observed value(two component models)
# example:1 - pbetamix(500:550/1000, knm_42_leaveAout)

pbetamix <- function(x, model){
  params <- beta_params(model)
  p1 <- model$flexmix@prior[1]
  p2 <- model$flexmix@prior[2]
  
  cdf <- p1*pbeta(x, shape1 = params[[1]][[1]], shape2 = params[[1]][[2]]) + 
    p2*pbeta(x, shape1 = params[[2]][[1]], shape2 = params[[2]][[2]])
  
  cdf
}

# similarly, for pdf
dbetamix <- function(x, model){
  params <- beta_params(model)
  p1 <- model$flexmix@prior[1]
  p2 <- model$flexmix@prior[2]
  
  pdf <- p1*dbeta(x, shape1 = params[[1]][[1]], shape2 = params[[1]][[2]]) + 
    p2*dbeta(x, shape1 = params[[2]][[1]], shape2 = params[[2]][[2]])
  
  pdf
}

# for random draw from the mixture beta distritbuion, 
# we use two steps procedure, first draw from each beta component, than draw from bernulli, then, combine them
rbetamix <- function(n, p1, a1, b1, a2, b2){
  c1 <- rbeta(n = n, shape1 = a1, shape2 = b1)
  c2 <- rbeta(n = n, shape1 = a2, shape2 = b2)
  p <- purrr::rbernoulli(n, p = p1)
  
  r <- ifelse(p == 1, yes = c1, no = c2)
  r
}


# A function to produce the quantile given a probabilty
# a temporal one can be used 
qbetamix <- function(x, model){
# objective function
  distance <- function(par, x, model){
    abs(pbetamix(par, model = model) - x)
  }
  
  result <- optimize(f = distance,
                     interval = c(0.01, 0.99),
                     x = x,
                     model = model)
  result$minimum
}

# summarize a list of betamixture models(reorder the components here by their means)
# example:km_bi_summary <- cv_summary(km_bi_list)
cv_summary <- function(x) {
  information_list <- lapply(x, FUN = function(x) {
    prior <- x$flexmix@prior
    params <- beta_params(x)
    mu_phi_params <- x$flexmix@components
    
    # order the components according to their mean: from small to large
    order <- order(unlist(lapply(params, FUN = function(x) x[[1]]/(x[[1]]+x[[2]]))))
    prior <- prior[order]
    params <- params[order]
    mu_phi_params <- mu_phi_params[order]
    names(params) <- NULL
    names(mu_phi_params) <- NULL
    
    # simlify the structure(this is updated on June 5/2020)
    mu_phi_params <- lapply(1:length(mu_phi_params), FUN = function(i){
      output <- mu_phi_params[[i]][[1]]@parameters
    })
    
    
    output <- list(prior=prior, params=params, order=order, mu_phi_params=mu_phi_params)
    output
  })
  information_list
}



# use one knm(two components), multiple km models(two components) to generate loglikelihood ratio
# example: loglr_generator(0.523, km_model = list(km_42_leaveAout), knm_model = list(knm_42_leaveAout))
loglr_generator <- function(testset, km_model, knm_model){
  km_cv_sum <- cv_summary(km_model)
  knm_cv_sum <- cv_summary(knm_model)
  
  loglr <- lapply(1:length(km_cv_sum), FUN = function(x){
    km_likelihood <- km_cv_sum[[x]][[1]][[1]]*dbeta(testset, shape1 = km_cv_sum[[x]][[2]][[1]][[1]],
                                                    shape2 = km_cv_sum[[x]][[2]][[1]][[2]]) +
      km_cv_sum[[x]][[1]][[2]]*dbeta(testset, shape1 = km_cv_sum[[x]][[2]][[2]][[1]],
                                     shape2 = km_cv_sum[[x]][[2]][[2]][[2]])
    knm_likelihood <- knm_cv_sum[[1]][[1]][[1]]*dbeta(testset, shape1 = knm_cv_sum[[1]][[2]][[1]][[1]],
                                                      shape2 = knm_cv_sum[[1]][[2]][[1]][[2]]) + 
      knm_cv_sum[[1]][[1]][[2]]*dbeta(testset, shape1 = knm_cv_sum[[1]][[2]][[2]][[1]],
                                      shape2 = knm_cv_sum[[1]][[2]][[2]][[2]])
    
    loglr <- log(km_likelihood/knm_likelihood)
  })
  
  loglr
}

# A function to simulate data set
# example:
# km_42_simulated <- simulate_bitamix(km_42_leaveAout)
# 
# ggplot(mapping = aes(x = ccf, y = ..density..)) + 
#   geom_line(aes(x = ccf, y = y), data = km_42_simulated, color = "green") + 
#   geom_line(aes(x = ccf, y = y1), data = km_42_simulated, color = "blue") + 
#   geom_line(aes(x = ccf, y = y2), data = km_42_simulated, color = "red")

simulate_bitamix <- function(model) {
  params <- beta_params(model)
  
  if (length(model$flexmix@prior) == 2) {
    p1 <- model$flexmix@prior[1]
    p2 <- model$flexmix@prior[2]
    # y: the betamix, y1: component 1, y2: component 2
    simulated <- data.frame(xlab = 1:99/100,
                            y1 =  dbeta(1:99/100, shape1 = params[[1]][[1]], shape2 = params[[1]][[2]]),
                            y2 =  dbeta(1:99/100, shape1 = params[[2]][[1]], shape2 = params[[2]][[2]]),
                            y = p1*dbeta(1:99/100, shape1 = params[[1]][[1]], shape2 = params[[1]][[2]]) +
                              p2*dbeta(1:99/100, shape1 = params[[2]][[1]], shape2 = params[[2]][[2]]))
  }
  
  if (length(model$flexmix@prior) == 1) {
    # y: the beta
    simulated <- data.frame(xlab = 1:99/100,
                            y = dbeta(1:99/100, shape1 = params[[1]][[1]], shape2 = params[[1]][[2]]))
  }

  simulated
}

## km land score generator========================
# The following two functions are used to extract corresponding ccf scores, 
# but note that these functions are written earlier, i.e. it may work with 4 bullets, 
# for knm, it consider cases from different barrels
# for within one barrel case, see the the 2-version of the functions
land_score_km_generator <- function(dir, FAUno = NULL, firstno = NULL, remove = NULL) {
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
    features_slc <- features_fau %>% select(bullet1, bullet2, landA, landB, ccf)
    
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
      select(ccf)
    
    output
  })
  
  output
}

reference <- list(n1 = data.frame(landA = 1:6, landB = 1:6),
                  n2 = data.frame(landA = 1:6, landB = rep(1:6, times = 2)[2:7]),
                  n3 = data.frame(landA = 1:6, landB = rep(1:6, times = 2)[3:8]),
                  n4 = data.frame(landA = 1:6, landB = rep(1:6, times = 2)[4:9]),
                  n5 = data.frame(landA = 1:6, landB = rep(1:6, times = 2)[5:10]),
                  n6 = data.frame(landA = 1:6, landB = rep(1:6, times = 2)[6:11]))

## knm land score generator========================
## Use with caution, it seems to work with special data
land_score_knm_generator <- function(dir, FAUno = NULL, firstno = NULL, remove = NULL) {
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
    features_slc <- features_fau %>% select(bullet1, bullet2, landA, landB, ccf)
    
    
    # reduced duplicated ones, but it doesn't affect the distribution
    output <- features_slc %>% 
      filter(as.numeric(factor(bullet1)) > as.numeric(factor(bullet2))) %>%
      select(ccf)
    
    output
  })
  
  output
}

# this works within a barrel to extract the knm land-land comparasions
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
    features_slc <- features_fau %>% select(bullet1, bullet2, landA, landB, ccf)
    
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
      select(ccf, bullet1, bullet2)
    
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
    features_slc <- features_fau %>% select(bullet1, bullet2, landA, landB, ccf)
    
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
      select(ccf, bullet1, bullet2)
    
    output
  })
  
  output
}
#onebarrel_km_temp <- land_score_km_generator3_RF(dir = "D:/LAPD-comp/")

# this works within a barrel to extract the km land-land RF comparisons
# read in km
land_score_km_generator3_RF <- function(dir, FAUno = NULL, firstno = NULL, remove = NULL) {
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
    features_fau$rfscore <- predict(bulletxtrctr::rtrees, newdata = features_fau, 
                                    type = "prob")[, 2]
    features_slc <- features_fau %>% select(bullet1, bullet2, landA, landB, rfscore, land1, land2)
    
    reference <- list(n1 = data.frame(landA = 1:6, landB = 1:6),
                      n2 = data.frame(landA = 1:6, landB = rep(1:6, times = 2)[2:7]),
                      n3 = data.frame(landA = 1:6, landB = rep(1:6, times = 2)[3:8]),
                      n4 = data.frame(landA = 1:6, landB = rep(1:6, times = 2)[4:9]),
                      n5 = data.frame(landA = 1:6, landB = rep(1:6, times = 2)[5:10]),
                      n6 = data.frame(landA = 1:6, landB = rep(1:6, times = 2)[6:11])) 
    
    features_sum <- features_slc %>% 
      group_by(bullet1, bullet2) %>%
      mutate(max_diagonal_score = which.max(compute_average_scores(as.numeric(landA), as.numeric(landB), rfscore))) %>%
      ungroup() %>%
      group_by(bullet1, bullet2, landA) %>%
      filter(as.numeric(landB) == reference[[max_diagonal_score[1]]][as.numeric(landA[1]),2]) %>%
      ungroup()
    
    
    # reduced duplicated ones, but it doesn't affect the distribution
    output <- features_sum %>% 
      filter(as.numeric(factor(bullet1)) > as.numeric(factor(bullet2))) %>%
      select(rfscore, bullet1, bullet2, land1, land2)
    
    output
  })
  
  output
}


# read in knm

# this works within a barrel to extract the knm land-land RF comparisons
land_score_knm_generator3_RF <- function(dir, FAUno = NULL, firstno = NULL, remove = NULL) {
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
    features_fau$rfscore <- predict(bulletxtrctr::rtrees, newdata = features_fau, 
                                    type = "prob")[, 2]
    features_slc <- features_fau %>% select(bullet1, bullet2, landA, landB, rfscore, land1, land2)
    
    reference <- list(n1 = data.frame(landA = 1:6, landB = 1:6),
                      n2 = data.frame(landA = 1:6, landB = rep(1:6, times = 2)[2:7]),
                      n3 = data.frame(landA = 1:6, landB = rep(1:6, times = 2)[3:8]),
                      n4 = data.frame(landA = 1:6, landB = rep(1:6, times = 2)[4:9]),
                      n5 = data.frame(landA = 1:6, landB = rep(1:6, times = 2)[5:10]),
                      n6 = data.frame(landA = 1:6, landB = rep(1:6, times = 2)[6:11])) 
    
    features_sum <- features_slc %>% 
      group_by(bullet1, bullet2) %>%
      mutate(max_diagonal_score = which.max(compute_average_scores(as.numeric(landA), as.numeric(landB), rfscore))) %>%
      ungroup() %>%
      group_by(bullet1, bullet2, landA) %>%
      filter(as.numeric(landB) != reference[[max_diagonal_score[1]]][as.numeric(landA[1]),2]) %>%
      ungroup()
    
    
    # reduced duplicated ones, but it doesn't affect the distribution
    output <- features_sum %>% 
      filter(as.numeric(factor(bullet1)) >= as.numeric(factor(bullet2))) %>%
      filter(as.numeric(factor(bullet1)) > as.numeric(factor(bullet2))|as.numeric(factor(landA)) > as.numeric(factor(landB))) %>%
      select(rfscore, bullet1, bullet2, land1, land2)
    
    output
  })
  
  output
}

# A function used in producing tables, extract parameters from cv_summary for two-component models. 
# work with list of models
# estimated_model <- extract_params(cv_summary(list(fau330_km_beta_no_A, fau330_knm_beta_no_A, 
#                                                   fau330_km_beta_A,fau330_knm_beta_A, 
#                                                   full_km_model, full_knm_model)))
# 
# rownames(estimated_model) <- c("fau330_km_beta_no_A", "fau330_knm_beta_no_A", "fau330_km_beta_A", 
#                                "fau330_knm_beta_A", "full_km_model", "full_knm_model")

extract_params <- function(summary_something) {
  
  temp1 <- lapply(summary_something, FUN = function(d){
    p1 <- d$prior[[1]]
    c1_alpha <- d$params[[1]]$alpha
    c1_beta <- d$params[[1]]$beta
    c1_mu <- plogis(d$mu_phi_params[[1]]$mean)
    c1_phi <- exp(d$mu_phi_params[[1]]$precision)
    
    c2_alpha <- d$params[[2]]$alpha
    c2_beta <- d$params[[2]]$beta
    c2_mu <- plogis(d$mu_phi_params[[2]]$mean)
    c2_phi <- exp(d$mu_phi_params[[2]]$precision)
    
    output <- c(p1, c1_alpha, c1_beta, c1_mu, c1_phi, c2_alpha, c2_beta, c2_mu, c2_phi)
    output
  })
  
  temp2 <- lapply(1:9, FUN = function(i){
    
    theparamvector <- lapply(temp1, FUN = function(d2){
      output <- d2[i]
    }) %>% unlist()
    
    theparamvector
  }) 
  
  names(temp2) <-  c("p1", "c1_alpha", "c1_beta", "c1_mu", "c1_phi", "c2_alpha", "c2_beta", "c2_mu", "c2_phi")
  
  temp3 <- as.data.frame(temp2)
  temp3
}









