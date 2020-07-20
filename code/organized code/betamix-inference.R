# Estimate the Fisher information matrix, apply the asymptotic normality of MLE
# We will make use of polygamma function (derivative of ln of gamma function)
# Full calculation is not presented (checked one more time before coding)
# The package betareg also provides standard error estimate: summary(model), vcov(model)
# AIC BIC can be calculated: logLik(model), AIC(gy2)

# For polygamma function
library(pracma)

# I_ab means the term in Fisher information matrix for partial derivative for alpha_1, then beta_1 (negative sign included)
# I_ab2 means the term in Fisher information matrix for partial derivative for alpha_1, then beta_2 (negative sign included)
# I_ba2 means the term in Fisher information matrix for partial derivative for beta_1, then alpha_2 (negative sign included)
# similar names for other therms
# since the two components in the betamix are treated equally, we have only coded for the first component but it is good for both
# since there is symmetry in beta distribution, we have only coded for alpha, but it is good for beta with some changes in use
# since Fisher information matrix is symmetric, we have only coded for I_ab, not I_ba etc.
# we use delta method for mu and phi parameterization inference

I_pp <- function(x, p1, a1, b1, a2, b2){
  f1 <- dbeta(x, shape1 = a1, shape2 = b1)
  f2 <- dbeta(x, shape1 = a2, shape2 = b2)
  sum((f1-f2)^2/(p1*f1 + (1-p1)*f2)^2)
}
I_pp(c(0.5, 0.6), p1 = 0.55, a1 = 2, b1 = 3, a2 = 3, b2 =2)

# This is not directly available, it needs to be derived from polygamma functions
# First derivative of gamma function
d1_gamma <- function(x){
  gamma(x)*psi(k = 0, z = x)
}

d2_gamma <- function(x){
  (psi(k = 1, z = x) + (d1_gamma(x))^2/gamma(x)^2)*gamma(x)
}

d1a_beta <- function(x, a, b){
  (1-x)^(b-1)*(1/gamma(b))^2*(1/gamma(a))^2*
    (d1_gamma(a+b)*x^(a-1)*gamma(a)*gamma(b) + 
       gamma(a+b)*log(x)*x^(a-1)*gamma(a)*gamma(b) -
       gamma(a+b)*x^(a-1)*d1_gamma(a)*gamma(b))
}

d2a_beta <- function(x, a, b){
  (1-x)^(b-1)*(1/gamma(b))^2*(-2*(1/gamma(a))^3*d1_gamma(a)*
                                (d1_gamma(a+b)*x^(a-1)*gamma(a)*gamma(b) + 
                                   gamma(a+b)*log(x)*x^(a-1)*gamma(a)*gamma(b) -
                                   gamma(a+b)*x^(a-1)*d1_gamma(a)*gamma(b)) +
                                (1/gamma(a))^2*
                                (d2_gamma(a+b)*x^(a-1)*gamma(a)*gamma(b) + 
                                   d1_gamma(a+b)*(log(x)*x^(a-1)*gamma(a)*gamma(b) + x^(a-1)*d1_gamma(a)*gamma(b)) +
                                   d1_gamma(a+b)*log(x)*x^(a-1)*gamma(a)*gamma(b) + 
                                   gamma(a+b)*(log(x)^2*x^(a-1)*gamma(a)*gamma(b) + log(x)*x^(a-1)*d1_gamma(a)*gamma(b)) -
                                   d1_gamma(a+b)*x^(a-1)*d1_gamma(a)*gamma(b) -
                                   gamma(a+b)*(log(x)*x^(a-1)*d1_gamma(a)*gamma(b) + x^(a-1)*d2_gamma(a)*gamma(b))))
}

I_aa <- function(x, p1, a1, b1, a2, b2) {
  f1 <- dbeta(x, shape1 = a1, shape2 = b1)
  f2 <- dbeta(x, shape1 = a2, shape2 = b2)
  
  sum(-p1/(p1*f1+(1-p1)*f2)^2*(d2a_beta(x, a1, b1)*(p1*f1+(1-p1)*f2)-d1a_beta(x, a1, b1)^2*p1))
}

d2ab_beta <- function(x, a, b) {
  log(1-x)*(1-x)^(b-1)*(1/gamma(b))^2*(1/gamma(a))^2*(d1_gamma(a+b)*x^(a-1)*gamma(a)*gamma(b) + 
                                                        gamma(a+b)*log(x)*x^(a-1)*gamma(a)*gamma(b) -
                                                        gamma(a+b)*x^(a-1)*d1_gamma(a)*gamma(b)) + 
    (1-x)^(b-1)*(-2*(1/gamma(b))^3*d1_gamma(b)*(1/gamma(a))^2*(d1_gamma(a+b)*x^(a-1)*gamma(a)*gamma(b) + 
                                                                 gamma(a+b)*log(x)*x^(a-1)*gamma(a)*gamma(b) -
                                                                 gamma(a+b)*x^(a-1)*d1_gamma(a)*gamma(b)) + 
                   (1/gamma(b))^2*(1/gamma(a))^2*(d2_gamma(a+b)*x^(a-1)*gamma(a)*gamma(b) + 
                                                    d1_gamma(a+b)*x^(a-1)*gamma(a)*d1_gamma(b) + 
                                                    d1_gamma(a+b)*log(x)*x^(a-1)*gamma(a)*gamma(b) + 
                                                    gamma(a+b)*log(x)*x^(a-1)*gamma(a)*d1_gamma(b) - 
                                                    d1_gamma(a+b)*x^(a-1)*d1_gamma(a)*gamma(b) -
                                                    gamma(a+b)*x^(a-1)*d1_gamma(a)*d1_gamma(b)))
}

# used d1b_beta(x, a1, b1) = d1a_beta(1-x, b1, a1), the former is not actually coded
I_ab <- function(x, p1, a1, b1, a2, b2) {
  f1 <- dbeta(x, shape1 = a1, shape2 = b1)
  f2 <- dbeta(x, shape1 = a2, shape2 = b2)
  
  sum(-p1/(p1*f1+(1-p1)*f2)^2*(d2ab_beta(x, a1, b1)*(p1*f1+(1-p1)*f2) - d1a_beta(x, a1, b1)*p1*d1a_beta(1-x, b1, a1)))
}

I_pa <- function(x, p1, a1, b1, a2, b2) {
  f1 <- dbeta(x, shape1 = a1, shape2 = b1)
  f2 <- dbeta(x, shape1 = a2, shape2 = b2)
  
  sum(-1/(p1*f1+(1-p1)*f2)^2*(d1a_beta(x, a1, b1)*(p1*f1+(1-p1)*f2) - (f1-f2)*p1*d1a_beta(x, a1, b1)))
}

I_pa2 <- function(x, p1, a1, b1, a2, b2) {
  f1 <- dbeta(x, shape1 = a1, shape2 = b1)
  f2 <- dbeta(x, shape1 = a2, shape2 = b2)
  
  sum(-1/(p1*f1+(1-p1)*f2)^2*(-d1a_beta(x, a2, b2)*(p1*f1+(1-p1)*f2)-(f1-f2)*(1-p1)*d1a_beta(x, a2, b2)))
}

I_aa2 <- function(x, p1, a1, b1, a2, b2) {
  f1 <- dbeta(x, shape1 = a1, shape2 = b1)
  f2 <- dbeta(x, shape1 = a2, shape2 = b2)
  
  sum(-1*(-1)*1/(p1*f1+(1-p1)*f2)^2*p1*d1a_beta(x, a1, b1)*(1-p1)*d1a_beta(x, a2, b2))
}


# Generated, the followings are modifications of the above for different symmetric cases=======================
# from I_pa
I_pb <- function(x, p1, a1, b1, a2, b2) {
  f1 <- dbeta(1-x, shape1 = b1, shape2 = a1)
  f2 <- dbeta(x, shape1 = a2, shape2 = b2)
  
  sum(-1/(p1*f1+(1-p1)*f2)^2*(d1a_beta(1-x, b1, a1)*(p1*f1+(1-p1)*f2) - (f1-f2)*p1*d1a_beta(1-x, b1, a1)))
}

# from I_pa2
I_pb2 <- function(x, p1, a1, b1, a2, b2) {
  f1 <- dbeta(x, shape1 = a1, shape2 = b1)
  f2 <- dbeta(1-x, shape1 = b2, shape2 = a2)
  
  sum(-1/(p1*f1+(1-p1)*f2)^2*(-d1a_beta(1-x, b2, a2)*(p1*f1+(1-p1)*f2)-(f1-f2)*(1-p1)*d1a_beta(1-x, b2, a2)))
}

# from I_aa2
I_ab2 <- function(x, p1, a1, b1, a2, b2) {
  f1 <- dbeta(x, shape1 = a1, shape2 = b1)
  f2 <- dbeta(1-x, shape1 = b2, shape2 = a2)
  
  sum(-1*(-1)*1/(p1*f1+(1-p1)*f2)^2*p1*d1a_beta(x, a1, b1)*(1-p1)*d1a_beta(1-x, b2, a2))
}

# from I_aa
I_bb <- function(x, p1, a1, b1, a2, b2) {
  f1 <- dbeta(1-x, shape1 = b1, shape2 = a1)
  f2 <- dbeta(1-x, shape1 = b2, shape2 = a2)
  
  sum(-p1/(p1*f1+(1-p1)*f2)^2*(d2a_beta(1-x, b1, a1)*(p1*f1+(1-p1)*f2)-d1a_beta(1-x, b1, a1)^2*p1))
}

# from I_aa2
I_ba2 <- function(x, p1, a1, b1, a2, b2) {
  f1 <- dbeta(1-x, shape1 = b1, shape2 = a1)
  f2 <- dbeta(x, shape1 = a2, shape2 = b2)
  
  sum(-1*(-1)*1/(p1*f1+(1-p1)*f2)^2*p1*d1a_beta(1-x, b1, a1)*(1-p1)*d1a_beta(x, a2, b2))
}

# from I_aa2
I_bb2 <- function(x, p1, a1, b1, a2, b2) {
  f1 <- dbeta(1-x, shape1 = b1, shape2 = a1)
  f2 <- dbeta(1-x, shape1 = b2, shape2 = a2)
  
  sum(-1*(-1)*1/(p1*f1+(1-p1)*f2)^2*p1*d1a_beta(1-x, b1, a1)*(1-p1)*d1a_beta(1-x, b2, a2))
}

# from I_aa
I_a2a2 <- function(x, p1, a1, b1, a2, b2) {
  p1 <- 1-p1
  
  f1 <- dbeta(x, shape1 = a2, shape2 = b2)
  f2 <- dbeta(x, shape1 = a1, shape2 = b1)
  
  sum(-p1/(p1*f1+(1-p1)*f2)^2*(d2a_beta(x, a2, b2)*(p1*f1+(1-p1)*f2)-d1a_beta(x, a2, b2)^2*p1))
}

# from I_ab
I_a2b2 <- function(x, p1, a1, b1, a2, b2) {
  p1 <- 1 - p1
  
  f1 <- dbeta(x, shape1 = a2, shape2 = b2)
  f2 <- dbeta(x, shape1 = a1, shape2 = b1)
  
  sum(-p1/(p1*f1+(1-p1)*f2)^2*(d2ab_beta(x, a2, b2)*(p1*f1+(1-p1)*f2) - d1a_beta(x, a2, b2)*p1*d1a_beta(1-x, b2, a2)))
}

# from I_aa
I_b2b2 <- function(x, p1, a1, b1, a2, b2) {
  p1 <- 1-p1
  
  f1 <- dbeta(1-x, shape1 = b2, shape2 = a2)
  f2 <- dbeta(x, shape1 = a1, shape2 = b1)
  
  sum(-p1/(p1*f1+(1-p1)*f2)^2*(d2a_beta(1-x, b2, a2)*(p1*f1+(1-p1)*f2)-d1a_beta(1-x, b2, a2)^2*p1))
}
  
# example and expected final form for our betamix================================
# the model is estimated in the model-full-data.r script, beta_params and cv_summary in the helper function
beta_params(rf_knm_full_model)
cv_summary(list(rf_knm_full_model))
p1 <- 0.7882412
a1 <- 1.245
b1 <- 2.439
a2 <- 9.081
b2 <- 7.814
x <- rf_knm_full_data$rfscore

i_11 <- I_pp(x, p1 = p1, a1 = a1, b1 = b1, a2 = a2, b2 = b2)
i_12 <- I_pa(x, p1 = p1, a1 = a1, b1 = b1, a2 = a2, b2 = b2)
i_13 <- I_pb(x, p1 = p1, a1 = a1, b1 = b1, a2 = a2, b2 = b2)
i_14 <- I_pa2(x, p1 = p1, a1 = a1, b1 = b1, a2 = a2, b2 = b2)
i_15 <- I_pb2(x, p1 = p1, a1 = a1, b1 = b1, a2 = a2, b2 = b2)
i_22 <- I_aa(x, p1 = p1, a1 = a1, b1 = b1, a2 = a2, b2 = b2)
i_23 <- I_ab(x, p1 = p1, a1 = a1, b1 = b1, a2 = a2, b2 = b2)
i_24 <- I_aa2(x, p1 = p1, a1 = a1, b1 = b1, a2 = a2, b2 = b2)
i_25 <- I_ab2(x, p1 = p1, a1 = a1, b1 = b1, a2 = a2, b2 = b2)
i_33 <- I_bb(x, p1 = p1, a1 = a1, b1 = b1, a2 = a2, b2 = b2)
i_34 <- I_ba2(x, p1 = p1, a1 = a1, b1 = b1, a2 = a2, b2 = b2)
i_35 <- I_bb2(x, p1 = p1, a1 = a1, b1 = b1, a2 = a2, b2 = b2)
i_44 <- I_a2a2(x, p1 = p1, a1 = a1, b1 = b1, a2 = a2, b2 = b2)
i_45 <- I_a2b2(x, p1 = p1, a1 = a1, b1 = b1, a2 = a2, b2 = b2)
i_55 <- I_b2b2(x, p1 = p1, a1 = a1, b1 = b1, a2 = a2, b2 = b2)

fisher_information_notation <- matrix(c(
  "pp", "pa", "pb", "pa2", "pb2",
  "n", "aa", "ab", "aa2", "ab2",
  "n", "n", "bb", "ba2", "bb2",
  "n", "n", "n", "a2a2", "a2b2",
  "n", "n", "n", "n", "b2b2"),
  byrow = T,
  nrow = 5)

fisher_information <- matrix(c(
  i_11, i_12, i_13, i_14, i_15,
  i_12, i_22, i_23, i_24, i_25,
  i_13, i_23, i_33, i_34, i_35,
  i_14, i_24, i_34, i_44, i_45,
  i_15, i_25, i_35, i_45, i_55),
  byrow = T,
  nrow = 5)

colnames(fisher_information) <- c("p1", "a1", "b1", "a2", "b2")
rownames(fisher_information) <- c("p1", "a1", "b1", "a2", "b2")
fisher_information

solve(fisher_information)

# delta method to get inference for mu and phi
# take phi_1 for example phi_1 = alpha_1 + beta_1
se_phi1 <- sqrt((4.6117 + 5.6788 + 5.6788 + 21.48046)*10^(-5))
se_phi1

summary(rf_knm_full_model)

exp(1.304228)*0.0052302

# almost same!
