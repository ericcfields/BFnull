#Example of Bayes factor results with small effect
#
#Author: Eric Fields
#Version Date: 22 February 2025

#Copyright (c) 2025, Eric Fields
#This code is free and open source software made available under the terms 
#of the CC BY 4.0 license
#https://creativecommons.org/licenses/by/4.0/


library(MBESS)
library(BayesFactor)

d.indsample.CI <- function(d, n1, n2, conf.level=0.95) {
  #Calculate the confidence interval for Cohen's d for an independent samples design
  #Formulas adapted from: https://effect-size-calculator.herokuapp.com/
  se <- sqrt(1/n1 + 1/n2)
  t <- d / se
  df <- n1 + n2 -2
  tCI <- conf.limits.nct(t, df, conf.level = conf.level)
  dL <- tCI$Lower.Limit * se
  dU <- tCI$Upper.Limit * se
  return(c(dL, dU))
}

create_mod <- function(dens_fun, lo, hi) {
  #Create a proper density function model from a starting function and an upper and lower bound
  
  #Modify density function to set values below lo and above hi to 0
  f <- function(x) dens_fun(x) * as.integer(x>lo)*as.integer(x<hi)
  
  #Re-scale to make a proper density (area under the curve = 1)
  K <- 1/integrate(f, lower=lo, upper=hi)$value
  return(function(x) f(x) * K)
  
}


################### EXAMPLE PARAMETERS ###################


#sample results for an independent samples experiment
M1 <- 0.0 #mean of condition 1
M2 <- 0.2 #mean of condition 2
s <- 1 #standard deviation in both conditions
n1 <- 40 #sample size in condition 1
n2 <- 40 #sample size in condition 2


################### FREQUENTIST T-TEST & DEFAULT BF ###################

#Cohen's d effect size
mean_diff <- M2 - M1
d <- mean_diff / s

#frequentist t-test
df <- n1 + n2 - 2
SE_diff <- sqrt( (s^2)/n1 + (s^2)/n2 )
t <- mean_diff / SE_diff
pval <- pt(abs(t), df, lower.tail=FALSE) * 2

#confidence intervals
CI80 <- d.indsample.CI(d, n1, n2, conf.level=0.80)
CI95 <- d.indsample.CI(d, n1, n2, conf.level=0.95)

#Bayes factor
H1_scale <- sqrt(2) / 2
lnBF <- ttest.tstat(t, n1, n2, rscale=H1_scale)$bf #ttest.stat function reports natural log of the BF
BF10 <- exp(lnBF) #BF in favor of alternative
BF01 <- 1 / BF10 #BF in favor of null

#report results
cat(sprintf("d = %.2f, n1 = %d, n2 = %d", d, n1, n2),
    sprintf("\nt(%d) = %.2f, p = %.3f", df, t, pval),
    sprintf("80 CI = [%.2f, %.2f]", CI80[1], CI80[2]),
    sprintf("95 CI = [%.2f, %.2f]", CI95[1], CI95[2]),
    sprintf("\nDefault alternative: Cauchy(0, 0.707)"),
    sprintf("defatul BF10 = %.2f", BF10),
    sprintf("default BF01 = %.2f", BF01),
    fill=2)


############################ CUSTOM PRIOR BF ############################

#Specify alternative
f <- function(delta) dnorm(delta, 0.2, 0.1) #starting function for alternative model
lo <- -Inf #probability for values of delta below this will be 0
hi <- Inf #probability for values of delta above this will be 0
alt_dens <- create_mod(f, lo, hi) #create the density function for the alternative model

#Likelihood under the null
L0 <- dt(t, df)

#Likelihood under the alternative
alt_likelihood_integrand <- function(delta, d, n1, n2) {
  t <- d / sqrt( 1/n1 + 1/n2 )
  df <- n1 + n2 - 2
  delta_t <- delta / sqrt( 1/n1 + 1/n2 )
  return(dt(t, df, delta_t) * alt_dens(delta))
}
L1 <- integrate(alt_likelihood_integrand, lower=lo, upper=hi, d=d, n1=n1, n2=n2)$value

#Bayes factor
BF10_custom <- L1/L0
BF01_custom <- L0/L1

#report results
cat(sprintf("\nCustom alternative"),
    sprintf("BF10 = %.2f", BF10_custom),
    sprintf("BF01 = %.2f", BF01_custom),
    fill=2)
