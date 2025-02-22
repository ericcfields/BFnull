#Bayes factor and confidence interval by sample size with a given observed effect size
#
#Author: Eric Fields
#Version Date: 22 February 2025

#Copyright (c) 2025, Eric Fields
#This code is free and open source software made available under the terms 
#of the CC BY 4.0 license
#https://creativecommons.org/licenses/by/4.0/


library(BayesFactor)
library(MBESS)

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


#Parameters
d <- 0.2 #effect size
H1_scale <- sqrt(2)/2 #scale factor for H1 Cauchy


sim_results <- data.frame()
for (n in seq(5, 1000)) {
  
  #Calculate confidence interval
  CI95 <- d.indsample.CI(d, n, n, conf.level=0.95)
  
  #Get the t-statistic
  t <- d / sqrt( 1/n + 1/n )
  
  #Bayes factor in favor of the null
  lnBF <- ttest.tstat(t, n, n, rscale=H1_scale)$bf ##The ttest.stat function reports the natural log of the BF
  BF10 <- exp(lnBF) #BF in favor of alternative
  BF01 <- 1 / BF10 #BF in favor of null
  
  #Add results to table
  sim_results[as.character(n), "BF10"] <- BF10
  sim_results[as.character(n), "CI_L"] <- CI95[1]
  sim_results[as.character(n), "CI_U"] <- CI95[2]
  
}
