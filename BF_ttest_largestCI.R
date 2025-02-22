#What is the largest CI upper bound for an effect where BF01 > 3
#for the independent samples design with equal sample sizes?
#
#Author: Eric Fields
#Version Date: 29 December 2024

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


#Parameters
cauchy_scale <- sqrt(2)/2
effect_sizes <- seq(0, 1.0, 0.01)
sample_sizes <- seq(3, 100, 1)

max_CI <- 0 #initialize
for (d in effect_sizes) {
  for (n in sample_sizes) {
    
    #Calculate confidence interval
    CI95 <- d.indsample.CI(d, n, n, conf.level=0.95)
    
    #Get the t-statistic
    t <- d / sqrt( 1/n + 1/n )
    
    #Bayes factor in favor of the null
    lnBF <- ttest.tstat(t, n, n, rscale=cauchy_scale)$bf ##The ttest.stat function reports the natural log of the BF
    BF10 <- exp(lnBF) #BF in favor of alternative
    BF01 <- 1 / BF10 #BF in favor of null
    
    #Update data if the BF01 > 3 and a new maximum CI upper bound is found
    if ((BF01 > 3) && max(abs(CI95)) > max_CI) {
      max_CI <- max(abs(CI95))
      max_d <- d
      max_n <- n
      max_BF <- BF01
    }
  }

}

#report results
cat("Largest CI upper bound when BF01 > 3",
    sprintf("d = %.2f, n = %d", max_d, max_n),
    sprintf("BF01 = %.2f", max_BF),
    sprintf("CI upper bound = %.2f", max_CI),
    fill=2)
