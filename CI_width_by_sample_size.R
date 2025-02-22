#CI for Cohen's d by sample size in independent samples and paired samples designs
#
#Author: Eric Fields
#Version Date: 19 October 2024

library(MBESS)

proj_dir <- "C:/Users/fieldsec/OneDrive - Westminster College/Documents/ECF/Research/Bayes factor simulations"
setwd(proj_dir)

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

d.pairedsample.CI <- function(dav, N, r, conf.level=0.95) {
  #Calculate the confidence interval for Cohen's d for a paired samples design
  #Formulas adapted from: https://effect-size-calculator.herokuapp.com/
  se <- sqrt( (2*(2 - 2*r)) / (N*2) )
  t <- dav / se
  df <- N - 1
  tCI <- conf.limits.nct(t, df, conf.level=conf.level)
  davL <- tCI$Lower.Limit * se
  davU <- tCI$Upper.Limit * se
  return(c(davL, davU))
}


#Sample sizes
ns <- c(20, 30, 40, 50, 75, 100, 200, 300, 400, 500, 750, 1000, 2000, 5000, 10000)

CI_table <- data.frame()
for (n in ns) {
  
  CI <- d.indsample.CI(0.0, n, n)
  CI_table[sprintf("%d", n), "independent"] <- CI[2]
  
  CI <- d.pairedsample.CI(0.0, n, 0.25)
  CI_table[sprintf("%d", n), "paired (r=0.2)"] <- CI[2]
  
  CI <- d.pairedsample.CI(0.0, n, 0.5)
  CI_table[sprintf("%d", n), "paired (r=0.5)"] <- CI[2]
  
  CI <- d.pairedsample.CI(0.0, n, 0.75)
  CI_table[sprintf("%d", n), "paired (r=0.8)"] <- CI[2]
  
}

write.csv(CI_table, file.path("output", "CI_width.csv"))
