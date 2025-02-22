#Simulations to explore frequentist properties of the Bayesian t-test
#
#Author: Eric Fields
#Version Date: 22 February 2025

#Copyright (c) 2025, Eric Fields
#This code is free and open source software made available under the terms 
#of the CC BY 4.0 license
#https://creativecommons.org/licenses/by/4.0/


library(BayesFactor)

proj_dir <- "C:/Users/fieldsec/OneDrive - Westminster College/Documents/ECF/Research/Bayes factor simulations"
setwd(proj_dir)


################### SIMULATION PARAMETERS ###################

#Number of simulated studies for each set of parameters
n_sim <- 1000

#Simulated data parameters
pop_d <- c(0.0, 0.2, 0.5, 0.8, 1.2) #Population Cohen's d
ns <- c(20, 40, 100, 250, 400) #Sample size in each group

#Test parameters
alpha <- 0.05 #Alpha level for frequentist t-test
H1_scales <- c(0.2, 0.4, sqrt(2)/2) #Scale factor for BF H1 distribution ("default" is sqrt(2)/2)


################### RUN SIMULTAIONS ###################

for (n in ns) {

  #Initialize results data frame
  sim_results <- data.frame()
  
  #Initialize BF and p-value arrays
  bf <- array(NaN, dim=c(n_sim, length(pop_d), length(H1_scales)))
  p_vals <- array(NaN, dim=c(n_sim, length(pop_d)))
  
  for (j in 1:length(pop_d)) {
    
    d <- pop_d[j]
    
    #Simulate data and conduct tests
    for (i in 1:n_sim) {
      
      #Generate data
      X <- rnorm(n, mean=0, sd=1)
      Y <- rnorm(n, mean=d, sd=1)
      
      #frequentist t-test
      t_results <- t.test(X, Y)
      p_vals[i, j] <- t_results$p.value
      
      #Bayesian t-test
      for (k in 1:length(H1_scales)) {
        bayes_results = ttestBF(x=X, y=Y, rscale=H1_scales[k])
        bf[i, j, k] <- as.vector(bayes_results)
      }
    }
    
    #Add results to table
    sim_results["t-test rej", sprintf("d = %.2f", d)] <- mean(p_vals[, j] <= alpha)
    for (k in 1:length(H1_scales)) {
      alt_text <- sprintf("H1 scale = %.3f", H1_scales[k])
      sim_results[sprintf("%s BF01>3", alt_text), sprintf("d = %.2f", d)] <- mean(bf[, j, k] < (1/3))
      sim_results[sprintf("%s BF01>10", alt_text), sprintf("d = %.2f", d)] <- mean(bf[, j, k] < (1/10))
      sim_results[sprintf("%s BF10>3", alt_text), sprintf("d = %.2f", d)] <- mean(bf[, j, k] > 3)
      sim_results[sprintf("%s BF10>10", alt_text), sprintf("d = %.2f", d)] <- mean(bf[, j, k] > 10)
    }
    
  }
  
  #Save simulation results to file
  write.csv(sim_results, file.path(proj_dir, "output", sprintf("BF_sim_results_n%d.csv", n)))

}

