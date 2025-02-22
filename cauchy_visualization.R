#Create visualization and probability table for the Cauchy distribution 
#used in default Bayes factor
#
#Author: Eric Fields
#Version Date: 22 February 2025

#Copyright (c) 2025, Eric Fields
#This code is free and open source software made available under the terms 
#of the CC BY 4.0 license
#https://creativecommons.org/licenses/by/4.0/


library(ggplot2)

proj_dir <- "C:/Users/fieldsec/OneDrive - Westminster College/Documents/ECF/Research/Bayes factor simulations"
setwd(proj_dir)

###### PLOT CAUCHY DISTRIBUTION

#Generate data to plot
Cohens_d <- seq(-3, 3, by=.0001)
cauchy_scale = sqrt(2) / 2
prob_density <- dcauchy(Cohens_d, location = 0, scale = cauchy_scale)
df <- data.frame(x = Cohens_d, y = prob_density)

#Plot
ggplot(df, aes(x = x, y = y)) +
  geom_line(linewidth=2, colour = "darkred") +
  xlab("Cohen's d") +
  ylab("probability density") +
  theme(text = element_text(size = 20)) +
  scale_x_continuous(breaks=c(-3, -2, -1, 0, 1, 2, 3))

#Save figure
ggsave(file.path("output", "Cauchy.707.tif"))


####### TABLE OF CAUCHY PROBABILITIES #######

cauchy_table <- data.frame()
for (d in c(0.2, 0.5, cauchy_scale, 0.8, 1.2, 2.0)) {
  cauchy_table[sprintf("%.3f", d), "probability"] <- pcauchy(d, location = 0, scale = cauchy_scale, lower.tail = FALSE) * 2
}
write.csv(cauchy_table, file.path("output", "cauchy_probabilities.csv"))
