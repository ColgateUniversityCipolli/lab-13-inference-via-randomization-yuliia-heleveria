################################################################################
# LAB 13 R CODE
# YULIIA HELEVERIA
# MATH 240 - SPRING 2025
################################################################################

################################################################################
# Load Libraries
################################################################################
library(tidyverse)
library(e1071)
library(boot)

################################################################################
# QUESTION 1
################################################################################

################################################################################
# Part a
################################################################################
#read the data file
dat.finches <- read_csv("zebrafinches.csv")
#separate further column
dat.further <- dat.finches$further

mu0 = 0
#compute t-statistics to get the t-value for the further data
t_furth <- t.test(dat.further, mu = mu0, alternative = "less")
t <- t_furth$statistic #extract t
n <- length(dat.further) #get n

#calculate skewness of data
skew <- skewness(dat.further)

#use Gaussian pdf and calculate the error
pdf.finches <- dnorm(t, mean = 0, sd = 1)
potential.error <- (skew/sqrt(n))*((2*t^2 +1)/6)*pdf.finches

################################################################################
# Part b
################################################################################
#generate a vector of t-values and compute pdf for each of them
t.vals <- seq(from = -10, to = 10, length.out = 1000)
pdf.t.vals <- dnorm(t.vals, mean = 0, sd = 1)

#calculate the error across of all t
t.potential.error <- (skew/sqrt(n))*((2*t.vals^2 +1)/6)*pdf.t.vals

#create a tibble to plot the errors for t-statistics
dat.error.plot <- tibble(t.vals, t.potential.error)

#plot the error for the t-statistics 
error.plot <- ggplot(dat.error.plot)+
  geom_line(aes(x= t.vals, y = t.potential.error))+ #plot the line for the errors
  theme_bw()+
  labs(title = "Edgeworth approximation for the error for p-value across t",
       x = "t-statistics",
       y = "Potential error")

################################################################################
# Part c
################################################################################
alpha <- 0.05
#calculate the critical value for the left-tailed test
t.crit <- qnorm(alpha, mean = 0, sd = 1)

#calculate pdf
pdf.sample <- dnorm(t.crit, mean = 0, sd = 1)

#calculate n 
n.sample <- (skew/(6*(0.10*alpha))*(2*t.crit^2+1)*pdf.sample)^2

################################################################################
# QUESTION 2
################################################################################

################################################################################
# Part a
################################################################################
n.resamples <- 1000

#separate closer and difference data
dat.closer <- dat.finches$closer
dat.diff <- dat.finches$diff

#shift the mean of the data to be 0 under null hypothesis
dat.closer.centered <- dat.closer - mean(dat.closer)

#get standard deviation
dat.closer.sd <- sd(dat.closer)

#store the approximation of T-statistics
t.stat.storage <- tibble(closer = rep(NA, n.resamples))

#perform resampling and calculate T-statistics
for (i in 1:n.resamples){
  resample <- sample(dat.closer, size = n, replace = T)
}
