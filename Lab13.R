#check if bootsraping is correbt

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
library(boot.pval)

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

#get standard deviation
dat.closer.sd <- sd(dat.closer)
dat.further.sd <- sd(dat.further)
dat.diff.sd <- sd(dat.diff)

#store the approximation of T-statistics
t.stat.storage <- tibble(closer = rep(NA, n.resamples),
                         further = rep(NA, n.resamples),
                         diff = rep(NA, n.resamples))

#perform resampling and calculate T-statistics
for (i in 1:n.resamples){
  #resample
  resample.closer <- sample(dat.closer, size = n, replace = T) 
  resample.further <- sample(dat.further, size = n, replace = T) 
  resample.diff <- sample(dat.diff, size = n, replace = T)
  #calculate T
  t.stat.closer <- mean(resample.closer)/(dat.closer.sd/sqrt(n)) 
  t.stat.further <- mean(resample.further)/(dat.further.sd/sqrt(n))
  t.stat.diff <- mean(resample.diff)/(dat.diff.sd/sqrt(n))
  #store the statistics
  t.stat.storage$closer[i] = t.stat.closer 
  t.stat.storage$further[i] = t.stat.further 
  t.stat.storage$diff[i] = t.stat.diff 
}

#store the shifted resamples in an object and shift the mean of the data to be 0 under null hypothesis
resamples.null.closer <- t.stat.storage$closer - mean(t.stat.storage$closer)
resamples.null.further <- t.stat.storage$further - mean(t.stat.storage$further)
resamples.null.diff <- t.stat.storage$diff - mean(t.stat.storage$diff)

#calculate the mean of resamples - should be 0 on average
mean.resample.closer <- mean(resamples.null.closer)
mean.resample.further <- mean(resamples.null.further)
mean.resample.diff <- mean(resamples.null.diff)

################################################################################
# Part b
################################################################################
nresamps <- 10000
resamples <- tibble(tstat.closer=rep(NA, nresamps),
                    tstat.further=rep(NA, nresamps),
                    tstat.diff=rep(NA, nresamps))
for(i in 1:nresamps){
  curr.resample.closer <- sample(x = dat.closer,
                          size = n,
                          replace = T)
  curr.resample.further <- sample(x = dat.further,
                                 size = n,
                                 replace = T)
  curr.resample.diff <- sample(x = dat.diff,
                                 size = n,
                                 replace = T)
  resamples$tstat.closer[i] <- mean(curr.resample.closer)/(dat.closer.sd/sqrt(n))
  resamples$tstat.further[i] <- mean(curr.resample.further)/(dat.further.sd/sqrt(n))
  resamples$tstat.diff[i] <- mean(curr.resample.diff)/(dat.diff.sd/sqrt(n))
}

# Shift the resamples to match the null
shifted.resamples.closer <- resamples$tstat.closer - abs(mean(resamples$tstat.closer))
shifted.resamples.further <- resamples$tstat.further + abs(mean(resamples$tstat.further))
shifted.resamples.diff.low <- resamples$tstat.diff - abs(mean(resamples$tstat.diff))
shifted.resamples.diff.high <- resamples$tstat.diff + abs(mean(resamples$tstat.diff))

# p-value: the proportion of observations that are at least as supportive of Ha
p.boot.closer <- mean(shifted.resamples.closer >= mean(resamples$tstat.closer))
mean(shifted.resamples.further <= mean(resamples$tstat.further))
mean(shifted.resamples.diff.low <= mean(resamples$tstat.diff))
mean(shifted.resamples.diff.high >= mean(resamples$tstat.diff))

# Helper Function
boot.t <- function(d, ind){
  mean(d[ind])/(sd)
}

#compute bootstrap p-value for shifted resamples
boots.closer <- boot(data=dat.closer, 
              statistic = boot_t_test(), 
              R = 30000)

#do bootsraping t_test
boot_t_test(dat.closer, mu= 0, R = 30000, alternative = "greater")

#calculate the p-value
t.closer <- t.test(resamples.null.closer, alternative = "greater", mu = 0)
p.closer.resample <- t.closer$p.value
