#check if bootsraping is correbt - do i need resample within resamples
#check p-value calculation
#do we need to graph anything
#how to compare in 2c

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
n.resamples <- 10000
mu0 <- 0

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
resamples.null.closer <- t.stat.storage$closer - mean(t.stat.storage$closer) + mu0
resamples.null.further <- t.stat.storage$further - mean(t.stat.storage$further) + mu0
resamples.null.diff <- t.stat.storage$diff - mean(t.stat.storage$diff) + mu0

#calculate the mean of resamples - should be 0 on average
mean.resample.closer <- mean(resamples.null.closer)
mean.resample.further <- mean(resamples.null.further)
mean.resample.diff <- mean(resamples.null.diff)

################################################################################
# Part b
################################################################################
#calculate t-statistics on the original sample
t.obs.closer <- mean(dat.closer)/(dat.closer.sd/sqrt(n))
t.obs.further <- mean(dat.further)/(dat.further.sd/sqrt(n))
t.obs.diff <- mean(dat.diff)/(dat.diff.sd/sqrt(n))

# p-value: the proportion of observations that are at least as supportive of Ha
p.val.closer <- mean(resamples.null.closer >= t.obs.closer) #right-tailed test
p.val.further <- mean(resamples.null.further <= t.obs.further) #left-tailed test
p.val.diff <- mean(abs(resamples.null.diff) >= abs(t.obs.diff)) #two-tailed test

#compare with the t-test p-values
t.test.closer <- t.test(dat.closer, mu = 0, alternative = "greater") 
t.test.further <- t.test(dat.further, mu = 0, alternative = "less") 
t.test.diff <- t.test(dat.diff, mu = 0, alternative = "two.sided")

#extract their p-values
t.test.closer.p <- t.test.closer$p.value
t.test.further.p <- t.test.further$p.value
t.test.diff.p <- t.test.diff$p.value

#create a comparison table
comparison.table <-tibble(
  "Case" = c("Closer", "Further", "Difference"),
  "Bootstrap" = c(p.val.closer, p.val.further, p.val.diff),
  "T-test" = c(t.test.closer.p, t.test.further.p, t.test.diff.p)
)

################################################################################
# Part c
################################################################################
#get 5th percentile of the shifted resamples
close.5th <- quantile(resamples.null.closer, probs = 0.05)
further.5th <- quantile(resamples.null.further, probs = 0.05)
diff.5th <- quantile(resamples.null.diff, probs = 0.05)

#calculate theoretical t-distribution value
df <- n-1
#approximate t0.05,n-1
t.critical.5th <- qt(0.05, df)

#comparison table
comparison.table.percentilr <-tibble(
  "Case" = c("Closer", "Further", "Difference"),
  "Bootstrap" = c(close.5th, further.5th, diff.5th),
  "Theoretical" = rep(t.critical.5th, 3)
)

################################################################################
# Part d - compute bootstrap percentile confidence intervals using resamples
################################################################################
mean.storage <- tibble(closer = rep(NA, n.resamples),
                         further = rep(NA, n.resamples),
                         diff = rep(NA, n.resamples))

#perform resampling
for (i in 1:n.resamples){
  #resample
  resample.closer <- sample(dat.closer, size = n, replace = T) 
  resample.further <- sample(dat.further, size = n, replace = T) 
  resample.diff <- sample(dat.diff, size = n, replace = T)
  #calculate and store the mean
  mean.storage$closer[i] = mean(resample.closer)
  mean.storage$further[i] = mean(resample.further)
  mean.storage$diff[i] = mean(resample.diff)
}
#confidence interval for closer data
ci.lower.closer <- quantile(mean.storage$closer, probs = 0.025)
ci.upper.closer <- quantile(mean.storage$closer, probs = 0.975)

#confidence interval for further data
ci.lower.further <- quantile(mean.storage$further, probs = 0.025)
ci.upper.further <- quantile(mean.storage$further, probs = 0.975)

#confidence interval for difference data
ci.lower.diff <- quantile(mean.storage$diff, probs = 0.025)
ci.upper.diff <- quantile(mean.storage$diff, probs = 0.975)

#compute t-test confidence intervals
t.test.closer <- t.test(dat.closer, mu = 0, alternative = "two.sided") 
t.test.further <- t.test(dat.further, mu = 0, alternative = "two.sided") 
t.test.diff <- t.test(dat.diff, mu = 0, alternative = "two.sided")

#extract t-test confidence intervals 
ci.closer.t.test <- t.test.closer$conf.int
ci.closer.t.test.lower <- ci.closer.t.test[1] #get t-test CI for closer data
ci.closer.t.test.upper <- ci.closer.t.test[2]

ci.further.t.test <- t.test.further$conf.int
ci.further.t.test.lower <- ci.further.t.test[1] #get t-test CI for further data
ci.further.t.test.upper <- ci.further.t.test[2]

ci.diff.t.test <- t.test.diff$conf.int
ci.diff.t.test.lower <- ci.diff.t.test[1] #get t-test CI for diff data
ci.diff.t.test.upper <- ci.diff.t.test[2]

#create a comparison table
ci.comparison.table <- tibble(
  "Case" = c("Closer", "Further", "Diff"),
  "Bootstrap CI Lower" = c(ci.lower.closer, ci.lower.further, ci.lower.diff),
  "T-test CI Lower" = c(ci.closer.t.test.lower, ci.further.t.test.lower, ci.diff.t.test.lower),
  "Bootstrap CI Upper" = c(ci.upper.closer, ci.upper.further, ci.upper.diff),
  "T-test CI Upper" = c(ci.closer.t.test.upper, ci.further.t.test.upper, ci.diff.t.test.upper)
)
