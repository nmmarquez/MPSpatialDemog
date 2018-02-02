################################################################################
#                                                                              #
# IDEM 156 - Spatial Demography Course 2018                                    #                             
# How non-independence affects the variance and the modelling results          #                          
# Or: The whole secret why we need spatial analysis techniques                 #                          
# Sebastian Klüsener, MPIDR                                                    #                          
#                                                                              #                          
################################################################################


# Erase all objects in workspace
rm(list=ls(all=TRUE))

# Create random sample
random <- rnorm(100, mean=40, sd=5)
random
var(random)

# Insert all observations 100 times
non.random <- c(rep(random,100))
# Variance is very similar compared to the one of the random data...
var(non.random)
# ...although the data consists of subgroups of data with no variance at all
var(sort(non.random)[1:100])
var(sort(non.random)[101:200])

# Create second random sample
random2 <- rnorm(100, mean=40, sd=5)

# Insert all observations 100 times
non.random2 <- c(rep(random2,100))

# Do a OLS-Model with one random sample as dependent variable and another 
# as covariate - Dependending on your random numbers, results will be in most 
# of the cases not be significant.
model.random <- lm(random ~ random2)

# Calculate the same model with the non-independent-sample
# This model has a much higher likeliness to create significant results.
model.non.random <- lm(non.random ~ non.random2)

# Observe, how the degrees of freedom (DF), standard errors, t-statistics
# F-statistics and r2-statistcs change if you calculate the model with 
# non-independent observations
summary(model.random)
summary(model.non.random)