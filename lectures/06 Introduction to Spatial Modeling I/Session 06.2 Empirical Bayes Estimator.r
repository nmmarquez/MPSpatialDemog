################################################################################
#                                                                              #
# IDEM 156 - Spatial Demography Course 2018                                    #
# Calculating and Mapping Global and Local Empirical Bayes Estimator           #
# Sebastian Kl?sener, MPIDR                                                    #
#                                                                              #
################################################################################

# Erase all objects in workspace
rm(list=ls(all=TRUE))

library(spdep)
library(maptools)
library(RColorBrewer)


################################################################################
# 1) Import data                                                               #
################################################################################

# Data: For the city districts of Auckland (New Zealand) we have
#       - the number of child deaths in the period 1977-1985 (M77_85)
#       - the number of children under 5 in the year 1981 (Und_81)
example(auckland)

# Calculating raw rates
raw.rates <- auckland$M77_85/(9*auckland$Und5_81)*1000
mean.rate <- sum(auckland$M77_85)/sum((9*auckland$Und5_81))*1000
plot(auckland$Und5_81, raw.rates, ylab="Mortality Rates", 
 xlab="Population at Risk", main="Raw Rates and Population at Risk")


################################################################################
# 2) Global Emprirical Bayes Rates                                             #
################################################################################

# Calculating Global Empirical Bayes Rates
res.glo <- EBest(auckland$M77_85, 9*auckland$Und5_81)

# Mean of raw rates
m <- mean.rate

# This plot illustrates, that Emprical Bayes Rates are shrunken to the mean 
# of the raw rates (the two lines in the graph)

par(mfrow=c(1,1))
plot(res.glo$raw*1000,res.glo$raw*1000, ylab="Empirical Bayes Rate", 
 xlab="Raw Rate", col="white", cex=1)
points(res.glo$raw*1000,res.glo$raw*1000, col="black",cex=1)
points(res.glo$raw*1000,res.glo$estmm*1000, col="red" ,cex=1)
abline(v=m)
abline(h=m)
legend("topleft", c("Raw Rate/ Raw Rate", "Raw Rate/ Global Empirical Bayes"), 
pch=1,col=c("black", "red"), bg="white")

# This plot illustrates, that population at risk has an impact on how the rates 
# are shrunken to the mean rate observed

par(mfrow=c(1,1))
plot(res.glo$raw*1000,res.glo$raw*1000, ylab="Empirical Bayes Rate", 
 xlab="Raw Rate", col="white", cex=1)
text(res.glo$raw*1000,res.glo$raw*1000,label=(auckland$Und5_81*9), 
 col="black",cex=0.8)
text(res.glo$raw*1000,res.glo$estmm*1000,label=(auckland$Und5_81*9), 
 col="red" ,cex=0.8)
abline(v=m)
abline(h=m)
legend("topleft", c("Raw Rate/ Raw Rate - Pop. at Risk", 
 "Raw Rate/ Global Empirical Bayes Rate - Pop. at Risk"), pch=1,
  col=c("black", "red"), bg="white")

################################################################################
# 3) Local Emprirical Bayes Rates                                              #
################################################################################

# Calculating Local Empirical Bayes Rates based on a Neighborhood Weight Matrix
res.loc <- EBlocal(auckland$M77_85, 9*auckland$Und5_81, auckland.nb)

par(mfrow=c(1,1))
plot(res.loc$raw*1000,res.loc$raw*1000, ylab="Local Empirical Bayes Rate", 
 xlab="Raw Rate", col="white", cex=1)
points(res.loc$raw*1000,res.loc$raw*1000, col="black",cex=1)
points(res.loc$raw*1000,res.glo$estmm*1000, col="red" ,cex=1)
points(res.loc$raw*1000,res.loc$est*1000, col="green3" ,cex=1)
abline(v=m)
abline(h=m)
legend("topleft", c("Raw Rate/ Raw Rate", 
 "Raw Rate/ Global Empirical Bayes Rate", 
 "Raw Rate/ Local Empirical Bayes Rate"), 
  pch=1,col=c("black", "red", "green3"), bg="white")

################################################################################
# 4) Comparing Maps                                                            #
################################################################################

# We use for all maps the same categorization bins and colors
brks <- c(-Inf,2,2.5,3,3.5,Inf)
cols <- brewer.pal(5,"Reds")

# Contrasting Raw Rates and Global Empirical Bayes Estimator
par(mfrow=c(1,2))                          
plot(auckland, col=cols[findInterval(raw.rates, brks, all.inside=TRUE)],
     border="grey25")
legend("bottomleft", fill=cols, legend=leglabs(brks), bty="n",cex=0.8)
title(main="Raw child mortality rates")

plot(auckland, col=cols[findInterval(res.glo$estmm*1000, brks, 
     all.inside=TRUE)],border="grey25")
legend("bottomleft", fill=cols, legend=leglabs(brks), bty="n",cex=0.8)
title(main="Global Empirical Bayes -\n child mortality rate")

# Contrast Global and Local Empirical Bayes Estimator
par(mfrow=c(1,2))
plot(auckland, col=cols[findInterval(res.glo$estmm*1000, brks, 
     all.inside=TRUE)],border="grey25")
legend("bottomleft", fill=cols, legend=leglabs(brks), bty="n",cex=0.8)
title(main="Global Empirical Bayes -\n child mortality rate")

plot(auckland, col=cols[findInterval(res.loc$est*1000, brks, 
                                     all.inside=TRUE)],border="grey25")
legend("bottomleft", fill=cols, legend=leglabs(brks), bty="n",cex=0.8)
title(main="Local Moment Estimator -\n child mortality rate")

# Contrast Raw Rate and Local Empirical Bayes Estimator
par(mfrow=c(1,2))                          
plot(auckland, col=cols[findInterval(raw.rates, brks, all.inside=TRUE)],
     border="grey25")
legend("bottomleft", fill=cols, legend=leglabs(brks), bty="n",cex=0.8)
title(main="Raw child mortality rates")

plot(auckland, col=cols[findInterval(res.loc$est*1000, brks, 
                                     all.inside=TRUE)],border="grey25")
legend("bottomleft", fill=cols, legend=leglabs(brks), bty="n")
title(main="Local Moment Estimator -\n child mortality rate",cex=0.8)


################################################################################
# 5) Local Moran's I on raw rates and rates smoothed by Empirical Bayes        #
#    Algorithm                                                                 #
################################################################################

# 5.1) Raw Rates

# Definitions
significance <- 0.05
plot.only.significant <- T
shape.shp <- auckland
n <- length(raw.rates)
varofint <- raw.rates
varofint.name <- "Raw Rates"
varlag <- lag.listw(nb2listw(auckland.nb), varofint)

# Calculate Lisa Test
lisa.raw <- localmoran(varofint,nb2listw(auckland.nb), alternative="two.sided")

# Get significance level
vec <- c(1:n)
vec <- ifelse(lisa.raw[,5] < significance, 1,0)

# Calculate Mean of Variable of interest and lagged value of variable of 
# interest
m.varofint <- mean(varofint)
m.varlag <- mean(varlag)

# Derive sector
sec <- c(1:n)
for (i in 1:n) {
    if (varofint[[i]]>=m.varofint & varlag[[i]]>=m.varlag) sec[i] <- 1
    if (varofint[[i]]<m.varofint & varlag[[i]]<m.varlag) sec[i] <- 2
    if (varofint[[i]]<m.varofint & varlag[[i]]>=m.varlag) sec[i] <- 3
    if (varofint[[i]]>=m.varofint & varlag[[i]]<m.varlag) sec[i] <- 4
}

# Define colors for sectors
sec.all <- sec
colors1 <- c(1:n)
for (i in 1:n) {
    if (sec.all[i]==1) colors1[i] <- "brown2"
    if (sec.all[i]==2) colors1[i] <- "royalblue3"
    if (sec.all[i]==3) colors1[i] <- "lightblue"
    if (sec.all[i]==4) colors1[i] <- "pink"
    if (sec.all[i]==0) colors1[i] <- "white"
}

# Mark all non-significant regions white
loc.m.data <- sec*vec
colors2 <- colors1
for (i in 1:n) {
    if (loc.m.data[i]==0) colors2[i] <- "white"
}

# Cluster map
par(mfrow=c(1,2))  
if (plot.only.significant==TRUE) 
   {plot(shape.shp, col=colors2, border="grey25")} else 
   {plot(shape.shp, col=colors1, border="grey25")}
legend("bottomleft",fill=c("brown2","royalblue3","lightblue","pink","white"),
       legend=c("High-High","Low-Low","Low-High","High-Low"),cex=0.8,bg="white")
title(paste("Significant Clusters\n",varofint.name))


# 5.2) Local Rates

# Definitions
significance <- 0.05
plot.only.significant <- T
shape.shp <- auckland
n <- length(res.loc$est*1000)
varofint <- res.loc$est*1000
varofint.name <- "Local Rates"
varlag <- lag.listw(nb2listw(auckland.nb), varofint)

# Calculate Lisa Test
lisa.raw <- localmoran(varofint,nb2listw(auckland.nb), alternative="two.sided")

# Get significance level
vec <- c(1:n)
vec <- ifelse(lisa.raw[,5] < significance, 1,0)

# Calculate Mean of Variable of interest and lagged value of variable of 
# interest
m.varofint <- mean(varofint)
m.varlag <- mean(varlag)

# Derive sector
sec <- c(1:n)
for (i in 1:n) {
  if (varofint[[i]]>=m.varofint & varlag[[i]]>=m.varlag) sec[i] <- 1
  if (varofint[[i]]<m.varofint & varlag[[i]]<m.varlag) sec[i] <- 2
  if (varofint[[i]]<m.varofint & varlag[[i]]>=m.varlag) sec[i] <- 3
  if (varofint[[i]]>=m.varofint & varlag[[i]]<m.varlag) sec[i] <- 4
}

# Define colors for sectors
sec.all <- sec
colors1 <- c(1:n)
for (i in 1:n) {
  if (sec.all[i]==1) colors1[i] <- "brown2"
  if (sec.all[i]==2) colors1[i] <- "royalblue3"
  if (sec.all[i]==3) colors1[i] <- "lightblue"
  if (sec.all[i]==4) colors1[i] <- "pink"
  if (sec.all[i]==0) colors1[i] <- "white"
}

# Mark all non-significant regions white
loc.m.data <- sec*vec
colors2 <- colors1
for (i in 1:n) {
    if (loc.m.data[i]==0) colors2[i] <- "white"
}

# Cluster map 
if (plot.only.significant==TRUE) 
  {plot(shape.shp, col=colors2, border="grey25")} else 
  {plot(shape.shp, col=colors1, border="grey25")}
legend("bottomleft",fill=c("brown2","royalblue2","lightblue","pink","white"),
       legend=c("High-High","Low-Low","Low-High","High-Low"),cex=0.8,bg="white")
title(paste("Significant Clusters\n",varofint.name))
