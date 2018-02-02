################################################################################
#                                                                              #                          
# IDEM 156 - Spatial Demography Course 2018                                    #                          
# Contrasting estimation techniques                                            #                          
# Sebastian Klüsener, MPIDR                                                    #                          
#                                                                              #                          
################################################################################

# Erase all objects in workspace
rm(list=ls(all=TRUE))

#install.packages("INLA", repos="https://inla.r-inla-download.org/R/stable", 
#                 dep=TRUE)

# Load libraries
library(RColorBrewer)
library(spdep)
library(rgdal)
library(maptools)
library(stats4)
library(INLA)
library(lrgs)

# Set your working directory
# If you need to change the working directory, you can use the following code:
# Set drive
main <- c("N")
# Set path to session folders
path <- ":/IMPRSD/IDEM 156 Spatial"
# Define session folder
path2 <- "/02 Sessions/08 Spatial Multi-level Modeling"
# Set working directory
setwd(paste(main,path,path2,sep=""))


################################################################################
# 1) Import and prepare data                                                   #
################################################################################

# Open shapefile
shape.shp <- readShapePoly('2004_06', IDvar="KREIS_KENN")
shapeogr.shp <- readOGR(".", "2004_06")
proj4string(shape.shp) <- proj4string(shapeogr.shp)
shape.shp <- spTransform(shape.shp, CRS("+proj=longlat +datum=WGS84"))

# Derive information on internal id and centroids. In this case we turn the id-
# values first into a numeric vector, which is easier to use, if we want to
# match data to the shapefile.
head(shape.shp$KREIS_KENN)
id <- as.numeric(paste(shape.shp$KREIS_KENN))/1000

# Define some standard objects
n <- length(id)
sq <- c(1:n)
sq_text <- paste(sq)
coords <- coordinates(shape.shp)

# Upload data file
mydata <- read.table(file="fertdata2005.csv",  sep=",", header=TRUE)

# Match data to shapefile
data.id <- mydata$ID

# Check, whether vector of data.ids and id of shapefile are of the same lenght
length(data.id) == length(id)

# Match by ID-column, which has to be similar in the shapefile and the datafile
o <- match(id, data.id)
mydata1 <- mydata[o,]
row.names(mydata1) <- shape.shp$KREIS_KENN
shape.shp <- spCbind(shape.shp, mydata1)

y <- mydata$AVAGE05
x1 <- mydata$UER05


################################################################################
# 2) Contrasting model estimations                                             #
################################################################################

# 2.1 Linear model
linear <- lm(y~x1)

# Here we check how long it to to calculate the model
# "user" refers to the time it took to execute our calculation instructions
# Linear models are usually very time efficient
system.time(lm(y~x1))

# Summary of the linear model
summary(linear)
         

# 2.2 Maximum likelihood estimation
LL <- function(beta0, beta1, mu, sigma) {
               R <- y - x1 * beta1 - beta0
               R <- suppressWarnings(dnorm(R, mu, sigma, log = TRUE))
               -sum(R)
}
# As beta0 start value I take the beta0 estimate from the linear model
mle_est  <- mle(LL, start = list(beta0 = 
                                   round(as.numeric(linear$coefficients[1]),1), 
                                 beta1 = 0, mu = 0, sigma=1))

# Maximum likelhood estimation is also quite fast
system.time(mle(LL, start = list(beta0 = 
                                   round(as.numeric(linear$coefficients[1]),1), 
                                 beta1 = 0, mu = 0, sigma=1)))

# Summary of the linear model
summary(mle_est)


# 2.3 Gibbs sampling

gibbs_est  <- Gibbs.regression(x1,y,NULL,Nsamples=100,trace='bsmt', fix='xy')

# The time it takes to calculate a regression with a Gibbs sample depends on
# the number of samples drawn in the Monte Carlo Markov Chain procedure
system.time(Gibbs.regression(x1,y,NULL,Nsamples=50,trace='bsmt', fix='xy'))
system.time(Gibbs.regression(x1,y,NULL,Nsamples=100,trace='bsmt', fix='xy'))
system.time(Gibbs.regression(x1,y,NULL,Nsamples=1000,trace='bsmt', fix='xy'))

b0 <- gibbs_est$B[seq(1,199,2)]
b1 <- gibbs_est$B[seq(2,200,2)]

gb_dat <- as.matrix(data.frame(b0,b1,b1))

# Source: Jeff Gill (c) 2002                                                                             
# "Bayesian Methods: A Social and Behavioral Sciences Approach"                                                        
# http://jgill.wustl.edu/BMSBSA/Chapter09/gibbs.s                              
plot.walk.multi <- function(walk.mat,sim.rm,col)  {
  walk.mat <- walk.mat[,-sim.rm]
  points(walk.mat[1,1],walk.mat[1,2])
  for(i in 1:(nrow(walk.mat)-1))  {
    segments(walk.mat[i,1],walk.mat[i,2],walk.mat[(i+1),1],walk.mat[i,2],col=col)
    segments(walk.mat[(i+1),1],walk.mat[i,2],walk.mat[(i+1),1],
             walk.mat[(i+1),2],col=col)
    
  }
}

# All iterations
par(mfrow=c(1,1))
plot(gb_dat,xlim=range(gb_dat[,1]),ylim=range(gb_dat[,3]),type="n",
     xlab="",ylab="")
mtext(outer=F,side=1,cex=1.1,"b0",line=2)
mtext(outer=F,side=2,cex=1.1,"b1",line=2.25)
plot.walk.multi(gb_dat[,],2,,col="black")


# Usually parameter estimates are worse in the initial period of the sampling.
# Thus one tends to exclude this "burn-in period" 
par(mfrow=c(1,2))
plot(gibbs_est$B[seq(1,199,2)],type="l",ylab="beta0",xlab="Sample")
abline(v=10,col="red")
plot(gibbs_est$B[seq(2,200,2)],type="l",ylab="beta1",xlab="Sample")
abline(v=10,col="red")

mean_b0_burn_in <- mean(gibbs_est$B[seq(1,9,2)])
mean_b1_burn_in <- mean(gibbs_est$B[seq(2,10,2)])

mean_b0_wt_burn_in <- mean(gibbs_est$B[seq(11,199,2)])
mean_b1_wt_burn_in <- mean(gibbs_est$B[seq(12,200,2)])

# All iterations
par(mfrow=c(1,1))
dat <- data.frame(b0,b0)
b0e <- c(unlist(split(dat,1:100))[-1])

dat <- data.frame(b1,b1)
b1e <- c(unlist(split(dat,1:100))[-200])

# Without burn-in period we tend to get closer to the linear model estimate
mean_b0_burn_in
mean_b0_wt_burn_in
linear$coef[1]

mean_b1_burn_in
mean_b1_wt_burn_in
linear$coef[2]

# Estimates
plot(gb_dat,xlim=range(gb_dat[,1]),ylim=range(gb_dat[,3]),type="n",
     xlab="",ylab="")
mtext(outer=F,side=1,cex=1.1,"b0",line=2)
mtext(outer=F,side=2,cex=1.1,"b1",line=2.25)
smoothScatter(b0e,b1e,add=T)
plot.walk.multi(gb_dat[,],2,col="grey50")
abline(v=mean_b0_wt_burn_in,col="orangered1",lwd=2)
abline(h=mean_b1_wt_burn_in,col="orangered1",lwd=2)


# 2.4 Laplace estimation (with INLA)
laplace <-  inla(y~f(x1, model="linear"),data=mydata)
system.time(inla(y~f(x1, model="linear"),data=mydata))

# Summary of the INLA model
summary(laplace)


# Comparison of the results
# All estimations techniques seem to model the relationship quite well
par(mfrow=c(1,1))
plot(x1, y,ylab="Average Age at Birth",xlab="Unemployement Rate",col="grey75")
abline(linear, col = "red",lty=2,lwd=2)
abline(mle_est@coef[c(1,2)], col = "blue",lty=5,lwd=2)
abline(mean_b0_wt_burn_in,mean_b1_wt_burn_in, col = "green2",lty=4,lwd=2)
abline(laplace$summary.fixed$mean, col = "brown",lty=3,lwd=2)




