################################################################################
#                                                                              #                         
# IDEM 156 - Spatial Demography Course 2018                                    #                          
# Calculating Global Indicators of Spatial Autocorrelation                     #                                                                                                                                 #
# Sebastian Kl?sener, MPIDR                                                    #                         
#                                                                              #                         
################################################################################

# Erase all objects in workspace
rm(list=ls(all=TRUE))

#install.packages(c("evd")),

# Load libraries
library(spdep)
library(maptools)
library(RColorBrewer)
library(evd)
library(rgdal)

# Set your working directory
# If you need to change the working directory, you can use the following code:
# Set drive


################################################################################
# 1) Import and prepare data                                                   #
################################################################################ 

shape.shp <- readShapePoly('2004_06', IDvar="KREIS_KENN")
shapeogr.shp <- readOGR(".", "2004_06")
proj4string(shape.shp) <- proj4string(shapeogr.shp)
shape.ll <- spTransform(shape.shp, CRS("+proj=longlat +datum=WGS84"))

plot(shape.shp)                   
names(shape.shp)

# Derive information on internal id and centroids. In this case we turn the id-
# values first into a numeric vector, which is easier to use, if we want to 
# match data to the shapefile.
head(shape.shp$KREIS_KENN)
is.vector(shape.shp$KREIS_KENN)
id <- as.numeric(paste(shape.shp$KREIS_KENN))/1000
head(id)
is.vector(id)

n <- length(id)
sq <- c(1:n)
sq_text <- paste(sq)

# Coordinates of UTM-projected shapefile
coords <- coordinates(shape.shp)
# Coordinates from reprojected Shapefile with LotLang-projection 
coords.ll <- coordinates(shape.ll)

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
names(shape.shp)

# In order to "attach" data from a shapefile attribute table, we have to call it
# with "@data"
attach(shape.shp@data)
               
# Create Neighborhood Weight Matrix 
nb.FOQ <- poly2nb(shape.shp, queen=TRUE)
summary(nb.FOQ)

# Alternative way to connect R|gen to Stralsund and Vorpommern
nb.FOQ[[350]]
nb.FOQ[[358]]
nb.FOQ[[362]]

nb.FOQ.cor <- nb.FOQ
nb.FOQ.cor[[362]]  <- as.integer(c(350, 358))
nb.FOQ.cor[[350]]  <- as.integer(c(358, 362))
nb.FOQ.cor[[358]]  <- as.integer(c(350,352,353,354,360,362))

# Create a nearest neighbor weigh matrix
l.5NN <- knearneigh(coords.ll, k=5, longlat=T)
nb.5NN <- knn2nb(l.5NN, row.names=id)
#plot.nb(nb.5NN, coords)

# Create a distance based matrix (distance between centroids)
# Distance of 20 km
d20km <- dnearneigh(coords.ll, 0, 20, row.names = id, longlat=T)
# Distance of 50 km
d50km <- dnearneigh(coords.ll, 0, 50, row.names = id,longlat=T)


################################################################################
# 2) Correlation and Spatial Auto-Correlation                                  #
################################################################################ 

# Normal correlation plot
y <- AVAGE05
x <- UER05
y.name <- c("Average Age at Birth 2005")
x.name <- c("Unemployment Rate 2005")

# Create color vector to highlight West and East German districts in different
# colors
color <- c(rep(0,n))
for (i in 1:n) {
color[i] <- if (DUM_E[i]==1) paste("blue") else paste("orange")
}
plot(x,y, col=color, xlab=x.name, ylab=y.name)
abline(v=mean(x), col="darkgrey", lty=2)                
abline(h=mean(y), col="darkgrey", lty=2)  
legend("topright", paste("Correlation",round(cor(x,y),2)))

# Correlation test
cortest <- cor.test(x,y)
cortest

# Plot of Average Age and Spatially lagged variable "Average Age"
y.name <- c("Spatially Lagged Average Age at Birth 2005")
x.name <- c("Average Age at Birth 2005")

# Calculate lagged value of way based on neighborhood definition
lag.y <- lag.listw(nb2listw(nb.FOQ.cor),y)

par(mfrow=c(1,1))
plot(y,lag.y, col=color, xlab=x.name, ylab=y.name)
abline(v=mean(y), col="darkgrey", lty=2)                
abline(h=mean(lag.y), col="darkgrey", lty=2)  
# We use the row-standardized weight matrix, for which the Moran'I Index 
# is the slope of the regression with lag.y as dependent variable and 
# x as covariate
abline(b <- lm(lag.y~y))

par(mfrow=c(1,1))
# For such a plot the moran.plot function exists in the package spdep. The 
# labels identify observations with high influence.
moran.plot(y,nb2listw(nb.FOQ.cor), col=color, 
           labels=GEN,xlab=x.name,ylab=y.name)
 

################################################################################
# 3) Create draft maps of some variables                                       #
################################################################################

# Check, which variables we have available
names(shape.shp)

# Shortly produce draft overview maps
var.data <- data.frame(SNM05, TFR05, UER05, DUM_E)
var.names <-c("Share Nonmarital Births 2005", "TFR 2005", 
"Unemployment Rate", "Dummy East")
m <- length(var.data)-1
colors <- matrix(nrow=n,ncol=m, data=rep(0,n*m))

par(mfrow=c(2,2))
# This code produces three maps with equal interval categorization
for (i in 1:m) {
varofint <- var.data[,i]

lb       <- 10
bins     <- seq(min(varofint), max(varofint), length=lb)
colpal   <- brewer.pal(lb-1, "YlOrRd")
colors[,i]   <- colpal[findInterval(varofint, bins, rightmost.closed=T)]
plot(shape.shp, col=colors[,i])
title(var.names[[i]], cex=1)
legend("bottomright",fill=colpal,legend=paste(round(bins[-length(bins)],2)
,"-",round(bins[-1],2)),cex=0.8, bg="white")
}
# This code produces one maps with a binary categorization (Dummy East)
colore <- c(rep(0,n))
for (i in 1:n) {
colore[i] <- if (DUM_E[i]==1) paste(colpal[[6]]) else paste(colpal[[2]]) 
}
lb       <- 10
colpal   <- brewer.pal(lb-1, "YlOrRd")
varofint <- var.data[,4]
plot(shape.shp, col=colore)
title(var.names[[4]], cex=1)
legend("bottomright",fill=c(colpal[[6]],colpal[[2]]), 
legend=paste(c("West","East")),cex=0.8, bg="white")


################################################################################
# 4) Calculate Moran's I                                                       #
################################################################################

#Create listwise objects
nb.FOQ.lw.W <- nb2listw(nb.FOQ.cor)
nb.FOQ.lw.B <- nb2listw(nb.FOQ.cor, style="B")
nb.d50km.lw.W <- nb2listw(d50km)
nb.5NN.lw.W <- nb2listw(nb.5NN)

# Create a vector of random number of the lenght n (in our case 439 district)
random <- rnorm(n, mean=0, sd=20)

################################################################################

# Contrast random variable with share nonmarital births
# Either activate (by removing the "#") the first two lines of code 
# to look at the random variable, or the third and fourth line to look at share 
# of nonmarital births.

# Random variable might be spatially autocorrelated by chance
var.of.int <- random
var.name <- c("Random")
# The share of nonmarital births is strongly spatially autocorrelated.
#var.of.int <- SNM05
#var.name <- c("SNM05")

################################################################################

# Calculate Moran's I based on first order queen adjacency
# Assumption of Normalization
mt.avag.nb.FOQ <- moran.test(var.of.int, nb.FOQ.lw.W, randomisation=F,
                             alternative="greater")
mt.avag.nb.FOQ

# Assumption of Randomization (adjusted for kurtosis)
mt.avag.nb.FOQ <- moran.test(var.of.int, nb.FOQ.lw.W, alternative="greater")
mt.avag.nb.FOQ

# Permutation test for Moran's I statistic
# Define number of permutations
per <- 9999
mt.avag.mc.nb.FOQ <- moran.mc(var.of.int, nb.FOQ.lw.W, per)
mt.avag.mc.nb.FOQ


# Plots of observed distribution of y and lagged y
x <- var.of.int
y <- lag.listw(nb2listw(nb.FOQ.cor),var.of.int)
xhist <- hist(x, plot=FALSE, breaks=n/20)
yhist <- hist(y, plot=FALSE, breaks=n/20)
top <- max(c(xhist$counts, yhist$counts))
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)
par(mar=c(3,3,1,1))
plot(x, y, xlab="", ylab="",pch=21, col="darkorange1", bg="darkorange1")
par(mar=c(0,3,1,1))
barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0,col="steelblue1",
 main=paste(var.name, "(x-axis) and spatially lagged",var.name,"(y-axis)"), 
 border="steelblue2")
par(mar=c(3,0,1,1))
barplot(yhist$counts, axes=FALSE, xlim=c(0, top),
 space=0, horiz=TRUE,col="steelblue1", border="steelblue2")


# Plot of the permutation test                                                        
par(mfrow=c(1,1), mar=rep(4,4))
hist(mt.avag.mc.nb.FOQ$res[c(1:per+1)],xlim=c(-1,1), col="steelblue1",
     border="steelblue2",breaks=seq(-1,1, by=0.01),
     xlab=paste("Moran I's obtained in ", per," permutations and observed (red)"), 
     main=paste(var.name, 
     "- Observed Moran's I (red)\n in contrast to values obtained by permutation"))
     abline(v=mt.avag.mc.nb.FOQ$res[c(per+1)], col="darkorange1",lwd=2)

# Contrast Moran's I based on different weigh matrices
mt.avag.nb.FOQ  <- moran.test(SNM05, nb.FOQ.lw.W)
mt.avag.nb.d50km <- moran.test(SNM05, nb.d50km.lw.W)
mt.avag.nb.5NN <- moran.test(SNM05, nb.5NN.lw.W)

mt.avag.nb.FOQ
mt.avag.nb.d50km
mt.avag.nb.5NN
                                            
# Contrast Moran's I scatterplots for different weight matrices
# These scatterplots are helpful to identify outliers
par(mfrow=c(2,2))
moran.plot(SNM05, nb.FOQ.lw.W, main="FOQ", 
 xlab="Share Nonmarital Births",
 ylab="Spatially Lagged Share Nonmarital Births", label=GEN)
moran.plot(SNM05, nb.d50km.lw.W, main="50km",
 xlab="Share Nonmarital Births",
 ylab="Spatially Lagged Share Nonmarital Births", label=GEN)
moran.plot(SNM05, nb.5NN.lw.W, main="5NN",
 xlab="Share Nonmarital Births",
 ylab="Spatially Lagged Share Nonmarital Births", label=GEN)

# Calculate Moran's I based on first order queen adjacency for all 
# four variables
mt.avag.nb.FOQ <- moran.test(AVAGE05, nb.FOQ.lw.W)
mt.snm.nb.FOQ <- moran.test(SNM05, nb.FOQ.lw.W)
mt.tfr.nb.FOQ <- moran.test(TFR05, nb.FOQ.lw.W)
mt.east.nb.FOQ <- moran.test(DUM_E, nb.FOQ.lw.W)

mt.avag.nb.FOQ
mt.snm.nb.FOQ
mt.tfr.nb.FOQ
mt.east.nb.FOQ

# Contrast Moran's I scatterplots 
par(mfrow=c(2,2))
moran.plot(SNM05, nb.FOQ.lw.W, main="Average Age at Birth 2005 (FOQ)",
xlab="Average Age",ylab="Spatially Lagged Share Nonmarital Births", label=GEN)
moran.plot(TFR05, nb.FOQ.lw.W, main="Share Nonmarital Births 2005",
xlab="Share Nonmarital Births",ylab="Spatially Lagged Share Nonmarital Births",
label=GEN)
moran.plot(UER05, nb.FOQ.lw.W, main="TFR 2005",
xlab="TFR",ylab="Spatially Lagged TFR", label=GEN)
moran.plot(DUM_E, nb.FOQ.lw.W, main="East",
xlab="East",ylab="Spatially Lagged East", label=GEN)


################################################################################
# 5) Application of Moran's I test in modelling procedures                     #
#    Simple Models testing for effects of unemployment on                      #
#    share non-marital births and total fertility rate                         #
#    (Attempt to control for spatial heterogeneity                             #
################################################################################

# Define variables amd models
y1 <- SNM05
y2 <- TFR05
x1 <- UER05
model1 <- y1 ~ x1
model2 <- y2 ~ x1

# Calculate Moran's I test for depdendent variables
mt.y1.nb.FOQ <- moran.test(y1, nb.FOQ.lw.W)
mt.y2.nb.FOQ <- moran.test(y2, nb.FOQ.lw.W)
mt.y1.nb.FOQ
mt.y2.nb.FOQ

# Calculate models
model1 <- lm(model1)
model2 <- lm(model2)
summary(model1)
summary(model2)

# Run lm.morantest to test, whether residuals are autocorrelated
mt.lmres.y1.nb.FOQ <- lm.morantest(model1, nb.FOQ.lw.W)
mt.lmres.y2.nb.FOQ <- lm.morantest(model2, nb.FOQ.lw.W)

# Contrast Moran's I of dependent variables and models
mt.y1.nb.FOQ
mt.lmres.y1.nb.FOQ
mt.y2.nb.FOQ
mt.lmres.y2.nb.FOQ

# Contrast Moran's I scatterplots of dependent variables and models
par(mfrow=c(2,2))
moran.plot(SNM05, nb.FOQ.lw.W, main="Share Nonmarital Births 2005",
xlab="Share Nonmarital Births",ylab="Spatially Lagged Share Nonmarital Births",
label=GEN)
moran.plot(model1$residuals, nb.FOQ.lw.W, 
main="Model Residuals SNM ~ Unemployment Rate", xlab="Residuals",
ylab="Spatially Lagged Residuals",label=GEN)
moran.plot(TFR05, nb.FOQ.lw.W, main="Total Fertility Rate 2005",
xlab="Total Fertility Rate",
ylab="Spatially Lagged Share Total Fertility Rate", label=GEN)
moran.plot(model2$residuals, nb.FOQ.lw.W, 
main="Model Residuals TFR ~ Unemployment Rate", xlab="Residuals",
ylab="Spatially Lagged Residuals",label=GEN)

# Contrast maps of dependent variables and maps of residuals
var.data <- data.frame(SNM05, TFR05,model1$residuals,model2$residuals)
var.names <-c("Share Nonmarital Births 2005", "TFR 2005", 
"Residuals Nonmarital Births Model", "Residuals TFR Model")
m <- length(var.data)
colors <- matrix(nrow=n,ncol=m, data=rep(0,n*m))

par(mfrow=c(2,2))
for (i in 1:m) {
varofint <- var.data[,i]
lb       = 5
mean <- mean(varofint)
sd <- sd(varofint)
min <- min(varofint)          
max <- max(varofint)
#if (i==1) 
bins = c(min,mean-sd, mean, mean+sd, max)
colpal   = brewer.pal(lb-1, "YlOrRd")
colors[,i]   = colpal[findInterval(varofint, bins, rightmost.closed=T)]
plot(shape.shp, col=colors[,i])
title(var.names[[i]], cex=1)
legend("bottomright",fill=colpal,legend=paste(round(bins[-length(bins)],2)
,"-",round(bins[-1],2)),cex=1, bg="white")
}

# Look at Cook's distance - Many outliers are situated in the former
# German-German border region
a <- moran.plot(SNM05,nb.FOQ.lw.W)
var.data <- data.frame(round(a$infmat[,5],6))
var.names <-c("Share Nonmarital Births 2005 \n Cook's Distance")
m <- length(var.data)
colors <- matrix(nrow=n,ncol=m, data=rep(0,n*m))

par(mfrow=c(1,1))
for (i in 1:m) {
varofint <- var.data[,i]
lb       = 10
mean <- mean(varofint)
sd <- sd(varofint)
min <- min(varofint)
max <- max(varofint)
bins     = seq(min(varofint), max(varofint), length=lb)
colpal   = brewer.pal(lb-1, "YlOrRd")
colors[,i]   = colpal[findInterval(varofint, bins, rightmost.closed=T)]
plot(shape.shp, col=colors[,i])
title(var.names[[i]], cex=1)
legend("bottomright",fill=colpal,legend=paste(round(bins[-length(bins)],2)
,"-",round(bins[-1],2)),cex=1, bg="white")
}
