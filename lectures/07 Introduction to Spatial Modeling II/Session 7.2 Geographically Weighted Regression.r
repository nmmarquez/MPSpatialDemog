################################################################################
#                                                                              #                          
# IDEM 156 - Spatial Demography Course 2018                                    #                          
# Geographically Weighted Regression                                           #                          
# Sebastian Kl?sener, MPIDR                                                    #                          
#                                                                              #                          
################################################################################

# Erase all objects in workspace
rm(list=ls(all=TRUE))

# Load libraries
library(spgwr)
library(RColorBrewer)
library(spdep)
library(rgdal)
library(maptools)


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

# Derive variable names
names(shape.shp)

# In order to "attach" data from a shapefile attribute table, we have to call it
# with "@data"
attach(shape.shp@data)
            
# Define position of east and west districts in table
west.dist <- c(1:326)
east.dist <-  c(327:439)
                 
# Create neighborhood weight matrices
# First order queen contiguity
nb.FOQ <- poly2nb(shape.shp, queen=TRUE)
summary(nb.FOQ)

# Code to connect R?gen
#edit.nb(nb.FOQ, polys=shape.shp)
nb.FOQ.cor <- nb.FOQ
nb.FOQ.cor[[362]]  <- as.integer(c(350, 358))
nb.FOQ.cor[[350]]  <- as.integer(c(358, 362))
nb.FOQ.cor[[358]]  <- as.integer(c(350,352,353,354,360,362))
lw.FOQ <- nb2listw(nb.FOQ.cor) 


################################################################################
# 2) Create sub-datasets of West and East Germany                              #
################################################################################

shape.west <-  shape.shp[c(1:326),]
plot(shape.west,col="seagreen3")
shape.east <-  shape.shp[c(327:439),]
plot(shape.east,col="orange",add=T)

nb.FOQw <- poly2nb(shape.west, queen=TRUE)
lw.FOQw <- nb2listw(nb.FOQw)
nb.FOQe <- poly2nb(shape.east, queen=TRUE)

# Conntect R?gen
coordse <- coordinates(shape.east)               
ne <- length(shape.east)                   
vec <- c(1:ne)
#plot.nb(nb.FOQe, coordse, add=TRUE) 
#text(coordse, paste(vec), col="red")
nb.FOQe.cor <- nb.FOQe
nb.FOQe.cor[[36]]  <- as.integer(c(24,32))
nb.FOQe.cor[[24]]  <- as.integer(c(32,36))
nb.FOQe.cor[[32]]  <- as.integer(c(24,26,27,28,34,36))
lw.FOQe <- nb2listw(nb.FOQe.cor)


################################################################################
# 3)   We again deal with the Simpson-Paradoxon, that both in East and West    #
#      Germany there is at the district level a negative association between   #
#      Share Nonmarital Births and Female Labour Force Participation, but, if  #
#      if we look at all German districts together, the association is         #
#      positive. We check, how the Geographically Weighted Regression performs #
#      on such our data in comparison to a linear model.                       #
# 3.1) Define Variables                                                        #
################################################################################
    
# Defining Variables      
# Dependent Variable: Share Non-Marital Births
y <- SNM05
# Covariate 1: Female Employment Rate
x1 <- FR155005                        
# Covariate 2: Unemployment Rate
x2 <- UER05
# Covariate 3: Income per Capita
x3 <- INCOME05/1000 
# Covatiate 4: Tourism
x4 <- log(TOUR05)
# Dummy East
dum_east <- DUM_E                            
                         
                         
################################################################################
# 3.2) Plot of our Simpson-Paradoxon                                           #
################################################################################                                                                                                   

# Sample West
ywest <- y[west.dist]
x1west <- x1[west.dist]

# Sample East
yeast <- y[east.dist]
x1east <- x1[east.dist]

# Scatter plot of y and x1
plot(x1,y, xlab="Female Employment Rate 2005", 
 ylab="Share Nonmarital Births 2005", 
 main="Share Nonmarital Births - Female Employment Rate")
# Points for Sample West and Sample East in different colors
points(x1west,ywest, col="orange2")
points(x1east,yeast, col="navyblue")
abline(coef(lm(y~x1)), lwd="1")

# Lines for Sample West and Sample East
abline(coef(lm(ywest~x1west)), col="red", lty=2)
abline(coef(lm(yeast~x1east)), col="red", lty=2)
legend("topleft", c("East","West"), pch=1, col=c("navyblue","orange2"), 
 bg="white")
                         
# Map of y (Share Non-marital Births 05) and x1 (Female Employment Rate)
par(mfrow=c(1,2))
var.data <- data.frame(y, x1)
var.names <-c("Share Nonmarital Births 2005", "Female Employment Rate 2005")
m <- length(var.data)
for (i in 1:m) {
    colors <- matrix(nrow=n,ncol=m, data=rep(0,n*m))
    varofint <- var.data[,i]
    mean <- mean(varofint)
    sd <- sd(varofint)
    min <- min(varofint)          
    max <- max(varofint)
    bins = c(min,mean-sd, mean-sd*0.5, mean, mean+sd*0.5, mean+sd, max)
    lb <- length(bins)
    med <- median(1:lb)
    colpal <- brewer.pal(length(bins), "PRGn")[-med]
    colors <- colpal[findInterval(varofint, bins, rightmost.closed=T)]
    plot(shape.shp, col=colors)
    title(var.names[[i]], cex=1)
    legend("bottomright",fill=colpal,
    legend=paste(round(bins[-length(bins)],2),"-",round(bins[-1],2)),cex=0.6,
    bg="white")
}


################################################################################
# 4)   Calculate GWR-Model                                                     #
# 4.1) Simple Model - Share Nonmarital Births ~ Female Employment Rate         #
################################################################################                                                                          

# Simple OLS-Model for comparison
model.west.and.east <- lm(y ~ x1)
summary(model.west.and.east)
lm.morantest(model.west.and.east,lw.FOQ)

# A simple "Geographically Weighted Regression", where we split up the sample
# into two regions and investigate, whether the estimates of x1 (in our case
# Female Employment Rate) differ.

# Running Model for West Germany
model.west <- lm(y[west.dist] ~ x1[west.dist])
# Running Model for East Germany
model.east <- lm(y[east.dist] ~ x1[east.dist])

# Summary and Moran's I test for the seperate models
summary(model.west)
summary(model.east)
lm.morantest(model.west,lw.FOQw)
lm.morantest(model.east,lw.FOQe)

# GWR-bandwith estimation based on given geographical co-ordinates of our n
# locations. In our case n is equal to 439 districts
# Code to calculate a global bandwith, which determines, which nearby neighbors
# j are included in the local regression. 
ger.bw.gl <- gwr.sel(y ~ x1, coords=coords, longlat=TRUE)

# An alternative way, where AIC is minimized
ger.bw.glaic <- gwr.sel(y ~ x1, coords=coords, method="AIC", longlat=TRUE)
#ger.bw.glaic <- 29.18101

# Code for adaptive bandwith, in which the number of neighboring observations j, 
# which are included in each local regression, are the same for each location. 
ger.bw.ad <- gwr.sel(y ~ x1,coords=coords, adapt=TRUE, longlat=TRUE)

# Calculating Geographically weighted regression
gwr.gl <-    gwr(y ~ x1, coords=coords, bandwidth=ger.bw.gl, longlat=TRUE)
gwr.glaic <- gwr(y ~ x1, coords=coords, bandwidth=ger.bw.glaic, longlat=TRUE)
gwr.ad  <-   gwr(y ~ x1, coords=coords, adapt=ger.bw.ad, longlat=TRUE)
gwr.gl
gwr.glaic
gwr.ad

# Obtained values are sitting in $SDF. 
gwr.gl$SDF[1:6,]

# Boxplot of the local beta estimates for Female Employment Rate
par(mfrow=c(1,1))
boxplot(gwr.gl$SDF$x1,gwr.glaic$SDF$x1,gwr.ad$SDF$x1, 
 names=c("GWR.GL.CV","GWR.GL.AIC","GWR.ADJ"), col="grey", 
 main="Local Beta-Estimates for Female Employment Rate")

# Add coefficients obtained in linear model
abline(h = model.west.and.east$coefficients[2], col="red", lwd="2")
text(0.6,model.west.and.east$coefficients[2]-0.3, "Beta-\nestimate\nlm-model", 
 col="red", cex=0.7)

# Map the Results - Contrast the beta estimates for the different bandwidth
par(mfrow=c(2,2))
var.data <- data.frame(y, gwr.gl$SDF$x1,gwr.glaic$SDF$x1,gwr.ad$SDF$x1 )
var.names <-c("Share Nonmarital Births 2005", 
  "Beta Female Emp. Rate GWR.GL.CV", "Beta Female Emp. Rate GWR.GL.AIC", 
  "Beta Female Emp. Rate GWR.ADF") 
m <- length(var.data)
colors <- matrix(nrow=n,ncol=m, data=rep(0,n*m))
for (i in 1:m) {
    varofint <- var.data[,i]
    mean <- mean(varofint)
    sd <- sd(varofint)
    min <- min(varofint)          
    max <- max(varofint)
    bins = c(min,mean-sd, mean-sd*0.5, mean, mean+sd*0.5, mean+sd, max)
    lb <- length(bins)
    med <- median(1:lb)
    colpal <- brewer.pal(length(bins), "PRGn")[-med]
    colors <- colpal[findInterval(varofint, bins, rightmost.closed=T)]
    plot(shape.shp, col=colors, bg="grey5",lty=0)
    title(var.names[[i]], cex=1)
    legend("bottomright",fill=colpal,
    legend=paste(round(bins[-length(bins)],2),"-",round(bins[-1],2)),cex=0.7,
    bg="white")
}

# Alternative Map-series, where we map all variables as well as the beta-
# estimate and the local R2
par(mfrow=c(2,2))
var.data <- data.frame(y, x1, gwr.gl$SDF$x1, gwr.gl$SDF$localR2)
var.names <-c("Share Nonmarital Births 2005", "Female Employment Rate 2005", 
 "Beta Female Emp. Rate 2005", "Local R2")
m <- length(var.data)
colors <- matrix(nrow=n,ncol=m, data=rep(0,n*m))
for (i in 1:m) {
    varofint <- var.data[,i]
    mean <- mean(varofint)
    sd <- sd(varofint)
    min <- min(varofint)          
    max <- max(varofint)
    bins = c(min,mean-sd, mean-sd*0.5, mean, mean+sd*0.5, mean+sd, max)
    lb <- length(bins)
    med <- median(1:lb)
    colpal <- brewer.pal(length(bins), "PRGn")[-med]
    colors <- colpal[findInterval(varofint, bins, rightmost.closed=T)]
    plot(shape.shp, col=colors, bg="grey5",lty=0)
    title(var.names[[i]], cex=1)
    legend("bottomright",fill=colpal,
    legend=paste(round(bins[-length(bins)],2),"-",round(bins[-1],2)),cex=0.7,
    bg="white")
}


################################################################################                                                                          
# 4.2) Second Model -                                                          #
#      Share Nonmarital Births ~ Female Employment Rate+ Dummy East            #
################################################################################                                                                          

# What happens if we control in addition for Dummy East
# As we adapt our mode, we need to determine new bandwidths
model.west.and.east2 <- lm(y ~ x1+dum_east)
summary(model.west.and.east)
summary(model.west.and.east2)

# As we adapt our mode, we need to determine new bandwidths
ger.bw2.gl <- gwr.sel(y ~ x1+dum_east, coords=coords, longlat=TRUE)
#ger.bw2.gl <- 33.1454
ger.bw2.ad <- gwr.sel(y ~ x1+dum_east,coords=coords, adapt=TRUE, longlat=TRUE)
#ger.bw2.ad <- 0.009530558

# How did the bandwidth change
tit <- c("ger.bw.gl","ger.bw2.gl","ger.bw.ad","ger.bw2.ad")
values <- round(c(ger.bw.gl,ger.bw2.gl,ger.bw.ad,ger.bw2.ad),3)
compare <- data.frame(values)
row.names(compare) <- tit 
compare
                                                                             
gwr.gl2 <- gwr(y ~ x1+dum_east, coords=coords, bandwidth=ger.bw2.gl, longlat=T)
gwr.ad2 <- gwr(y ~ x1+dum_east, coords=coords, adapt=ger.bw2.ad, longlat=T)

# Boxplot of beta-estimates for female employment rate
par(mfrow=c(1,2))
bp <- boxplot(gwr.gl$SDF$x1,gwr.ad$SDF$x1, ylim=c(-3,3),
       names=c("GWR.GL - y~x1","GWR.AD - y~x1"), col="grey", 
       main="Local Beta-Estimates for Female Emp. Rate")
# Coefficient obtained in linear model
abline(h = model.west.and.east$coefficients[2], col="red", lwd="2")
text(0.55,model.west.and.east$coefficients[2]-0.3, "Beta-\nestimate\nlm-model", 
 col="red", cex=0.7, ylim=c(-3,3))

boxplot(gwr.gl2$SDF$x1,gwr.ad2$SDF$x1, ylim=c(-3,3),
 names=c("GWR.GL - y~x1+dE","GWR.AD - y~x1+dE"), 
 col="grey", main="Local Beta-Estimates for Female Emp. Rate")
# Coefficient obtained in linear model
abline(h = model.west.and.east2$coefficients[2], col="red", lwd="2")
text(0.55,model.west.and.east2$coefficients[2]+0.55, 
 "Beta-\nestimate\nlm-model", col="red", cex=0.7)

# Boxplot of beta estimates for Dummy East
par(mfrow=c(1,1))
bp <- boxplot(gwr.gl2$SDF$dum_east, gwr.ad2$SDF$dum_east,
       names=c("GWR.GL.CV - y~x1+dE","GWR.ADJ - y~x1+dE"), col="grey", 
       main="Local Beta-Estimates for Dummy East")
# Coefficient obtained in linear model
abline(h = model.west.and.east2$coefficients[3], col="red", lwd="2")
text(0.55,model.west.and.east2$coefficients[3]+2.2, "Beta-\nestimate\nlm-model", 
 col="red", cex=0.7, ylim=c(-3,3))

# Map of y, beta estimates of female labour force participation and dummy east, 
# and local r-squared 
par(mfrow=c(2,2))
var.data <- data.frame(y, gwr.gl2$SDF$x1, gwr.gl2$SDF$dum_east, 
 gwr.gl2$SDF$localR2)
var.names <-c("Share Nonmarital Births 2005", "Beta Female Emp. Rate 2005",
 "Beta Dummy East 2005", "Local R2")
m <- length(var.data)
colors <- matrix(nrow=n,ncol=m, data=rep(0,n*m))
for (i in 1:m) {
    varofint <- var.data[,i]
    mean <- mean(varofint)
    sd <- sd(varofint)
    min <- min(varofint)          
    max <- max(varofint)
    bins = c(min,mean-sd, mean-sd*0.5, mean, mean+sd*0.5, mean+sd, max)
    lb <- length(bins)
    med <- median(1:lb)
    colpal <- brewer.pal(length(bins), "PRGn")[-med]
    colors <- colpal[findInterval(varofint, bins, rightmost.closed=T)]
    plot(shape.shp, col=colors, bg="grey5",lty=0)
    title(var.names[[i]], cex=1)
    legend("bottomright",fill=colpal,
    legend=paste(round(bins[-length(bins)],2),"-",round(bins[-1],2)),cex=0.7,
    bg="white")
}

################################################################################                                                                          
# 4.3) Third Model -                                                           #
#      Share Nonmarital Births ~ Female Employment Rate+Unemployment Rate+     #
#                                Income+Tourism                                #
################################################################################  

# Model with four variables
model.west.and.east3 <- lm(y ~ x1+x2+x3+x4)

ger.bw3.gl <- gwr.sel(y ~ x1+x2+x3+x4, coords=coords, longlat=TRUE)
#ger.bw3.gl <- 38.94795
gwr.gl3 <- gwr(y ~ x1+x2+x3+x4, coords=coords, bandwidth=ger.bw3.gl, longlat=T)

# Boxplots of beta estimates
par(mfrow=c(2,2))
bp <- boxplot(gwr.gl3$SDF$x1, col="grey", 
       main="Loc. Beta-Estimates for Fem. Emp. Rate", cex.main=1)    
abline(h = model.west.and.east3$coefficients[2], col="red", lwd="2")
bp <- boxplot(gwr.gl3$SDF$x2, col="grey", 
       main="Loc. Beta-Estimates for Unemp. Rate", cex.main=1)    
abline(h = model.west.and.east3$coefficients[3], col="red", lwd="2")
bp <- boxplot(gwr.gl3$SDF$x3, col="grey", 
       main="Local Beta-Estimates for Income", cex.main=1)
abline(h = model.west.and.east3$coefficients[4], col="red", lwd="2")
bp <- boxplot(gwr.gl3$SDF$x4, col="grey", 
       main="Local Beta-Estimates for Tourism", cex.main=1)
abline(h = model.west.and.east3$coefficients[5], col="red", lwd="2")

par(mfrow=c(2,2))
var.data <- data.frame(gwr.gl3$SDF$x1, gwr.gl3$SDF$x2, gwr.gl3$SDF$x3,
 gwr.gl3$SDF$x4)
var.names <-c("Beta Female Employment Rate 2005", "Beta Unemployment Rate 2005",
 "Beta Income 2005", "Beta Tourism 2005")

m <- length(var.data)
colors <- matrix(nrow=n,ncol=m, data=rep(0,n*m))
for (i in 1:m) {
    varofint <- var.data[,i]
    mean <- mean(varofint)
    sd <- sd(varofint)
    min <- min(varofint)          
    max <- max(varofint)
    bins = c(min,mean-sd, mean-sd*0.5, mean, mean+sd*0.5, mean+sd, max)
    lb <- length(bins)
    med <- median(1:lb)
    colpal <- brewer.pal(length(bins), "PRGn")[-med]
    colors <- colpal[findInterval(varofint, bins, rightmost.closed=T)]
    plot(shape.shp, col=colors, bg="grey5",lty=0)
    title(var.names[[i]], cex=1)
    legend("bottomright",fill=colpal,
    legend=paste(round(bins[-length(bins)],2),"-",round(bins[-1],2)),cex=0.7,
    bg="white")
}
