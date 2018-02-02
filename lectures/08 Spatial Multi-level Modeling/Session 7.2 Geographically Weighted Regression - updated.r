################################################################################
#                                                                              #                          
# IDEM 156 - Spatial Demography Course 2018                                    #                          
# Geographically Weighted Regression                                           #                          
# Sebastian Klüsener, MPIDR                                                    #                          
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

# Set your working directory
# If you need to change the working directory, you can use the following code:
# Set drive
main <- c("C")
# Set path to session folders
path <- ":/ownCloud/01 MPI/120 Spatial Demography 2018"
# Define session folder
path2 <- "/02 Sessions/07 Introduction to Spatial Modeling II"
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

# Coordinates longlat-projection
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

# Code to connect Rügen
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

# Conntect Rügen
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
# locations. In our case n is equal to 439 districts.
# Code to calculate a global bandwith, which determines, in which areas around
# region i neighbors are included in the regression
ger.bw.gl <- gwr.sel(y ~ x1, coords=coords, longlat=TRUE)

# An alternative way, where AIC is minimized
ger.bw.glaic <- gwr.sel(y ~ x1, coords=coords, method="AIC", longlat=TRUE)

# Code for an adaptive bandwith, in which the number of neighboring observations j, 
# which are included in each local regression, are the same for each location. 
ger.bw.ad <- gwr.sel(y ~ x1,coords=coords, adapt=TRUE, longlat=T)
ger.bw.adaic <- gwr.sel(y ~ x1,coords=coords, method="AIC",adapt=TRUE, longlat=T)

# Calculating Geographically weighted regression
# With longlat coordinates
gwr.gl <-    gwr(y ~ x1, coords=coords, bandwidth=ger.bw.gl, longlat=TRUE)
gwr.glaic <- gwr(y ~ x1, coords=coords, bandwidth=ger.bw.glaic, longlat=TRUE)
gwr.ad  <-   gwr(y ~ x1, coords=coords, adapt=ger.bw.ad, longlat=TRUE)
gwr.adaic  <-   gwr(y ~ x1, coords=coords, adapt=ger.bw.adaic, longlat=TRUE)

# Running a gwr with an adaptive bandwidth of 1 (all observations included) 
# returns almost the estimates of the global model.
gwr.ad  <-   gwr(y ~ x1, coords=coords, adapt=1, longlat=TRUE)

gwr.gl
gwr.glaic
gwr.ad
gwr.adaic

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
# As we adapt our model, we need to determine new bandwidths
model.west.and.east2 <- lm(y ~ x1+dum_east)
summary(model.west.and.east)
summary(model.west.and.east2)

# As we adapt our model, we need to determine new bandwidths
ger.bw2.gl <- gwr.sel(y ~ x1+dum_east, coords=coords, longlat=TRUE)

ger.bw2.ad <- gwr.sel(y ~ x1+dum_east,coords=coords, adapt=TRUE, longlat=TRUE)

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

################################################################################                                                                          
# 4.4) Chech adaptive vs. global bandwidth                                     #
################################################################################  

# Modified gwr function allowing to extract the weights by distance
gwr_ch <- function (formula, data = list(), coords, bandwidth, gweight = gwr.Gauss, 
                    adapt = NULL, hatmatrix = FALSE, fit.points, longlat = NULL, 
                    se.fit = FALSE, weights, cl = NULL, predictions = FALSE, 
                    fittedGWRobject = NULL, se.fit.CCT = TRUE) 
{
  timings <- list()
  .ptime_start <- proc.time()
  this.call <- match.call()
  p4s <- as.character(NA)
  Polys <- NULL
  if (is(data, "SpatialPolygonsDataFrame")) 
    Polys <- as(data, "SpatialPolygons")
  if (is(data, "Spatial")) {
    if (!missing(coords)) 
      warning("data is Spatial* object, ignoring coords argument")
    coords <- coordinates(data)
    p4s <- proj4string(data)
    if (is.null(longlat) || !is.logical(longlat)) {
      if (!is.na(is.projected(data)) && !is.projected(data)) {
        longlat <- TRUE
      }
      else {
        longlat <- FALSE
      }
    }
    data <- as(data, "data.frame")
  }
  if (is.null(longlat) || !is.logical(longlat)) 
    longlat <- FALSE
  if (missing(coords)) 
    stop("Observation coordinates have to be given")
  if (is.null(colnames(coords))) 
    colnames(coords) <- c("coord.x", "coord.y")
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "weights"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  dp.n <- length(model.extract(mf, "response"))
  weights <- as.vector(model.extract(mf, "weights"))
  if (!is.null(weights) && !is.numeric(weights)) 
    stop("'weights' must be a numeric vector")
  if (is.null(weights)) 
    weights <- rep(as.numeric(1), dp.n)
  if (any(is.na(weights))) 
    stop("NAs in weights")
  if (any(weights < 0)) 
    stop("negative weights")
  y <- model.extract(mf, "response")
  x <- model.matrix(mt, mf)
  lm <- lm.wfit(x, y, w = weights)
  lm$x <- x
  lm$y <- y
  gTSS <- c(cov.wt(matrix(y, ncol = 1), wt = weights, method = "ML")$cov * 
              dp.n)
  if (hatmatrix) 
    se.fit <- TRUE
  if (hatmatrix) 
    predictions <- TRUE
  if (missing(fit.points)) {
    fp.given <- FALSE
    fittedGWRobject <- NULL
    predictions <- TRUE
    fit.points <- coords
    colnames(fit.points) <- colnames(coords)
    predx <- x
  }
  else fp.given <- TRUE
  griddedObj <- FALSE
  if (is(fit.points, "Spatial")) {
    if (predictions) {
      t1 <- try(slot(fit.points, "data"), silent = TRUE)
      if (class(t1) == "try-error") 
        stop("No data slot in fit.points")
      predx <- try(model.matrix(delete.response(mt), fit.points))
      if (class(predx) == "try-error") 
        stop("missing RHS variable in fit.points")
      if (ncol(predx) != ncol(x)) 
        stop("new data matrix columns mismatch")
    }
    Polys <- NULL
    if (is(fit.points, "SpatialPolygonsDataFrame")) {
      Polys <- as(fit.points, "SpatialPolygons")
      fit.points <- coordinates(fit.points)
    }
    else {
      griddedObj <- gridded(fit.points)
      fit.points <- coordinates(fit.points)
    }
  }
  else {
    if (predictions && fp.given) 
      stop("predictions not available for matrix fit points")
  }
  n <- NROW(fit.points)
  rownames(fit.points) <- NULL
  if (is.null(colnames(fit.points))) 
    colnames(fit.points) <- c("x", "y")
  if (predictions) {
    if (nrow(predx) != nrow(fit.points)) 
      stop("new data matrix rows mismatch")
    fit.points <- cbind(fit.points, predx)
  }
  fit_are_data <- isTRUE(all.equal(fit.points, coords, check.attributes = FALSE))
  input_predictions <- predictions
  if (fit_are_data && !predictions) {
    predictions <- TRUE
    input_predictions <- FALSE
    predx <- x
    fit.points <- cbind(fit.points, predx)
  }
  m <- NCOL(x)
  if (NROW(x) != NROW(coords)) 
    stop("Input data and coordinates have different dimensions")
  if (missing(bandwidth) && is.null(adapt)) 
    stop("Bandwidth must be given for non-adaptive weights")
  if (!is.null(adapt)) {
    stopifnot(is.numeric(adapt))
    stopifnot((adapt >= 0))
    stopifnot((adapt <= 1))
  }
  else {
    stopifnot(length(bandwidth) == 1)
  }
  if (missing(bandwidth)) 
    bandwidth <- NULL
  lhat <- NA
  yhat <- NULL
  if (!is.null(fittedGWRobject)) {
    yhat <- fittedGWRobject$SDF$pred
  }
  GWR_args <- list(fp.given = fp.given, hatmatrix = hatmatrix, 
                   longlat = longlat, bandwidth = bandwidth, adapt = adapt, 
                   se.fit = se.fit, predictions = predictions, se.fit.CCT = se.fit.CCT, 
                   fit_are_data = fit_are_data)
  timings[["set_up"]] <- proc.time() - .ptime_start
  .ptime_start <- proc.time()
  if (!is.null(cl) && length(cl) > 1 && fp.given && !hatmatrix) {
    if (requireNamespace("parallel", quietly = TRUE)) {
      l_fp <- lapply(parallel::splitIndices(nrow(fit.points), 
                                            length(cl)), function(i) fit.points[i, , drop = FALSE])
      parallel::clusterEvalQ(cl, library(spgwr))
      varlist <- list("GWR_args", "coords", "gweight", 
                      "y", "x", "weights", "yhat")
      env <- new.env()
      assign("GWR_args", GWR_args, envir = env)
      assign("coords", coords, envir = env)
      assign("gweight", gweight, envir = env)
      assign("y", y, envir = env)
      assign("x", x, envir = env)
      assign("weights", weights, envir = env)
      assign("yhat", yhat, envir = env)
      parallel::clusterExport(cl, varlist, env)
      res <- parallel::parLapply(cl, l_fp, function(fp) .GWR_int(fit.points = fp, 
                                                                 coords = coords, gweight = gweight, y = y, x = x, 
                                                                 weights = weights, yhat = yhat, GWR_args = GWR_args))
      parallel::clusterEvalQ(cl, rm(varlist))
      rm(env)
      df <- list()
      df$df <- as.data.frame(do.call("rbind", lapply(res, 
                                                     function(x) x$df)))
      bw <- do.call("c", lapply(res, function(x) x$bw))
      results <- NULL
    }
    else {
      stop("parallel not available")
    }
  }
  else {
    df <- .GWR_int(fit.points = fit.points, coords = coords, 
                   gweight = gweight, y = y, x = x, weights = weights, 
                   yhat = yhat, GWR_args = GWR_args)
    if (!fp.given && hatmatrix) 
      lhat <- df$lhat
    bw <- df$bw
    results <- NULL
  }
  timings[["run_gwr"]] <- proc.time() - .ptime_start
  .ptime_start <- proc.time()
  if (predictions && !input_predictions) 
    predictions <- FALSE
  if (!fp.given && hatmatrix) {
    v1 <- sum(diag(lhat))
    B2 <- t(lhat) %*% lhat
    v2 <- sum(diag(B2))
    edf <- dp.n - 2 * v1 + v2
    B1 <- t(diag(dp.n) - lhat) %*% (diag(dp.n) - lhat)
    rss <- c(t(y) %*% B1 %*% y)
    delta1 <- sum(diag(B1))
    sigma2 <- rss/delta1
    odelta2 <- sum(diag(B1)^2)
    delta2 <- sum(diag(B1 %*% B1))
    nu1 <- sum(diag(B2))
    sigma2.b <- rss/dp.n
    AICb.b <- 2 * dp.n * log(sqrt(sigma2.b)) + dp.n * log(2 * 
                                                            pi) + (dp.n * ((n + v1)/(dp.n - 2 - v1)))
    AICh.b <- 2 * dp.n * log(sqrt(sigma2.b)) + dp.n * log(2 * 
                                                            pi) + dp.n + v1
    AICc.b <- 2 * dp.n * log(sqrt(sigma2.b)) + dp.n * log(2 * 
                                                            pi) + dp.n * ((delta1/delta2) * (dp.n + nu1))/((delta1^2/delta2) - 
                                                                                                             2)
    results <- list(v1 = v1, v2 = v2, delta1 = delta1, delta2 = delta2, 
                    sigma2 = sigma2, sigma2.b = sigma2.b, AICb = AICb.b, 
                    AICh = AICh.b, AICc = AICc.b, edf = edf, rss = rss, 
                    nu1 = nu1, odelta2 = odelta2, n = dp.n)
    timings[["postprocess_hatmatrix"]] <- proc.time() - .ptime_start
    .ptime_start <- proc.time()
  }
  if ((!fp.given || fit_are_data) && is.null(fittedGWRobject)) {
    localR2 <- numeric(n)
    if (is.null(adapt)) {
      bw <- bandwidth
      bandwidthR2 <- rep(bandwidth, n)
    }
    else {
      bandwidthR2 <- gw.adapt(dp = coords, fp = fit.points[, 
                                                           1:2, drop = FALSE], quant = adapt, longlat = longlat)
      bw <- bandwidthR2
    }
    if (any(bandwidth < 0)) 
      stop("Invalid bandwidth")

# Add object to store weights and distances
    list_weight <- list()
    list_dist <- list()
    for (i in 1:n) {
      dxs <- spDistsN1(coords, fit.points[i, 1:2], longlat = GWR_args$longlat)
      if (any(!is.finite(dxs))) 
        dxs[which(!is.finite(dxs))] <- .Machine$double.xmax/2
      w.i <- gweight(dxs^2, bandwidthR2[i])
      w.i <- w.i * weights
# Store weight and distances
      list_weight[[i]] <- w.i
      list_dist[[i]] <- dxs
      if (any(w.i < 0 | is.na(w.i))) 
        stop(paste("Invalid weights for i:", i))
      RSS <- sum(w.i * (y - df$df[, "pred"])^2)
      yss <- sum(w.i * (y - weighted.mean(y, w.i))^2)
      localR2[i] <- 1 - (RSS/yss)
    }
    df$df <- cbind(df$df, localR2)
    timings[["postprocess_localR2"]] <- proc.time() - .ptime_start
    .ptime_start <- proc.time()
  }
  if (se.fit) {
    EDFS <- NULL
    normSigmaS <- NULL
    EDF <- NULL
    normSigma <- NULL
    if (fp.given && !is.null(fittedGWRobject)) {
      if (fittedGWRobject$hatmatrix) {
        EDF <- fittedGWRobject$results$edf
        normSigma <- sqrt(fittedGWRobject$results$rss/EDF)
        EDFS <- fittedGWRobject$results$n - fittedGWRobject$results$v1
        normSigmaS <- sqrt(fittedGWRobject$results$rss/EDFS)
      }
    }
    if (!fp.given && hatmatrix) {
      EDFS <- results$n - results$v1
      normSigmaS <- sqrt(results$rss/EDFS)
      EDF <- results$edf
      normSigma <- sqrt(results$rss/EDF)
    }
    ses <- grep("_se", colnames(df$df))
    senms <- colnames(df$df)[ses]
    betase <- df$df[, ses, drop = FALSE]
    df$df[, ses] <- NA
    if (predictions) {
      pred.se <- df$df[, "pred.se", drop = FALSE]
      df$df[, "pred.se"] <- NA
    }
    if (!is.null(EDF)) {
      betaseEDF <- normSigma * sqrt(betase)
      colnames(betaseEDF) <- paste(senms, "EDF", sep = "_")
      df$df[, ses] <- normSigmaS * sqrt(betase)
      df$df <- cbind(df$df, betaseEDF)
      if (predictions) {
        pred.se_EDF <- normSigma * sqrt(pred.se)
        df$df[, "pred.se"] <- normSigmaS * sqrt(pred.se)
        df$df <- cbind(df$df, pred.se_EDF)
      }
    }
    else {
      warning("standard errors set to NA, normalised RSS not available")
    }
    timings[["postprocess_SE"]] <- proc.time() - .ptime_start
    .ptime_start <- proc.time()
  }
  df <- as.data.frame(df$df)
  if (predictions) 
    fit.points <- fit.points[, 1:2, drop = FALSE]
  row.names(fit.points) <- row.names(df)
  SDF <- SpatialPointsDataFrame(coords = fit.points, data = df, 
                                proj4string = CRS(p4s))
  if (griddedObj) {
    gridded(SDF) <- TRUE
  }
  else {
    if (!is.null(Polys)) {
      df <- data.frame(SDF@data)
      rownames(df) <- sapply(slot(Polys, "polygons"), function(i) slot(i, 
                                                                       "ID"))
      SDF <- SpatialPolygonsDataFrame(Sr = Polys, data = df)
    }
  }
  timings[["final_postprocess"]] <- proc.time() - .ptime_start
  z <- list(SDF = SDF, lhat = lhat, lm = lm, results = results, 
            bandwidth = bw, adapt = adapt, hatmatrix = hatmatrix, 
            gweight = deparse(substitute(gweight)), gTSS = gTSS, 
            this.call = this.call, fp.given = fp.given, timings = do.call("rbind", 
                                                                          timings)[, c(1, 3)])
  class(z) <- "gwr"
  invisible(z)
# List to store weights and distances
  retlist <- list()
  retlist[[1]] <- list_weight
  retlist[[2]] <- list_dist
  return(retlist)
}
# Important to link the new adjusted funtion to the namespace of the library spgwr
# so that it has access to hidden functions.
environment(gwr_ch) <- asNamespace('spgwr')

# Adaptive weight is indeed adaptive, but includes more then 2 regions.
ger.bw.ad*439

# 2-dimensional plot of the weighting kernel 
par(mfrow=c(2,2))
gwr.ad  <-   gwr_ch(y ~ x1, coords=coords, adapt=ger.bw.ad, longlat=TRUE)

plot(gwr.ad[[2]][[1]],round(gwr.ad[[1]][[1]],4),main="Region 1 - adaptive",
     ylab="Weight for region j",xlab="Distance to region 1 in km")
abline(v=ger.bw.gl,col="orangered1")
plot(gwr.ad[[2]][[2]],round(gwr.ad[[1]][[2]],4),main="Region 2 - adaptive",
     ylab="Weight for region j",xlab="Distance to region 2 in km")
abline(v=ger.bw.gl,col="orangered1")

# Also fixed bandwidht gives regions beyong the bandwidht of ger.bw.gl still weight 
gwr.gl  <-   gwr_ch(y ~ x1, coords=coords, bandwidth=ger.bw.gl, longlat=TRUE)
gwr.gl  <-   gwr_ch(y ~ x1, coords=coords, bandwidth=100, longlat=TRUE)

plot(gwr.gl[[2]][[1]],round(gwr.gl[[1]][[1]],4),main="Region 1 - fixed (100km)",
     ylab="Weight for region j",xlab="Distance to region 1 in km")
abline(v=100,col="orangered1")
abline(h=0.5,col="orangered1")
plot(gwr.gl[[2]][[2]],round(gwr.gl[[1]][[2]],4),main="Region 2 - fixed (100km)",
     ylab="Weight for region j",xlab="Distance to region 2 in km")
abline(v=100,col="orangered1")
abline(h=0.5,col="orangered1")


# Plot of distances between districts (considers distrance from  regions a to b,
# as well as distance from b to a)
dist <- unlist(gwr.gl[[2]])

# Erase all distance 
dist <- dist[dist>0]
  
# Density plot of the distances
par(mfrow=c(1,1))
plot(density(unlist(gwr.gl[[2]])),main="Distances between the 439 German districts",xlab="Kilometers")
abline(v=ger.bw.gl,col="red")