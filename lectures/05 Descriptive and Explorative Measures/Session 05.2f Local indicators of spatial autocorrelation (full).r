################################################################################
#                                                                              #                         
# IDEM 156 - Spatial Demography Course 2018                                    #                          
# Calculating and Mapping Local Indicators of Spatial Autocorrelation          #                                                                                                                                 #
# Sebastian Kl?sener, MPIDR                                                    #                         
#                                                                              #                         
################################################################################

# Erase all objects in workspace
rm(list=ls(all=TRUE))

# Load libraries
library(spdep)
library(RColorBrewer)
library(maptools)
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

# Create neighborhood weight matrices
nb.FOQ <- poly2nb(shape.shp, queen=TRUE)
summary(nb.FOQ)

# Code to connect R?gen
#edit.nb(nb.FOQ, polys=shape.shp)
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
# 2) Create draft maps of some variables                                       #
################################################################################

# Check, which variables we have available
names(shape.shp)

# Shortly produce draft overview maps
var.data <- data.frame(SNM05, TFR05, AVAGE05, UER05)
var.names <-c("Share Nonmarital Births 2005", "TFR 2005", 
 "Average Age at Birth", "Unemployment Rate 2005")
m <- length(var.data)
colors <- matrix(nrow=n,ncol=m, data=rep(0,n*m))

par(mfrow=c(2,2))
# This code produces four maps with equal interval categorization
for (i in 1:m) {
    varofint <- var.data[,i]
    lb <-  5
    bins <- seq(min(varofint), max(varofint), length=lb)
    colpal <- brewer.pal(lb-1, "YlOrRd")
    colors[,i] <- colpal[findInterval(varofint, bins, rightmost.closed=T)]
    plot(shape.shp, col=colors[,i])
    title(var.names[[i]], cex=1)
    legend("bottomright",fill=colpal,
    legend=paste(round(bins[-length(bins)],2),"-",round(bins[-1],2)),
     cex=0.8, bg="white")
}


################################################################################            
# 3) Calculate LISA                                                            #
################################################################################

################################################################################
################################################################################
# In this section we set definitions for the rest of the code

# Which variables are available in the shape.file attribute table
names(shape.shp)

# Define variable of interest
varofint <- AVAGE05

# Define name of the variable of interest 
varofint.name <- c("Average Age at Birth 2005")

# Define neighborhood matrix (type nb) 
# (choice in this case nb.FOQ.cor, d50km, nb.5NN)
nb <- nb.FOQ.cor

# Define weight style (W=row-standardized)
ws <- c("W")

# Define significance level for the cluster maps 
# (0.0001,0.001,0.01,0.05)
significance <- 0.05

# p-adjustment method (can be "none", "bonferroni", "holm",
# "hochberg","hommel","fdr")
p.ad.meth <- c("none")

# Should the cluster map show only regions belonging to significent clusters, 
# or all regions
plot.only.significant <- TRUE


################################################################################
################################################################################

# Transfer nb-object in listwise object
listw <- nb2listw(nb, style=ws)

# Create lagged values of variable of interest
varlag <- lag.listw(listw, varofint)

# Map of Variable and Lagged Variable
par(mfrow=c(1,2))                                      
var.data <- data.frame(varofint, varlag)
var.names <-c(varofint.name, paste("Spatially Lagged\n", varofint.name))
m <- length(var.data)
colors <- matrix(nrow=n,ncol=m, data=rep(0,n*m))
for (i in 1:m) {
    varofint1 <- var.data[,i]
    mean <- mean(varofint1)
    sd <- sd(varofint1)
    min <- min(varofint1)          
    max <- max(varofint1)
    bins <- c(min,mean-sd, mean, mean+sd, max)
    lb <- c(length(bins)-1)
    colpal <- brewer.pal(length(bins-1), "PiYG")
    colors <- colpal[findInterval(varofint1, bins, rightmost.closed=T)]
    plot(shape.shp, col=colors)
    title(var.names[[i]], cex=1)
    legend("bottomright",fill=colpal,legend=paste(round(bins[-length(bins)],2),
    "-", round(bins[-1],2)),cex=0.8)
}

# Calculate Lisa Test
lisa.FOQ <- localmoran(varofint,listw, alternative="two.sided",
                       p.adjust.method=p.ad.meth)


# Get significance level
vec <- c(1:n)
vec <- ifelse(lisa.FOQ[,5] < significance, 1,0)

    
################################################################################            
# 4) Plot LISA-maps                                                            #
################################################################################

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
par(mfrow=c(1,1),mar=c(5, 4, 4, 2) + 0.1)
if (plot.only.significant==TRUE) 
   {plot(shape.shp, col=colors2, border="grey25",lwd=0.7)} else 
   {plot(shape.shp, col=colors1, border="grey25",lwd=0.7)}
   legend("bottomright",fill=c("brown2","royalblue3","lightblue",
                               "pink","white"),
          legend=c("High-High","Low-Low","Low-High","High-Low"),
          border = "grey25", cex=0.8, bg="white", 
          box.col="white")
    title(paste("Significant Clusters\n",varofint.name))

# Significance map
sig.data <- data.frame(lisa.FOQ[,5])
m <- length(sig.data)
colors <- matrix(nrow=n,ncol=m, data=rep(0,n*m))
for (i in 1:m) {
    signific <- sig.data[,i]
    lb <- 5
    bins <- c(0,0.0001,0.001,0.01,0.05,1)
    colpal <- c(rev(brewer.pal(lb+1, "YlGn")[-c(1:2)]),"white")
    colors[,i] <- colpal[findInterval(signific, bins, rightmost.closed=T)]
    colors[vec==0,i] <- c("white")
    plot(shape.shp, col=colors[,i], border="grey25", lwd=0.7)
    title("Significance Level")    
    colpalad <- colpal
    colpalad[c(which(bins==significance):length(colpalad))] <- c("white")
    binsp <- c("0","0.0001","0.001","0.01","0.05",1)
    legend("bottomright",fill=colpalad,
           legend=paste(binsp[-length(bins)],"-",binsp[-1]), border="grey25",
           cex=0.8, bg="white", box.col="white")
}

################################################################################
# 4) Plot results by number of neighbors                                       #           
################################################################################

# Get information on the number of neighbors of each region i
nnei <- card(nb)

# Plot - Local Moran's I
par(mfrow=c(1,1))
plot(nnei, lisa.FOQ[,1], main="Local Moran's I", xlab="Number of Neighbors", 
 ylab="Local Moran's I")
lines(sort(unique(nnei)),aggregate(lisa.FOQ[,1], list(nnei), mean)[,2],
 col="red")

# Plot - Expected Local Moran's I by number of neighbors
plot(nnei, lisa.FOQ[,2])
plot(nnei, lisa.FOQ[,2], main="Expected Local Moran's I", 
 xlab="Number of Neighbors", ylab="Expected Local Moran's I")
lines(sort(unique(nnei)),aggregate(lisa.FOQ[,2], list(nnei), mean)[,2],
 col="red",xlab="Number of Neighbors",ylab="Expected Local Moran's I")

# Plot - Variance of Local Moran's I by number of neighbors
plot(nnei, lisa.FOQ[,3], main="Variance of Local Moran's I", 
 xlab="Number of Neighbors", ylab="Variance of Local Moran's I")
lines(sort(unique(nnei)),aggregate(lisa.FOQ[,3], list(nnei), mean)[,2],
 col="red",xlab="Number of Neighbors",ylab="Variance of Local Moran's I")

# Plot - Standard deviate of Local Moran's I by number of neighbors
plot(nnei, lisa.FOQ[,4], main="Standard deviate of Local Moran's I", 
 xlab="Number of Neighbors", ylab="Standard Deviate")
lines(sort(unique(nnei)),aggregate(lisa.FOQ[,4], list(nnei), mean)[,2],
 col="red",xlab="Number of Neighbors",
 ylab="Standard deviate of Local Moran's I")

# Plot - P-value of Local Moran'I by number of neighbors
plot(nnei, lisa.FOQ[,5], xlab="Number of Neighbors", 
 ylab="P-values of Local Moran's I")
lines(sort(unique(nnei)),aggregate(lisa.FOQ[,5], list(nnei), mean)[,2],
 col="red")

# Map showing Local Moran'I and Variance of Local Moran's I
par(mfrow=c(1,2))                                      
var.data <- data.frame(lisa.FOQ[,1], lisa.FOQ[,3])
var.names <-c("Local Moran's I","Variance of Local Moran's I")
m <- length(var.data)
#colors <- matrix(nrow=n,ncol=m, data=rep(0,n*m))
for (i in 1:m) {
    varofint1 <- var.data[,i]
    mean <- mean(varofint1)
    sd <- sd(varofint1)
    min <- min(varofint1)          
    max <- max(varofint1)
    bins <- c(min,mean-sd*0.5, mean, mean+sd*0.5, max)
    lb <- c(length(bins)-1)
    colpal <- brewer.pal(length(bins-1), "Reds")
    colors3 <- colpal[findInterval(varofint1, bins, rightmost.closed=T)]
    plot(shape.shp, col=colors3)
    title(var.names[[i]], cex=1)
    legend("bottomright",fill=colpal,legend=paste(round(bins[-length(bins)],2),
    "-", round(bins[-1],2)),cex=0.8)
}


################################################################################            
# 5) LISA-plots                                                                #
################################################################################

# Get information on the number of neighbors of each region i
nnei <- card(nb)

# Plot with all observations highlighted according to Cluster-membership
par(mfrow=c(1,1),mar=rep(4,4))
moran.plot(varofint, listw, col="grey50", bg=colors1, labels=FALSE, pch=21,
 xlab=varofint.name, ylab=paste("Spatially Lagged",varofint.name))
title("Cluster membership")
legend("topleft", pt.bg="lightblue", col="grey50", pch=21, legend="Low-High")
legend("topright", pt.bg="brown2", col="grey50", pch=21, legend="High-High")
legend("bottomleft", pt.bg="royalblue3", col="grey50", pch=21, 
       legend="Low-Low")
legend("bottomright", pt.bg="pink", col="grey50", pch=21, legend="High-Low")

# Lisa-Plot with all non-significant observations in white color
moran.plot(varofint,listw, col="grey50", bg=colors2, labels=FALSE, pch=21,
 xlab=varofint.name, ylab=paste("Spatially Lagged",varofint.name))
title("LISA plot - Cluster membership")
legend("topleft", pt.bg="lightblue", col="grey50", pch=21, legend="Low-High")
legend("topright", pt.bg="brown2", col="grey50", pch=21, legend="High-High")
legend("bottomleft", pt.bg="royalblue3", col="grey50", pch=21, legend="Low-Low")
legend("bottomright", pt.bg="pink", col="grey50", pch=21, legend="High-Low")
                        
# Lisa-Plot with point-size varying by number of neighbors
moran.plot(varofint,listw, col="grey50", bg=colors2, labels=FALSE, pch=21, 
           cex=nnei/2, xlab=varofint.name,
           ylab=paste("Spatially Lagged",varofint.name))
title("LISA plot - Cluster membership\n point size varies by number of neighbors j")
legend("topleft", pt.bg="lightblue", col="grey50", pch=21, legend="Low-High")
legend("topright", pt.bg="brown2", col="grey50", pch=21, legend="High-High")
legend("bottomleft", pt.bg="royalblue3", col="grey50", pch=21, 
       legend="Low-Low")
legend("bottomright", pt.bg="pink", col="grey50", pch=21, legend="High-Low")

library(ggplot2)


# Plot with colors according to the significance level, point size varying by
# number of neighbors 
moran.plot(varofint,listw, col="grey50",bg=colors[,1], labels=FALSE, pch=21, 
           cex=nnei/2, xlab=varofint.name, 
           ylab=paste("Spatially Lagged",varofint.name))
title("LISA plot - Significance\n point size varies by number of neighbors j")
legend("topleft", legend="Low-High")
legend("topright", legend="High-High")
legend("bottomleft", legend="Low-Low")
legend("bottomright", legend="High-Low")