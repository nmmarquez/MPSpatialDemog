################################################################################
#                                                                              #                         
# IDEM 156 - Spatial Demography Course 2018                                    #                          
# Introduction to the package spdep - Neighborhood weight matrices             #                                                                                                                                 #
# Sebastian Klüsener, MPIDR                                                    #                         
#                                                                              #                         
################################################################################


# Erase all objects in workspace
rm(list=ls(all=TRUE))

# Load libraries
library(spdep)
library(rgdal)
library(maptools)
library(RColorBrewer)

# Set your working directory
# If you need to change the working directory, you can use the following code:
# Set drive
main <- c("C")
# Set path to session folders
path <- ":/ownCloud/01 MPI/120 Spatial Demography 2018"
# Define session folder
path2 <- "/02 Sessions/04 Basic Principles and Challenges of Spatial Analysis"
# Set working directory
setwd(paste(main,path,path2,sep=""))


################################################################################
# 1)     Example: Columbus, Ohio                                               #
#        Contiguity-based weight matrices                                      #
#        Rook vs. Queen                                                        #
# 1.1)   Import data and preparation                                           #
################################################################################

# Import shapefile
columbus <- 
readShapePoly(system.file("etc/shapes/columbus.shp", package="spdep")[1])
plot(columbus)

# Derive coordinates of centroids from shapefile
coords <- coordinates(columbus)

# Derive total number of units
n <- length(columbus)

# Derive id, according to which the attributes are sorted
id <-columbus$POLYID

# Plot city quarters of Columbus, Ohio and polygon centroids
plot(columbus , border="darkgrey")
points(coords)

# Plot city quarters of Columbus, Ohio and order of file
n <- length(columbus)
plot(columbus , border="darkgrey")
text(coords, paste(id), cex=0.8, col="red")


################################################################################
# 1.2)  Create neighborhood weight matrix                                      #
################################################################################

# Create neighborhood weight matrix 
# (First order queen definition, first order rook definition)
# This function is calles "poly2nb" which is the abbreviation for "polygon 'to' 
# neighborlist"
nb.FOQ <- poly2nb(columbus, queen=TRUE, row.names=columbus$POLYID)
nb.FOR <- poly2nb(columbus , queen=FALSE, row.names=columbus$POLYID)
nb.FOQ
nb.FOR

# Plot first order queen weight matrix
plot(columbus, border="darkgrey")
plot.nb(nb.FOQ, coords, add=TRUE)
text(coords-0.07, paste(id), cex=0.8, col="red")

# Look at neighborhood file in Detail
nb.FOQ[1:n]

# Export neighborhood matrix (might be helpful, if it takes long to compute it)
write.nb.gal(nb.FOQ,'cor.col.nb.GAL', ind=columbus$POLYID)

# Import neighborhood matrix)
nb.FOQ1 <- read.gal('cor.col.nb.GAL')

# Test to check, whether exported and imported file are identical
difference.nb <- diffnb(nb.FOQ, nb.FOQ1)
difference.nb

# Contrast differences between rook and queen neighborhood matrices
difference.nb <- diffnb(nb.FOQ, nb.FOR)
plot(columbus, border="darkgrey")
plot.nb(difference.nb, coords, add=TRUE, col="darkred")
text(coords-0.07, paste(id), cex=0.8, col="red")


# Plot map and highlight region i and its neighbors j
reg_i <- c(15)
for (i in 1:length(reg_i)) {
    region <- reg_i[i]
    neighbors <- as.vector(nb.FOQ[[reg_i]])
    color <- rep("white",n)
    color[region] <- "green"
    color[neighbors] <-"lightgreen"
    plot(columbus, border="darkgrey", col=color)
    plot.nb(nb.FOQ, coords, add=TRUE)
    text(coords-0.07, paste(id), cex=0.8, col="black")
}

################################################################################
# 1.3)   Change neighborhood weight matrix                                     #
#        Delete connection between two regions with a coded procedure          #
################################################################################

# a) Definitions for automized code to delete neighborhood relations

# Define nb-File
nb <- nb.FOQ
# Define shapefile
shape <- columbus
# Define regions to be deleted (enter regions in pairs, whose neighborhood  
# link should be deleted)
d.b.r.m <- c(35, 40, 31, 39)

# Test, whether neighborhood weight file is symmetric:
# if a is neighbor of b, b is also a neighbor of a 
is.symmetric.nb(nb)

# Automatized code to delete neighborhood relations
coords <- coordinates(shape)
n <- length(shape)
len.d <- length(d.b.r.m)/2
d.b.r <- matrix(ncol=2,nrow=len.d, data=d.b.r.m, byrow=TRUE)
coords <- coordinates(shape)
for (i in 1:len.d) {
    neigh1 <- nb[[d.b.r[i,1]]]  
    neigh2 <- nb[[d.b.r[i,2]]]  
    er1 <- neigh1[-(which(neigh1==d.b.r[i,2]))]
    er2 <- neigh2[-(which(neigh2==d.b.r[i,1]))]
    nb[[d.b.r[i,1]]] <- er1
    nb[[d.b.r[i,2]]] <- er2
}

# Test, whether neighborhood weight file is still symmetric
is.symmetric.nb(nb)


# Map results
color <- rep("white",n)
color[d.b.r.m] <- "green"
plot(shape, border="darkgrey", col=color) 
plot.nb(nb, coords, add=TRUE)
text(coords-0.07, paste(id), cex=0.8, col="darkred")

# b) Definitions for automized code to create neighborhood relations
# Define nb-File
nb <- nb
# Define shapefile
shape <- columbus
# Define pairs of regions to be connected 
c.b.r.m <- c(35, 40, 31, 39,2,21)

# Test, whether neighborhood weight file is still symmetric
is.symmetric.nb(nb)

# Automatized code to create neighborhood relations
coords <- coordinates(shape)
n <- length(shape)
len.c <- length(c.b.r.m)/2
c.b.r <- matrix(ncol=2,nrow=len.c, data=c.b.r.m, byrow=TRUE)
for (i in 1:len.c) {
    neigh1 <- nb[[c.b.r[i,1]]]  
    neigh2 <- nb[[c.b.r[i,2]]]  
    cr1 <- unique(sort(c(neigh1,c.b.r[i,2])))
    cr2 <- unique(sort(c(neigh2,c.b.r[i,1])))
    nb[[c.b.r[i,1]]] <- as.integer(cr1)
    nb[[c.b.r[i,2]]] <- as.integer(cr2)
# This code erases the zero in case the region did not have a neighbor before
    if (nb[[c.b.r[i,1]]][1]==0) nb[[c.b.r[i,1]]] <- nb[[c.b.r[i,1]]][-1] 
    if (nb[[c.b.r[i,2]]][1]==0) nb[[c.b.r[i,2]]] <- nb[[c.b.r[i,2]]][-1]
}
# Test, whether neighborhood weight file is still symmetric
is.symmetric.nb(nb)

color <- rep("white",n)
color[c.b.r.m] <- "green"
plot(shape, border="darkgrey", col=color) 
plot.nb(nb, coords, add=TRUE)
text(coords-0.07, paste(id), cex=0.8, col="darkred")


################################################################################
# 1.4 Contiguity-based weight matrices (Second or higher order)                #
################################################################################

# First Order Queen
nb.FOQ <- poly2nb(columbus, queen=TRUE, row.names=columbus$POLYID)

# Second Order Queen with function nblag. This function creates a list with
# n elements depending on the number of orders (in our case 2 - the first
# containing the first order neighbors, the second the second order neighbors.
nb.SOQ <- nblag(nb.FOQ, 2)

# First Order Neighbors (in total 236)
nb.SOQ[[1]]
# Second Order Neighbors (in total 416)
nb.SOQ[[2]]

# Combine first and second order lags into one nb-file containing both first and
# second order links (in total 236+416=652)
nb.FSOQ <- nblag_cumul(nb.SOQ)
nb.FSOQ


################################################################################
# 2) Turn neighborhood weight file into a listwise object or a spatial         #
#    neighbor object (might be neccessary for calculations)                    #         
################################################################################


# Erase all objects in memory
rm(list = ls(all = TRUE))

# Import shapefile
shape.shp <- readShapePoly('simple', IDvar="ID")
shapeogr.shp <- readOGR(".", "simple")
is.projected(shape.shp)
is.projected(shapeogr.shp)
proj4string(shape.shp) <- proj4string(shapeogr.shp)

# Derive matrix with centroids
coords <- coordinates(shape.shp)

# Derive total number of units
n <- length(shape.shp)

# Derive id, according to which the attributes are sorted
id <-shape.shp$ID

# Create a vector with numbers from 1 to total number of regions in the 
# shapefile
sq <- c(1:n)

# Turn numbers into text
sq_text <-paste(c(1:n))


# Creating first order queen weight matrix
nb.FOQ <- poly2nb(shape.shp, queen=TRUE,row.names=id)
summary(nb.FOQ)

# Again the map for illustrations
plot(shape.shp, border="darkgrey") 
test <- text(coords, sq_text, cex=2, col="red") 

# Based on row standardised weights (style="W")
# The next two lines produce an error, as there is one region, which has no 
# neighbors
poly2nb(shape.shp, queen=TRUE,row.names=id)
nb.FOQ.lw.W <- nb2listw(nb.FOQ, style="W")

# Here we use the file, where Rügen has been linked to two neighboring regions
# Row-standardized (W)
nb.FOQ.snap <- poly2nb(shape.shp, queen=TRUE,row.names=id, snap=2000)
nb.FOQ.lw.W <- nb2listw(nb.FOQ.snap, style="W")
nb.FOQ.lw.W$weights
# Just Binary, not row-standarized (style="B")
nb.FOQ.lw.B <- nb2listw(nb.FOQ.snap, style="B")
nb.FOQ.lw.B$weights

# Turn into weight matrices
# Binary weight matrix
wij <- listw2mat(nb.FOQ.lw.W)
cij <- listw2mat(nb.FOQ.lw.B)

# Matrix vs. spatial neighbour object -> saving computation time by just 
# creating a list of all connections
nb.FOQ.lw.sn <- listw2sn(nb.FOQ.lw.W)
nb.FOQ.lw.sn

# Contrast length of neighborhood matrix and spatial neighbor object                   
len.matrix <- n*n                   
len.spatial.neigh.object <-length (nb.FOQ.lw.sn$from)

len.matrix
len.spatial.neigh.object


################################################################################
# 3) Weight matrices based on nearest neighbors and distance                   #
################################################################################

shape.ll <- spTransform(shape.shp, CRS("+proj=longlat +datum=WGS84"))
coords.ll <- coordinates(shape.ll)
spcoords <- SpatialPoints(coords.ll,proj4string=CRS("+proj=longlat +datum=WGS84"))

 
# 2.4.1) k-nearest neighbors
l.5NN <- knearneigh(coords.ll, k=5,longlat=T)
nb.5NN <- knn2nb(l.5NN)

# Such a neighbors list is usually not symmetric!
is.symmetric.nb(nb.5NN)

# 2.4.2) Creating a distance based matrix (distance between centroids)

#Distance of 20 km
d20km <- dnearneigh(coords.ll, 0, 20, longlat=T)

# Distance of 50 km
d50km <- dnearneigh(coords.ll, 0, 50, longlat=T)

# Contrast all neighborhood weight matrices in plots
par(mfrow=c(2,2))
plot(shape.ll, border="darkgrey")
plot.nb(nb.FOQ, coords.ll, add=TRUE)
title("First Order Queen")
plot(shape.ll, border="darkgrey")
plot.nb(nb.5NN, coords.ll, add=TRUE)
title("5 Nearest Neighbors")
plot(shape.ll, border="darkgrey")
plot.nb(d20km, coords.ll, add=TRUE )
title("Distance 20 km")
plot(shape.ll, border="darkgrey")
plot.nb(d50km, coords.ll, add=TRUE)
title("Distance 50 km")


################################################################################
# 4) Weight matrices weighted by distance                                      #
################################################################################

# 2.5.1) Inverse distance
# Create a distance-based list, where all regions i neighbor all regions j
# To be on the save side, take a distance beyond 22000 km.
nb.all <- dnearneigh(coords.ll, 0, 100000,longlat=T) 
nb.all

# nbdists <- Function to calculate for a neighborlist the distances between
# the centroids of regions which are defined as neighbors 
dist <- nbdists(nb.all,coords.ll,longlat=T)
dist

# In order to receive the inverse distance, we take the just created list
# of distances between the centroids of neighboring regions a lapply-
# functions. The lapply performs for each element (in our case the distance 
# between two regions i and j) the same calculation. We take the inverse 
# distance by dividing 1/distance

# Here we do it manually for the first element
head(dist)
# First distance in the list
fdist <- dist[[1]][1]
finvdist <- 1/fdist
finvdist

# Here we apply the lapply function for the first element
invdist <- lapply(dist, function(x) 1/x)
# Check, whether we came to the same result for the first element
head(invdist)

# Turn neighbor object in listwise object
inv.listw <- nb2listw(nb.all, glist=invdist, style="B") 

# In order to check whether lists are symmetrical (they should be!)
col.mat <- listw2mat(inv.listw)
all.equal(col.mat, t(col.mat), check.attributes=FALSE) 


# 2.5.2) Inverse squared distance
# Code very similar to to the inverse distance. We just change the formula in
# the lapply-function

invsqrdist <- lapply(dist, function(x) 1/x^2)
invsqr.listw <- nb2listw(nb.all, glist=invsqrdist, style="B") 
col.mat.sqr <- listw2mat(invsqr.listw) 

# 2.5.3) Double-power distance weight
# Flexible function which combines a finite bandwidth with bell shaped
# taper function
# Definition of finite bandwidth (weight is 0 beyong bandwidth)
bandwidth <- 50                        
# Here we define a quadratic distance function for k
# As an alternative, we could also set k to 3 or 4. 
k <- 2

# These two lines of code check which of the distances between regions are above
# the threshold distance. They are stored in the vector "er"
col.distd <- unlist(listw2mat(nb2listw(nb.all, glist=dist, style="B") ))
er <- which(col.distd>bandwidth)

# Code very similar to to the inverse distance. We just change the formula in
# the lapply-function
dpdwdist <- lapply(dist, function(x) (1-(x/bandwidth)^2)^2)
dpdw.listw <- nb2listw(nb.all, glist=dpdwdist, style="B") 

# Here we transfer the listwise document to a matrix
col.mat.dpdw <- listw2mat(dpdw.listw) 
# All weights for which the distance between regions is above to zero. For this
# we use the vector "er" which we have created above
col.mat.dpdw[er] <- 0
# Retransform matrix to listwise element.
dpdw1.listw <- mat2listw(col.mat.dpdw)

# Plot of distance weigthed matrices (for region 3 and its neighbors)
#png(file="Distances.png",width = 2400, height = 2400, res=300)
par(mfrow=c(2,2))
region <- 2
maxd <- max(unlist(dist))
plot(0:50,rep(1,51),type="l",xlim=c(1,maxd),ylim=c(0,1.5),
 xlab="Distance in km",ylab="Weight",main="Bandwidth=50 km")
lines(50:maxd,rep(0,maxd-49))
lwd50 <-nb2listw(d50km,style="B")
col.mat.pl <- listw2mat(lwd50)
points(dist[[region]],col.mat.pl[region,-region],col="red")

region <- 2
plot(0:maxd,1/0:maxd,type="l",, xlab="Distance in km",ylab="Weight",
 ylim=c(0,0.2), main="Inverse Distance")
points(dist[[region]],col.mat[region,-region],col="red")

region <- 2
plot(0:maxd,1/(0:maxd)^2,type="l",, xlab="Distance",ylab="Weight",ylim=c(0,0.1),
 main="Inverse Squared Distance")
points(dist[[region]],col.mat.sqr[region,-region],col="red")

region <- 2
plot(0:maxd,c((1-(0:bandwidth/bandwidth)^2)^2,rep(0,maxd-bandwidth)),type="l",
 xlab="Distance",ylab="Weight",ylim=c(0,1),xlim=c(0,maxd),
 main="Double-Power Distance Weight\nbw=50, k=2")
points(dist[[region]],col.mat.dpdw[region,-region],col="red")
#dev.off()


################################################################################
# 5)     Create map showing neighbor connections                               #
#          Sometimes shapefiles are corrupt, which might result in wrong       #
#          neighborhood weightfiles. These maps might help to check the file.  #
################################################################################

# Derive number of neighbors per unit
nom.of.neigh <- c(rep(0,n))
nom.of.neigh1 <- nom.of.neigh
for (i in 1:n) {
  nom.of.neigh[i] <- as.integer(length(nb.FOQ[[i]]))
} 
# for-loop, which is correct, but creates warnings:
for (i in 1:n) {
    if (as.integer(nb.FOQ[[i]])==0) nom.of.neigh1[i] <- 0
    else nom.of.neigh1[i] <- 1
}
# Same function without warnings:
for (i in 1:n) {
    if (mean(as.integer(nb.FOQ[[i]])==0)) nom.of.neigh1[i] <- 0
    else nom.of.neigh1[i] <- 1
} 
non <- nom.of.neigh*nom.of.neigh1

# Create color scheme for map, in which regions with no neighbor are highlighted 
# in red
nonc <- nom.of.neigh*nom.of.neigh1
for (i in 1:n) {
    if (nonc[i]<1) nonc[i] <- "red"
    else nonc[i] <- "white" 
} 

# This is an example where the nonc-information is matched to the shape-file. 
# However, this is not neccessary, as these can be kept seperate as long as the
# order of the shapefile-attribute table and the order of the derived data
# is not changed.
shape.shp.nb <- spCbind(shape.shp, nonc)

par(mfrow=c(2,2))
# First plot - First order queen neighborhood weight file
plot(shape.ll, border="darkgrey")
title("First order queen neighborhood weight matrix")
plot.nb(nb.FOQ, coords.ll, add= TRUE, col="red")

# Second plot - Map, in which regions without neighbor are highlighted
plot(shape.shp.nb, col=nonc)
title("Regional units without neighbor")
legend("topleft",fill=c("white","red"), 
       legend=paste(c("neighors","no neighbors")), cex=1, bg="white")

# Create map, in which regions are colored according to number of neighbors
non10 <- non
non10[non10>8] <-8
colpal   = brewer.pal(9, "YlOrRd")
color <- rep(0,n)

for (i in 1:n) {
    color[i] <- colpal[non10[i]+1] 
} 

# Map indicating the number of neighbors
plot(shape.shp.nb, col=color)
title("Neighbors per regions")
legend("topleft",fill=colpal, legend=c(paste(c(0:7)),">7"),
 cex=1, bg="white")
