################################################################################
#                                                                              #                         
# IDEM 156 - Spatial Demography Course 2018                                    #                          
# Advanced ways to construct spatial weight matrices                           #
# Sebastian Klüsener, MPIDR                                                    #                         
#                                                                              #                         
################################################################################

# Erase all objects in workspace
rm(list=ls(all=TRUE))

# Load libraries
library(spdep)
library(maptools)
library(rgdal)
library(mapproj)
library(raster)
library(maps)

# Please insert path to your working directory
# Set drive
main <- c("C")
# Set path to session folders
path <- ":/ownCloud/01 MPI/120 Spatial Demography 2018"
# Define session folder
path2 <- "/02 Sessions/11 Spatial Simulations"
# Set working directory
setwd(paste(main,path,path2,sep=""))


################################################################################
# 1) Import and prepare data                                                   #
################################################################################ 

# Download administrative boundaries data for Russia from GADM.org
# Level defines the level of the administrative hierarchy. Download specifies
# whether also a local copy should be stored in your working directory
# Level 1, Country is chosen via ISO3 code
shape.ll <- getData("GADM", country="RUS", level=1,download=TRUE)
head(shape.ll@data)
# Reprojected version more suitable for mapping Russia
crs1 <- c("+proj=laea +lat_0=90 +lon_0=90 +x_0=0 +y_0=0 +datum=WGS84")
crs2 <- c("+units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
crs_rus <- paste(crs1,crs2)
shape.proj <- spTransform(shape.ll,CRS=crs_rus)
plot(shape.proj)

# Load shapefile with all countries of the world
# Could also be obtained from GADM, but the download would take some
# time
world.shp <- readOGR(".", "country")
world.proj <- spTransform(world.shp,CRS=crs_rus)

rus <- which(world.proj$GMI_CNTRY=="RUS")
world.proj_bg <- world.proj[,-rus]  

tiff(file="Russia.tif",width = 3200, height = 2400,
     compression="lzw", res=300)
     plot(shape.proj,col="white",bg="grey95",border="grey25",lwd=0.5)
     plot(world.proj_bg,col="grey75",border="grey25",lwd=0.5,add=T)  
     plot(shape.proj,col="snow1",border="grey25",lwd=0.5,add=T)
     llgridlines(world.proj,lty = 1,lwd=0.5,norths=c(seq(0,90,5)),
                easts=c(seq(0,180,10),seq(-170,-10,10)),plotLabels=F)
dev.off()

# Define ID
id <- shape.shp$ID_1

# Centroids of the reprojected Shapefile with longlat-projection (will be used
# for distance calculatons)
coords.ll <- coordinates(shape.ll)

################################################################################
# 2) Combining two weight matrices                                             #
################################################################################ 

# Distance of 500 km
d500km <- dnearneigh(coords.ll, 0, 500, longlat=T)

# 10 regions do not have a neighbor in a distance 500km-file
summary(d500km)

# In this example you believe that spatial dependency of the process of interest
# is not extending beyond a distance that is beyond 500kms. Thus, you are 
# hesitant to extent the distance value beyond 500km. But you would like to have
# at least for each region one neighbor.

# To achieve this you create a weight file which just contains the first nearest 
# neighbors. 
l.1NN <- knearneigh(coords.ll, k=1, longlat=T)
nb.1NN <- knn2nb(l.1NN)

# For those regions that have a neighbor within 1000km distance, you know that their 
# nearest neighbor is already included in the d1000km weight matrix. Thus running a 
# a union.nb command the two matrices you just include the nearest neighbor of those that
# do not have a nearest neighbor within 1000 km.
united.nb <- union.nb(d500km,nb.1NN)

summary(d500km)
summary(united.nb)

# Other commands
# To obtain links two weight matrices have in common. With this we are back to
# d500km 
intersected.nb <- intersect.nb(d500km,united.nb)
summary(intersected.nb)

# To obtain differences between two weight matrices:
difference.nb <- setdiff.nb(d500km,united.nb)
summary(difference.nb)

# Now you could transfer this neighbor file in a weight matrix.


################################################################################
# 3) Obtaining an explanatory variable which contains combined weighted infos  #
#    from region i and region j.                                               #
################################################################################ 

# Load shapefile for Germany
shape.shp <- readShapePoly('German_Districts_0406')
shapeogr.shp <- readOGR(".", "German_Districts_0406")
proj4string(shape.shp) <- proj4string(shapeogr.shp)
# Reprojected version of the shapefile which has a longlat projection
shape.ll <- spTransform(shape.shp, CRS("+proj=longlat +datum=WGS84"))

# E.g. not only population density in region i as such, but also in surrounding 
# regions j might be important for your outcome variable, but you do not want 
# to add both population density and lagged population density to your model
# due to multicollinearity concerns. Thus you decide to create a variable, 
# which provides information on Population Density both in the area i and 
# neighboring regions j

# Neighbor-file in which Rügen is connected
nb.FOQ.cor <- read.gal('FOQ_German_Districts_corrected.nb')

# Here we have to use the "W"-specification for the matrix to get the average
# value of the region
nb.FOQ.lw.W <- nb2listw(nb.FOQ.cor, style="W",zero.policy=T)

# Population Density
PD <- shape.shp$PD05

# Spatially lagged values of population density
lag.PD <- lag.listw(nb2listw(nb.FOQ.cor),PD)

# Weight you want to give to population density in region i
wi <- 0.6

# Weight you want to give to population density in region j
wj <- 0.4

# Weighted exploratory variable considering information from i and j
PDij.weight <- (wi*PD)+(wj*lag.PD)


################################################################################
# 4) Connectivity matrices                                                     #
################################################################################ 

shape1.shp <- readShapePoly('simple')
shapeogr1.shp <- readOGR(".", "simple")
proj4string(shape1.shp) <- proj4string(shapeogr1.shp)
# Reprojected version of the shapefile which has a longlat projection
shape1.ll <- spTransform(shape1.shp, CRS("+proj=longlat +datum=WGS84"))

n <- length(shape1.ll)

# Centroids of the reprojected Shapefile with LotLang-projection (will be used
# for distance calculatons)
coords1.ll <- coordinates(shape1.ll)

# NB-object in which all are neighbors
nb.all <- dnearneigh(coords1.ll, 0, 100000,longlat=T) 
nb.all

# Turn nb-object in listw-object.
lw.all <- nb2listw(nb.all,style="W")

# Turn listw-object in weight matrix
wij <- listw2mat(lw.all)

# Weighted connectivity information
connectwij <- as.matrix(read.table("Connectivity.csv",sep=",",head=F))

# Connect connectivities and turn in listwise object
wij <- connectwij
lw.con <- mat2listw(wij)

# Assign connectivity weihts to lw.all
lw.allcon <- lw.all

for (i in 1:n) {
     lw.allcon$weights[[i]] <- lw.con$weights[[i]]
}

# Compare
lw.all$weights
lw.allcon$weights
