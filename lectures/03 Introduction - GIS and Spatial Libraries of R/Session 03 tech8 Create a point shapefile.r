################################################################################
#                                                                              #                          
# IDEM 156 - Spatial Demography Course 2018 - Technical File 7                 #                          
# Generate a point shapefile                                                   #                          
# Sebastian Kl?sener, MPIDR                                                    #                          
#                                                                              #                          
################################################################################


# Erase all objects in workspace
rm(list=ls(all=TRUE))

# Load libraries
library(spdep)
library(rgdal)
library(maptools)


################################################################################
# 1) Import and prepare data                                                   #
################################################################################ 

# Open file with data and information on longitude and latitude
# E.g., obtained from http://itouchmap.com/latlong.html 
mydata <- read.table("TopDestGermany.csv", sep=",", header=T)
names(mydata)

# Define a data frame or matrix with the coordinates
mycoordinates <- data.frame(mydata$East,mydata$North)


################################################################################
# 2) Create point shapefile                                                    #
################################################################################ 

# Create a spatial point shapefile and attach to it your data as 
# Spatial Object data frame. As long as you derived the coordinates 
# from the same files and have not changed their order, you can be sure,
# that the coordinates and the data will be in the same row-wise order. 
pointshape.shp <- SpatialPointsDataFrame(mycoordinates, data=mydata)

# This shapefile has no projection information attached.
proj4string(pointshape.shp) 
  
# Define projection (in this case you know that the shapefile is based
# on simple information on the longitude and latitude of locations)
proj4string(pointshape.shp) <- CRS("+proj=longlat") 
proj4string(pointshape.shp)

# Save your new shapefile
writeOGR(pointshape.shp, ".",'mypointshapefile', driver="ESRI Shapefile")


################################################################################
# 3) Plot of your point shapefile                                              #
################################################################################ 

# Open shapefile of Germany as background map
shapeogr.shp <- readOGR(".", "2008_w_fert_data")

# This shapefile is based on a UTM-projection
proj4string(shapeogr.shp)
  
# Change projection of that shapefile to projection of point shapefile 
# (longlat)
shape.ll <- spTransform(shapeogr.shp, proj4string(pointshape.shp))

# Plot map
par(mfrow=c(1,1), mar=c(1,1,2.5,1))
# First plot is to define the plot window and to mark the borders of 
# Germany in bold
plot(shape.ll, border="grey75", lwd=1)
# Here we add a title
title("My self-generated Point Shapefile")
# Here we erase all internal borders of Germany by plotting shape.ll in
# white with the line type of the borders being "0"
plot(shape.ll, ,col="grey98",lty=0, add=TRUE)
# Plot of our point shapefile
plot(pointshape.shp, add=TRUE, pch=16, col="firebrick1",cex=0.8)
# Extract point coordinates
coord <- coordinates(pointshape.shp) 
# Add text description to our points
text(coord,c(paste(mydata$Tourist)),pos=c(1,4,3,3,3,1), cex=0.5)


################################################################################
# 4) Export Map                                                                #
################################################################################ 

# Export as png-file
png(file="Pointshapefile.png",width = 1200, height = 1200, res=300)
   par(mfrow=c(1,1), mar=c(1,1,2.5,1))
   # First plot is to define the plot window and to mark the borders of 
   # Germany in bold
   plot(shape.ll,border="grey75", lwd=1)
   # Here we add a title
   title("My self-generated Point Shapefile")
   # Here we erase all internal borders of Germany by plotting shape.ll in
   # white with the line type of the borders being "0"
   plot(shape.ll, col="grey98",lty=0, add=TRUE)
   # Plot of our point shapefile
   plot(pointshape.shp, add=TRUE, pch=16, col="firebrick2",cex=0.8)
   # Extract point coordinates
   coord <- coordinates(pointshape.shp) 
   # Add text description to our points
   text(coord,c(paste(mydata$Tourist)),pos=c(1,4,3,3,3,1), cex=0.5)
dev.off()
