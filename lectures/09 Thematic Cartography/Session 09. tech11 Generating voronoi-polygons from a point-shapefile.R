################################################################################
#                                                                              #                          
# IDEM 156 - Spatial Demography Course 2018 - Technical File 13                #                                                    
# Creation of a Voronoi-polygon shapefile from a point-shapefile               #
# Tgese are also referred to as Thiessen polygons
# Voronoi Function by Carson Farmer (slightly adapted by S. Klüsener)          #
# Source: http://carsonfarmer.com/2009/09/voronoi-polygons-with-r/             #
#                                                                              #                          
################################################################################

# Erase all objects in workspace
rm(list=ls(all=TRUE))

# Load libraries
library(spdep)
library(deldir)
library(maptools)
library(rgeos)
library(rgdal)

# Set your working directory
# If you need to change the working directory, you can use the following code:
# Set drive
main <- c("N")
# Set path to session folders
path <- ":/IMPRSD/IDEM 156 Spatial"
# Define session folder
path2 <- "/02 Sessions/09 Thematic Cartography"
# Set working directory
setwd(paste(main,path,path2,sep=""))


################################################################################
#                                                                              #
# 1) Preparation                                                               #
#                                                                              #
################################################################################ 


# Code for the function creating Voronoi-Polygons out of point-shapefiles
# Voronoi Function by Carson Farmer (slightly adapted by S. Klüsener)                                            
# Source: http://carsonfarmer.com/2009/09/voronoi-polygons-with-r/ 
voronoipolygons <- function(x) {
  require(deldir)
  require(sp)
  if (.hasSlot(x, 'coords')) {
    crds <- x@coords  
  } else crds <- x
  z <- deldir(crds[,1], crds[,2])
  w <- tile.list(z)
  polys <- vector(mode='list', length=length(w))
  for (i in seq(along=polys)) {
    pcrds <- cbind(w[[i]]$x, w[[i]]$y)
    pcrds <- rbind(pcrds, pcrds[1,])
    polys[[i]] <- Polygons(list(Polygon(pcrds)), ID=as.character(i))
  }
  SP <- SpatialPolygons(polys)
  voronoi <- SpatialPolygonsDataFrame(SP, 
             data=data.frame(x@data, 
             row.names=sapply(slot(SP, 'polygons'), 
             function(x) slot(x, 'ID'))))
}

# Countries
world <- readShapePoly("CNTRY92",IDvar="FIPS_CODE")
worldogr <- readOGR(".", "CNTRY92")
proj4string(world) <- proj4string(worldogr)

# Load the point shapefile 'CITIES':
cities <- readShapePoints("CITIES")
citiesogr <- readOGR(".", "CITIES")
proj4string(cities) <- proj4string(citiesogr)


################################################################################
#                                                                              #
# 2) Creation of Voronoi-Shapefile                                             #
#                                                                              #
################################################################################ 


# Create voronoifile for all cities of the world
voronoifile <- voronoipolygons(cities)
plot(voronoifile)

# Just for cities in Ghana, Togo, Benin and Nigeria
cities_reduced <- cities[which(cities$COUNTRY=="Ghana"|
                                 cities$COUNTRY=="Nigeria"|
                                 cities$COUNTRY=="Togo"|
                                 cities$COUNTRY=="Benin"),]

voronoifile_AF.shp <- voronoipolygons(cities_reduced)
proj4string(voronoifile_AF.shp) <- proj4string(cities_reduced)



# Example plot
# First line defines the plot window which is focused on the are between 
# Ghana and Nigeria
plot(world[which(world$NAME=="Ghana"|world$NAME=="Nigeria"),],bg="#F0FFFF")
# After the extent of the plor window has been defined, the whole world 
# is plotted.
plot(world,border="grey80",col="white",add=T)
# Plot of the Voronoifile
plot(voronoifile_AF.shp,add=T)
# Alternative - Voronoifile for all cities of the world
#plot(voronoifile,add=T)
# Cities for which voroinoi-polygons are created
points(cities[which(cities$COUNTRY=="Ghana"|
                    cities$COUNTRY=="Nigeria"|
                    cities$COUNTRY=="Togo"|
                    cities$COUNTRY=="Benin"),],col="red",pch=19)


################################################################################
#                                                                              #
# 2) Cut out parts of Voronoi-Shapefile situated in the sea                    #
#                                                                              #
################################################################################ 



# Take out all parts of the voronoifile which are situated in the sea
# First merge all countries of the world in one polygon without any national 
# borders, which then just shows the coastal-line
world_without_borders <- unionSpatialPolygons(world, rep(1,length(world)))
# Use a cookie-cutter function to cut out all
voronoi_land_AF_pr.shp <- gIntersection(voronoifile_AF.shp, 
                                        world_without_borders, byid = TRUE)
# Add SpatialPolygonsDataframe
voronoi_land_AF.shp <- SpatialPolygonsDataFrame(voronoi_land_AF_pr.shp,
                                                voronoifile_AF.shp@data,
                                                match.ID=F)


png(file="Voronoi_map.png",width = 1200, height = 1200, res=300)
# Plot with the improved voronoi-shapefile
# First line defines the plot window which is focused on the are between 
# Ghana and Nigeria
par(mar=rep(0,4))
plot(voronoi_land_AF.shp,border="white",bg="#F0FFFF")
# After the extent of the plor window has been defined, the whole world is 
# plotted
plot(world_without_borders,border="grey80",col="grey90",add=T)
# Plot of the Voronoifile
plot(voronoi_land_AF.shp,border="grey50",add=T)
# Cities for which voroinoi-polygons are created
points(cities[which(cities$COUNTRY=="Ghana"|
                    cities$COUNTRY=="Nigeria"|
                    cities$COUNTRY=="Togo"|
                    cities$COUNTRY=="Benin"),],col="red",pch=19)
dev.off()


################################################################################
#                                                                              #
# 3) Save Voronoi-shapefile                                                    #
#                                                                              #
################################################################################ 

# Export Voronoi-Polygons
writeOGR(voronoi_land_AF.shp,".",'voronoi_polygons',
         driver="ESRI Shapefile",overwrite_layer=T)
