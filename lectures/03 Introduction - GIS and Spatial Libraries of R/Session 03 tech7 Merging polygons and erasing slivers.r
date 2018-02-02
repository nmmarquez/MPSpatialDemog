################################################################################
#                                                                              #                         
# IDEM 156 - Spatial Demography Course 2018                                    #                          
# Merge regions with unionSpatialPolygons-function and erase slivers           #
# Sebastian Kl?sener, MPIDR                                                    #                         
#                                                                              #                         
################################################################################

# Erase all objects in workspace
rm(list=ls(all=TRUE))

# Load libraries
library(spdep)
library(maptools)
library(rgdal)
library(rgeos)


################################################################################
# 1) Import and prepare data                                                   #
################################################################################ 

# Set drive
main <- c("N")
# Set path to session folders
path <- ":/IMPRSD/IDEM 156 Spatial"
# Define session folder
path2 <- "/02 Sessions/03 Introduction - GIS and Spatial Libraries of R"
# Set working directory
setwd(paste(main,path,path2,sep=""))

# Open shapefile and attach projection information to it.
shape.shp <- readShapePoly('2004_06', IDvar="KREIS_KENN")
shapeogr.shp <- readOGR(".", "2004_06")
proj4string(shape.shp) <- proj4string(shapeogr.shp)

# Plot shapefile
par(mfrow=c(1,1))
plot(shape.shp,border="grey25")             
names(shape.shp)

# Derive information on internal id and centroids. In this case the ids have 
# been loaded into R as factors. As these ids are in fact all numbers, we turn
# them with the paste and as.numeric commands into numbers, which are easier 
# to handle, if we want to match data to the shapefile.
head(shape.shp$KREIS_KENN)
is.vector(shape.shp$KREIS_KENN)
is.factor(shape.shp$KREIS_KENN)
id <- as.numeric(paste(shape.shp$KREIS_KENN))/1000


################################################################################
# 2) Perform union-operation to merge regions of shapefile                     #
################################################################################ 

# Define merge vector (in our case first two numbers of ID)
idunion <- floor(id/1000)

# Alternative merge vector, in which we just join two regions 
# (with ids 1053 and 2000)
idunion2 <- id
join <- which(id==1053)
idunion2[join] <- 2000

# Perform union procedure to combine polygons based on the information we 
# prepared in our union-vector.
shape.merged <- unionSpatialPolygons(shape.shp, idunion,threshold=1e-03)

# The merged polygon does not contain a data attribute table.
shape.merged@data

# But the ids from the merge vector are available through the names function.
ids <- names(shape.merged)

# Create data frame and name row-names of data frame according to ID.
data <- data.frame(ids)
row.names(data) <- ids

# Connect the data frame to the shapefile
shape.merged <- SpatialPolygonsDataFrame(shape.merged, data)


################################################################################
# 3) Erase slivers from shapefile                                              #
################################################################################ 

# Plot shows that some parts of the internal district boundaries within the
# German states remained, when we performed the union-procedure
par(mfrow=c(1,1),mar=rep(0,4))
plot(shape.merged)

# Derive n
n <- length(shape.merged)


# Define Resultlist in which we store the polygons
resultlist <- list()


# Function to decompose region-wise the regions into their polygons 
# (might be one or several polygons (e.g. mainland, islands, slivers).
# Detect polygons which are not closed (ringDir==-1) and erase these.
# Then store remaining polygons region-wise in list.

for (i in 1:n) {
  # Decomposition of region into its polygons
  pol <- list()
  pol <- shape.merged[i,]@polygons
  lenpol <- length(pol)
  pls <- pol[[1]]
  
  # Extract ringDir-information polygon-wise
  npol <- length(pls@plotOrder)
  nlist <- list()
  for (j in 1:npol) {
    nlist[[j]] <- pls@Polygons[[j]]@ringDir
  }
  # Check whether ringDir is -1 for one of the polygons, which
  # implies the polygon is not closed and probably sliver
  number <- unlist(nlist)
  num <- which(nlist==-1)
  
  # If there are unclosed polygons, they are erased here
  if (length(num)!=0) {red <- pls@Polygons[-num]
                       polys1 <- Polygons(red,ids[i])}
  
  # If not, we just keep the existing polygons
  if (length(num)==0) {polys1 <- Polygons(pls@Polygons,ids[i])}
  
  # Store polygons region-wise in list
  resultlist[[i]] <- polys1
}

# Create Spatial Polygons from resultlist.
spg <- SpatialPolygons(resultlist)
plot(spg)

# We have to match the data frame again
spg@data

# Again we have to take care that the row-names are named according to the ids.
data <- data.frame(ids)
row.names(data) <- ids
shape.clean <- SpatialPolygonsDataFrame(spg, data)

# Check, whether IDs are linked correctly.
plot(spg)
text(coordinates(shape.clean),ids)

# Check, whether file you want to export has projection information included.
proj4string(shape.clean)

# As we did not reproject our data, we can take this information from the 
# shapefile object that we initially loaded.
proj4string(shape.clean) <- proj4string(shape.shp)

# Save cleaned shapefile
writeOGR(shape.clean,".",'shape_cleaned',
         driver="ESRI Shapefile",overwrite_layer=T)