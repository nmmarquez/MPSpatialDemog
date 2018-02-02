################################################################################
#                                                                              #                          
# IDEM 156 - Spatial Demography Course 2018 - Technical File 15                #                          
# Combining Raster Data with Point and Polygon Shapefiles                      #
#                                                                              #
# Sebastian Kl?sener, MPIDR                                                    #                          
#                                                                              #                          
################################################################################

# Erase all objects in workspace
rm(list=ls(all=TRUE))

# Install libraries
#install.packages(c("raster","rgdal","rasterVis","rgl"))

# Open libraries
library(sp)
library(rgdal)
library(rgeos)
library(raster)
library(rasterVis)
library(rgl)
library(maptools)
library(classInt)

# Set your working directory
# If you need to change the working directory, you can use the following code:
# Set drive

################################################################################
#                                                                              #
# 1) For this example we use elevation data from the SRTM 90m resolution       #
#    digitial elevation database:                                              #
#    http://www.cgiar-csi.org/data/srtm-90m-digital-elevation-database-v4-1    #
#    Data at more detailed resolution would, for example, be available here:   #
#    https://www.ngdc.noaa.gov/mgg/global/                                     #  
#                                                                              #
################################################################################ 

# Fetch Dadministrative boundaries data for Nepal from GADM.org
# Level defines the level of the administrative hierarchy. Download specifies
# whether also a local copy should be stored in your working directory
# Level 1, Country is chosen via ISO3 code
shape.shp <- getData("GADM", country="NPL", level=2,download=TRUE)
head(shape.shp@data)

altitude.shp <- getData('alt', country='NPL', mask=F,download=T)

plot(altitude.shp,col=terrain.colors(255))
plot(shape.shp,add=T,border="grey25")

## Interactive 3D plot of the elevation file
plot3D(altitude.shp)

# Write our raster object to file
writeRaster(altitude.shp,"nepal.tif",overwrite=T)

# Read our raster object from file
nepal <- raster("nepal.tif")

# To cookie-cut the elevation data for Nepal out of 
# the elevation file, we first need to turn the administrative 
# shapefile in a raster object, using the raster detail of the 
# elevation file
nepal_mask_raster <- rasterize(shape.shp, altitude.shp)
# Rasterized polygon shapfile
plot(nepal_mask_raster)

# Only keep the points of the elevation raster, which are also in the
# rasterized polygon shapefile
nepal.ele.mask <- mask(altitude.shp, nepal_mask_raster)
plot(nepal.ele.mask,col=terrain.colors(100))
   
# Check elevation profile to choose better color scheme
hist(nepal.ele.mask,main="Elevation Profile of Nepal")
brks <- c(seq(0,2000,100),seq(2200,6000,200),7000,9000)
abline(v=brks,col="red",lty=2)
# Plot with improved color scheme
plot(nepal.ele.mask,col=terrain.colors(length(brks)-1),breaks=brks, legend=F)

# What is the median elevation in Sudan?
summary(nepal.ele.mask)

# Aggregate raster to create a raster object with a lower resolution (larger 
# cells)
nepal.ele.mask.agg <- aggregate(nepal.ele.mask,fact=10,fun=mean)
length(nepal.ele.mask)
length(nepal.ele.mask.agg)
plot(nepal.ele.mask.agg,breaks=brks,col=c("transparent",
                                          terrain.colors(length(brks)-2)),
                                          legend=F)

# Remove original raster file from environment
rm(altitude.shp)


################################################################################
#                                                                              #
# 2) We want to derive for the Nepalese regions information on the elevation   #
#   (mean, median, standard deviation, terrain ruggedness). For this we need   #
#   to perform an overlay operation between the shapefile and the rasterfile.  #    
#                                                                              #
################################################################################ 

# Now we can extract values for all the regions. For this we first create an
# empty list in which we store the documents.
terrain <- list()

# Derive number of regions
l.ob <- length(shape.shp)

#extract terrain
for (i in 1:l.ob) {   
    terrain[[i]] <-  extract(nepal.ele.mask,shape.shp[i,])
}

# In the list of 26 regions the raster points with the information
# are stored in one list. In order to derive summary statistics
# with the lapply command, we first unlist this internal lists.
t <- lapply(terrain,unlist)

# Now we can extract the summary statistics with the lapply command.
# mean, median, min, max, sd
mean.ele <- lapply(t,mean)
med.ele <- lapply(t,median)
min.ele <- lapply(t,min)
max.ele <- lapply(t,max)
sd.ele <- lapply(t,sd)

# Generate dataset with outcomes. This you could, e.g. join, with the shapefile
# by using the spCbind or other commands
ele.dat <- data.frame(shape.shp$ID_2,paste(shape.shp$NAME_2,
                                           shape.shp$HASC_2,sep="_"),
                     unlist(mean.ele), unlist(med.ele),
                     unlist(min.ele), unlist(max.ele),
                     unlist(sd.ele))
colnames(ele.dat) <- c("ID", "Name", "Mean","Median","Minimum",
                       "Maximum","SDev")

# Mean vs median
plot(ele.dat$Median,ele.dat$Mean,xlab="Median Elevation",ylab="Mean Elevation")


# Derive terrain ruggedness
# To calculate terrain ruggedness, we need the exact positions of the raster
# points. Thus, we cannot use extract to extract the values. Instead we have to
# create for each region a seperate master file.
ele_reg <- list()
for (i in 1:l.ob) {
  # Crop out area around location 
  loc.reg <- crop(nepal.ele.mask,shape.shp[i,])
  # Mask pixels outside the buffer region
  ele_reg [[i]] <- mask(loc.reg,shape.shp[i,])
  # Reports us were we currently are
  print(i)
}
# Plot of one region
plot(ele_reg[[2]],col=terrain.colors(length(brks)-1),breaks=brks, legend=F)

# Calculate ruggedness index and st
rug <- list()
sdev <- list()

# Focal moving window (weight): a window of three rows by three columns 
# focused on the central raster point
f <- matrix(1, nrow=3, ncol=3)

# Length of all regions
l.reg <- length(ele_reg)

# Extract ruggedness information
for (i in 1:l.reg) {
  loc <- ele_reg[[i]]
  l.loc <- length(loc)
  # This function is calculating the absolute difference between
  # the cental point in the 
  TRI <- focal(loc, w=f, fun=function(x, ...) sum(abs(x[-5]-x[5]))/8,
               pad=TRUE, padValue=NA)
  # Mean of ruggedness
  rug[[i]] <- mean(TRI[1:l.loc],na.rm=T)
  # Standard deviation information
  sdev[[i]] <- sd(loc[1:l.loc],na.rm=T)
  print(i)
}

# Terrain Ruggedness
ele.dat$TRI <- unlist(rug)
ele.dat$SDev1 <- unlist(sdev)

# Contrasting Terrain Ruggedness with standard deviation
plot(ele.dat$TRI,ele.dat$SDev)


################################################################################
#                                                                              #
# 3) Analogous Procedure for a point dataset                                   #
#                                                                              #
################################################################################ 

# Load Library
library(RJSONIO)

cityname <- c("Kathmandu", "Pokhara Lekhnath","Lalitpur")
countrycode <- c("NPL","NPL","NPL")
mydata <- data.frame(cityname,countrycode)
colnames(mydata) <- c("CityLong","CountryCode")

# Extracting geocodes (latlong-projected)
# Based on a code suggested by Jochem
# http://stackoverflow.com/questions/13905098/how-to-get-the-longitude-and-latitude-coordinates-from-a-city-name-and-country-i
# Cities
nrow <- nrow(mydata)
counter <- 1
mydata$lon[counter] <- 0
mydata$lat[counter] <- 0
while (counter <= nrow){
  CityName <- gsub(' ','%20',mydata$CityLong[counter]) #remove space for URLs
  CountryCode <- mydata$Country[counter]
  url <- paste(
    "http://nominatim.openstreetmap.org/search?city="
    , CityName
    , "&countrycodes="
    , CountryCode
    , "&limit=9&format=json"
    , sep="")
  x <- fromJSON(url)
  if(is.vector(x)){
    mydata$lon[counter] <- x[[1]]$lon
    mydata$lat[counter] <- x[[1]]$lat    
  }
  counter <- counter + 1
}
mydata

# Define a dataframe or matrix with the coordinates
mycoordinates <- data.frame(as.numeric(paste(mydata$lon)),
                            as.numeric(paste(mydata$lat)))

# Create a Spatial Point Shapefile and attach to it your data as 
# Spatial Object dataframe. As long as you derived the coordinates 
# from the same files and have not changed their order, you can be sure,
# that the coordinates and the data will be in the same rowwise order 
pointshape.shp <- SpatialPointsDataFrame(mycoordinates, 
                                         data=mydata)
# This file is similar as the raster file also in longlat projection
proj4string(pointshape.shp) <- proj4string(nepal.ele.mask)

# Check location of cities
plot(nepal.ele.mask,col=terrain.colors(length(brks)-1),breaks=brks, 
     legend=F)
plot(pointshape.shp,add=T,pch=21,bg="red",col="black")

# Here we are just extracting the value for the raster point might:
cities.ele <- list()
l.cit <- length(pointshape.shp)
for (i in 1:l.cit) {   
  cities.ele[[i]] <-  extract(nepal.ele.mask,pointshape.shp[i,])
}
# This procedure is likely to create bias for several reasons:
# - the point-shapefile and the raster file might not perfectly match (also
#   depending on the decimals used to derive location information for the 
#   shapefile)
# - areas of big cities are likely to extend beyong one raster point (also 
#   depends on the detail of the raster)

# A possible solution is to include a buffer of 10000 m around the point. This 
# we do inside the extract function, which in contrast to the gBuffer-function
# that we used last week allows us to define correctly distance in meters also
# for our latlong-projected data. In the extract function we can also directly
# specify a function to generate summary statistics for the raster points within 
# the buffer (e.g. median, standard deviation).
cities.buf.med.ele <- list()
cities.buf.sd.ele <- list()
l.cit <- length(pointshape.shp)
for (i in 1:l.cit) {   
  cities.buf.med.ele[[i]] <-  extract(nepal.ele.mask,pointshape.shp[i,],buffer=10000,
                                      fun=median)
  cities.buf.sd.ele[[i]] <-  extract(nepal.ele.mask,pointshape.shp[i,],buffer=10000,
                                     fun=sd)
}

# Generate dataset with outcomes. This you could, e.g. join, with the shapefile
# by using the spCbind or other commands
cities.ele.dat <- data.frame(mydata,unlist(cities.ele),
                             unlist(cities.buf.med.ele),
                             unlist(cities.buf.sd.ele))
colnames(cities.ele.dat) <- c("Name","CNTR","Lon","Lat","Ele.Exact","Ele.Median.Buf",
                              "Ele.Sd.Buf")
cities.ele.dat
