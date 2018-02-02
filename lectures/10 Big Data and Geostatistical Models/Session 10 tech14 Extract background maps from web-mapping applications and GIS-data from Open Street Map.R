################################################################################
#                                                                              #                          
# IDEM 156 - Spatial Demography Course 2018 - Technical File 14                #                          
# - Extract background maps from web-mapping applications such as              #  
#   OpenStreetMap                                                              #
# - Exract GIS-data from OpenStreetMap                                         #                                                                             
# Sebastian Kl?sener, MPIDR                                                    #                          
#                                                                              #                          
################################################################################

# Erase all objects in workspace
rm(list=ls(all=TRUE))

library(OpenStreetMap)
library(osmar)
library(scales)
library(maptools)
library(rgdal)


# Please intert path to your working directory
# Set your working directory
# If you need to change the working directory, you can use the following code:
# Set drive

################################################################################
# 1) Using the OpenStreetMap-library to derive background maps                 #
#    on which you can plot your shapefiles with data                           #
################################################################################

# Example of an extraction of an OSM map using the OpenStreetMap-library.
# To extract the map your first define the minimum and maximum longitude and 
# latitude-values of your plot. In this example, we chose values that define the 
# area in and aroundGermany
min_latitude <- 47
max_latitude <- 56
min_longitude <- 5
max_longitude <- 16

# Different background maps are available.
nm <- c("osm", "maptoolkit-topo", "bing", "stamen-toner",
        "stamen-watercolor","stamen-terrain","esri","esri-topo")
      

# Dividing your plot area in eight sections to plot examples of the different 
# background maps
par(mfrow=c(2,4))
# Using for-loop to print the eight different backgrounds defined in nm-vector
for(i in 1:length(nm)){
   map <- openmap(c(max_latitude,min_longitude),
                 c(min_latitude,max_longitude),
                 minNumTiles=3,type=nm[i])
   plot(map)
}

# Open Shapefile of Germany and assign projection information to it.
shape.shp <- readShapePoly('2004_06.shp', IDvar="KREIS_KENN")
shapeogr.shp <- readOGR(".", "2004_06")
proj4string(shape.shp) <- proj4string(shapeogr.shp)

# Reproject that shapefile to longlat-information
shape.ll <- spTransform(shape.shp, CRS("+proj=longlat +datum=WGS84"))


# We chose one background map which we liked and map now just one plot of 
# Germany using this background map.
min_latitude <- 47
max_latitude <- 56
min_longitude <- 5
max_longitude <- 16

# The higher the minNumTiles (resolution), the higher the detail, at which pixels
# are extracted and the sharper the map. However, be aware that in increasing the
# detail you are also increasing the time it takes to upload this map into R.
par(mfrow=c(1,1))
map <- openmap(c(max_latitude,min_longitude),
               c(min_latitude,max_longitude),
               minNumTiles=8,type="bing")

# The background map is coming in longlat-projection. However, in order to be 
# able to plot the background files in combination with shapefiles which we 
# uploaded in R with functions such as readShapePoly or readOGR, we need to run 
# the openproj-command over the background-map object. 
mappr <- openproj(map,projection="+proj=longlat")

# Now we can do combined plots of the background map and our shapefiles in 
# longlat projection:
plot(mappr)
plot(shape.ll,border="grey75",col="transparent",add=T)

# As the map of Germany does not look very nice in longlat-projection, we use 
# instead the Cartesian projection in which our German shapefile is originally
# projected. This is an UTM-projection which quite well preserves area, shapes,
# and distances for spatial objects located within Germany
# Here we reproject the longlat-map into the projection of our shapeogr.shp
# uploaded in line 53, which represents the original projection of our Germany
# shapefile. The resulting maps looks much nicer.
map_coord <- openproj(map,projection=proj4string(shapeogr.shp))
plot(map_coord)
plot(shape.shp,add=T,border="grey75",col="transparent")

# Another way to make even longlat-projected maps of Germany look as we commonly 
# know them, is to first plot the longlat-projected shape.ll-object, which is then 
# in contrast to the longlat-projected background map automatically plotted in a 
# mercator-based projection.
plot(shape.ll)
plot(mappr,add=T)
plot(shape.ll, add=T,border="grey75",col="transparent")


# Now we would like to extend the area of the background map so that it fills 
# the whole background of the plot area. For this we just extend the area 
# plotted by modifying our min and max latitude and longitude values a bit.
min_latitude <- 47-7
max_latitude <- 56+4
min_longitude <- 5-4
max_longitude <- 16+4

par(mfrow=c(1,1))
map <- openmap(c(max_latitude,min_longitude),
               c(min_latitude,max_longitude),
               minNumTiles=8,type="mapquest-aerial")
map_coord <- openproj(map,projection=proj4string(shapeogr.shp))

# Depending on our plot window size, now the whole background should be filled 
# by the background map. If not, we would further need to extend our plot 
# window by re-specifying the min and max longitude and latitude values in 
# lines 111-114 and run the code again from these lines until line 130.
plot(shape.shp)
plot(map_coord,add=T)
# The alpha command allows to define the transparency of a fll- or background-
# color. You can choose values between 0 (not printed) up to 1 (full colors).
plot(shape.shp,add=T,border="grey75",col=alpha("white",0.15))


################################################################################
# 2) Fetching of geodata from OpenStreetMap                                    #
################################################################################

# Now instead of loading a pixel background map, we want to fetch polygons, 
# line- or point objects from OpenStreetMap and import them into R. 

# If we map Germany completely, we would not see which polygon-files are 
# available
par(mfrow=c(1,1))
map <- openmap(c(max_latitude,min_longitude),
               c(min_latitude,max_longitude),
               minNumTiles=8,type="osm")
plot(map)

# But the further we zoom in by defining very small distances between the 
# minimum and maximum longitudes and latitudes or our plotting window, we get 
# an overview over available  elements such as buildings, streets, public 
# transport, infrastructure etc.
# Enclosed an example for the area around the MPIDR
map <- openmap(c(54.096,12.108),
               c(54.092,12.115),
               minNumTiles=8,type="osm")
plot(map)

# Now we are using functions of the osmar library to fetch spatial data out of 
# OpenStreetMap into R.
# Defining the api we want to use (in this case OpenStreetMap API)
api <- osmsource_api()

# Defining the area in which we want to fetch data. For this you can either use
# the corner_bbox-command in which you are againdefining min- and max-values for 
# the longitudes and latitudes (?corner_bbox) or the center_bbox command where 
# you are fetching data in an area around a central location. In this example we
# are catching data in box of 50*50 which is centered on the MPIDR.
bb <- center_bbox(12.111, 54.094037, 25, 25)

# Code to derive data for our defined fetching area.
mpidr <- get_osm(bb, source = api)

# We have just captured all geodata in the fetching area which we defined.
# The plot looks like the MPIDR-building. So everything seems to be fine.
plot(mpidr)

# Now we want to isolate in the spatial information which we fetched all those
# points that belong to buildings. In this example it is just the building of the
# MPIDR. With the following lines of code we are turning it into a polygon 
# shapefile object.
bg_ids <- find(mpidr, way(tags(k == "building")))
bg_ids <- find_down(mpidr, way(bg_ids))
bg <- subset(mpidr, ids = bg_ids)
bg_poly <- as_sp(bg, "polygons")
bg.shp <- bg_poly
# As we know that this is the MPIDR building, we assign MPIDR as name to it
bg.shp$Name <- c("MPIDR")
# But it also has a unique ID:
bg.shp$id

# Here we first plot the complete fetched data and then the newly derived
# shapefile of the MPIDR over it.
plot(mpidr)
plot(bg.shp,add=T,col="grey90")


################################################################################
# 3) Example for a overlay/intersection operation, in which we intersect a     #
#    polygon-shapefile with a point-shapefile                                  #
################################################################################

################################################################################
# 3.1) Preparation of the data for the point-shapefile                         #
################################################################################

# Location dots of the course participants in the building.
perdat <- read.table("Persons_comp.csv",sep=",", head=T)
perdat <- perdat[perdat$Print=="TRUE",]
coordinates <- data.frame(perdat[,2],perdat[,3])
# Turning this location information in a point-shapefile.
per.shp <- SpatialPointsDataFrame(coordinates, data=perdat)
# Plotting the shapefile
plot(mpidr)
plot(bg.shp,add=T,col="grey90")
plot(per.shp,add=T,col="red",pch=20)

# The new shapefile we created does not have projection information connected to
# it. 
proj4string(per.shp)
# The bg.shp which we created from our OSM-data is in longlat projection, which
# is the same projection in which also the participant location dots are that
# we used to create the point-shapefile.
proj4string(bg.shp)
# Thus, we can just assign the projection information of bg.shp to per.shp.
proj4string(per.shp) <- proj4string(bg.shp)


################################################################################
# 3.2) Performing the overlay/intersect operation                              #
################################################################################

# Here we are intersecting points with polygons to see which points are within 
# which polygons (in case some of the points are overlapping with the polygons). 
# In this example, we are checking which points are in the MPIDR-building. With 
# the over-operation we are deriving a data-frame in which the rows represent 
# the participants, while it contains in the columns information with which 
# polygons the participants points are operlapping
over <- over(per.shp, bg.shp)
over

# In this case, all points are in the MPIDR-building, and we assign one column
# of the  over-Object to per.shp Pointshapefile of the participants.
per.shp$bg <- over$Name
# Now we just need to see which of the participants are in the building...
present <- which(per.shp$bg=="MPIDR")
# and print there name based on information from the column giving the 
# participants names "$Person".
# Who's is in the house?
print(paste(per.shp$Person[present]))


################################################################################
#                                                                              #
# 4) Example where we fetch elements from a larger area, derive different kinds#
#    of shapefiles (buildings, ways, landuse) and connect address data to the  #
#    buildings.                                                                #
#                                                                              #
################################################################################

################################################################################
# 4.1) Fetching spatial objects, turning them in shapefiles and linking        #
# attributes from OpenStreetMap to these shapefiles (e.g. addresses)           #
################################################################################

# In case we want to increase the area from which we fetch data, we just 
# increase the width and height of the area from which we are fetching.
#api <- osmsource_api()
# Rostock
bb <- center_bbox(12.111, 54.094037, 400, 400)
# Bolschoi Theatre, Moscow
#bb <- center_bbox(37.618183, 55.760664, 400,400) 
# Area in New York
#bb<- center_bbox(-74.005941, 40.712784, 250, 250)
#areaplus <- get_osm(bb, source = api)
# Area in New Dehli
#bb<- center_bbox(77.220724, 28.632244, 250, 250)
#areaplus <- get_osm(bb, source = api)
# Follina, Italy
#bb<- center_bbox(12.118912, 45.952778, 250, 250)
# Munich
#bb<- center_bbox(11.581981, 48.135125, 250, 250)
# Barcelona
#bb<- center_bbox(2.173403, 41.385064, 250, 250)


areaplus <- get_osm(bb, source = api)


# Webpage with list of potentially available keys:
# http://wiki.openstreetmap.org/wiki/Category:En_key_descriptions
overview_elements <- data.frame((areaplus$way$tags))
# Upper hierachy
table(overview_elements$k)
# Lower hierachy
table(overview_elements$v)
# Some example, what information could be derived about the 
# buildings or in terms of landuse in the area. This could be lined 
overview_elements[overview_elements$k=="building",]
overview_elements[overview_elements$k=="highway",]
overview_elements[overview_elements$k=="name",]
overview_elements[overview_elements$k=="addr:postcode",]
overview_elements[overview_elements$k=="addr:street",]
overview_elements[overview_elements$k=="addr:housenumber",]
overview_elements[overview_elements$k=="addr:city",]
overview_elements[overview_elements$k=="landuse",]
overview_elements[overview_elements$k=="building:height",]
overview_elements[overview_elements$k=="building:material",]
#overview_elements[overview_elements$k=="amenity",]
#overview_elements[overview_elements$k=="website",]
# overview_elements[overview_elements$k=="email",]


# Derive polygon-shapefiles of buildings
bg_ids <- find(areaplus, way(tags(k == "building")))
bg_ids <- find_down(areaplus, way(bg_ids))
bg <- subset(areaplus, ids = bg_ids)
bgplus.shp <- bg_poly <- as_sp(bg, "polygons")

# Derive line-shapefiles of streets and ways
hw_ids <- find(areaplus, way(tags(k == "highway")))
hw_ids <- find_down(areaplus, way(hw_ids))
hw <- subset(areaplus, ids = hw_ids)
hwplus.shp <- hw_line <- as_sp(hw, "lines")

# Derive polygon-shapefiles of public transport infrastructure
rw_ids <- find(areaplus, way(tags(k == "railway")))
rw_ids <- find_down(areaplus, way(rw_ids))
rw <- subset(areaplus, ids = rw_ids)
rwplus.shp <- rw_line <- as_sp(rw, "lines")

# Derive polygon-shapefiles of postal_code data
pc_ids <- find(areaplus, way(tags(k == "postal_code")))
pc_ids <- find_down(areaplus, way(pc_ids))
pc <- subset(areaplus, ids = pc_ids)
pcplus.shp <- pcw_line <- as_sp(pc, "lines")

# Derive information on landuse and link it to the polygon-shapefile
# showing the extent of areas with the same landuse (e.g. residential.
# commercial)
lu_ids <- find(areaplus, way(tags(k == "landuse")))
lu_ids <- find_down(areaplus, way(lu_ids))
lu <- subset(areaplus, ids = lu_ids)
luplus.shp <- as_sp(lu, "polygons")
lud <- data.frame(overview_elements[overview_elements$k=="landuse",])
o <- match(luplus.shp$id,lud$id)
lud1 <- lud[o,]
luplus.shp$landuse <-lud1$v  

# Assigning colors to some landuse values
col <- rep(NA,length(luplus.shp))
col[luplus.shp$landuse=="grass"] <- c("lightgreen")
col[luplus.shp$landuse=="commercial"] <- c("aliceblue")
col[luplus.shp$landuse=="construction"] <- c("lightgrey")
col[luplus.shp$landuse=="residential"] <- c("mistyrose")

# Derive information on addresses and links this information to
# the buildings-shapefile
str <- overview_elements[overview_elements$k=="addr:street",]
housn <- overview_elements[overview_elements$k=="addr:housenumber",]
city <- overview_elements[overview_elements$k=="addr:postcode",]
pco <- overview_elements[overview_elements$k=="addr:city",]
os <- match(bgplus.shp$id,str$id)
str1 <- str[os,]
bgplus.shp$street <- str1$v
oh <- match(bgplus.shp$id,housn$id)
housn1 <- housn[oh,]
bgplus.shp$housn <- housn1$v
oc <- match(bgplus.shp$id,city$id)
city1 <- city[oc,]
bgplus.shp$city <- city1$v
opco <- match(bgplus.shp$id,pco$id)
pco1 <- pco[opco,]
bgplus.shp$pco <- pco1$v

################################################################################
# 4.2) Outcome plot                                                            #
################################################################################

# Now we have elements from a larger area around the MPIDR building.
plot(bgplus.shp,col="white",lty=0)
# Landuse information
plot(luplus.shp,col=col,lty=0,add=T)
#plot(openproj(map),add=T)
# Buildings
plot(bgplus.shp,col="red",add=T)
# Streets
plot(hwplus.shp,col="grey75",add=T)
# Tram/Railway
plot(rwplus.shp,col="black",lwd=1,add=T)

# Here we display the addresses of the building, in case they are available.
# Otherwise, NAs would be shown.
centroids <- coordinates(bgplus.shp)
# In case you would like to have the house number before the street name,
# you would need to switch "bgplus.shp$street" and "bgplus.shp$housn".
text(centroids[,1],centroids[,2],paste(bgplus.shp$street,bgplus.shp$housn,sep=" "),cex=0.5)

# In case there is no address data available. you could use the unique 
# building-IDs to it.
bgplus.shp$id

# Save building-shapefile with address data.
writeOGR(bgplus.shp, ".",'extracted_buildings', driver="ESRI Shapefile")

# Save landuse-shapefile with landuse data.
writeOGR(luplus.shp, ".",'extracted_landuse', driver="ESRI Shapefile")
