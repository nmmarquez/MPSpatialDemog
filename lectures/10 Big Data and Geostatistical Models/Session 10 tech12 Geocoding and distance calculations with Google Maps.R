################################################################################
#                                                                              #                          
# IDEM 156 - Spatial Demography Course 2018 - Technical File 12                #                          
# Derive Location points and calculate distances via Google Maps               #                          
# Please note that you agree to Google's terms of use, which might             #
# interfere with your research plans:                                          #
# https://developers.google.com/maps/terms?hl=de-AT                            #                                                                             
#                                                                              #
# Sebastian Kl?sener, MPIDR                                                    #                          
#                                                                              #                          
################################################################################

# Erase all objects in workspace
rm(list=ls(all=TRUE))


# Install libraries
#install.packages(c("ggmap","osmar","plyr","OpenStreetMap","plyr","scales"))

# Open libraries
library(ggmap)
library(OpenStreetMap)
library(maptools)
library(rgdal)
library(rgeos)
library(plyr)
library(scales)

################################################################################
#                                                                              #
# 1) Extracting information from google maps through API                       #
#                                                                              #
################################################################################ 

# Limited to 2500 request per day in each category)

# Geocodes 
geocodeQueryCheck()
# Road distances/travel times
distQueryCheck()
# Travel routes
routeQueryCheck()


################################################################################ 
#                                                                              #  
# 1.1) Extract geocode for a city or an address                                #
# 1.1.1) For a city                                                            #
#                                                                              #
################################################################################ 

rostock.geo <- geocode("Rostock")
rostock.geo

# After running the query, we have 2499 left for today
geocodeQueryCheck()

# What is according to google the "center" of Rostock?
# For this we use the quickmap command with zoom defining how far we want to 
# zoom in.
qmap("Rostock", zoom = 18,maptype="satellite") 

# Let us zoom a little bit out. As an alternative to the qmap-function, we use 
# now a combination of the get_googlemap and ggmap-functions that provide more
# options to specify commands.
mapros <- get_googlemap("Rostock", zoom = 14,maptype="satellite", 
                        markers=rostock.geo,scale=2) 
ggmap(mapros,extent ="device")  

 
# Different background maptypes are available 
# qmap("Rostock", zoom = 14,maptype="road")  
# qmap("Rostock", zoom = 14,maptype="hybrid") 
# qmap("Rostock", zoom = 14,maptype="terrain") 
# Terrain is more interesting in mountain areas
# qmap("Matterhorn", zoom = 14,maptype="terrain") 
# Even artistic maps from Stamen
# qmap("Moscow", zoom = 14,maptype="toner",source="stamen") 
# This last option seems currently not to be supported and returns an error:
# qmap("Moscow", zoom = 14,maptype="watercolor",source="stamen") 

################################################################################ 
#                                                                              #
# 1.1.2) For an address                                                        #
#                                                                              #
################################################################################ 

# Be aware that spelling errors/places with similar names
# might cause many errors. However, small deviations in the
# spelling do not necessarily cause errors.
# Example for the address of the Administration School
nes.geo <- geocode("100 Novaya Street,Skolkovo,Moscow")
nes.geo
nes.geo1 <- geocode("100 Novaya Street,Skolkovo")
nes.geo1
nes.geo2 <- geocode("100 Novaya Street,Skolovo,Moscow")
nes.geo2
nes.geo3 <- geocode("100 Novaa Street,Skolovo,Moscow")
nes.geo3

maploc <- get_googlemap("Novaya Street 100,Skolkovo,Moscow", zoom = 18,maptype="satellite", 
                        markers=nes.geo,scale=2) 
ggmap(maploc,extent ="device")


stp.geo <- geocode("Saint Petersburg, Russia")
stp.geo
qmap(location = "Saint Petersburg, Russia", zoom = 14, source = 'osm')

ber.geo <- geocode("Berlin")
moscow.geo <- geocode("Moscow")
  
# Calculate spatial distance (spherical) between "NES", "Moscow",
# and "Saint Petersburg"
spdist <- data.frame(rbind(nes.geo,moscow.geo,stp.geo,ber.geo))
rownames(spdist) <- c("NES","Moscow","Saint Petersburg","Berlin")

#distances <- read.table("Distance.csv",sep=",",head=T)
distances <- spDists(as.matrix(spdist),longlat=T)
rownames(distances) <- rownames(spdist)
colnames(distances) <- rownames(spdist)
distances
#write.table(distances,"Distance.csv",sep=",",row.names=F)

################################################################################ 
#                                                                              #
# 1.2) Derive on ground distances and travel times between cities (depending   #
#      on existent road infrastructure)                                         #
#                                                                              #
################################################################################

from <- c("Berlin", "Saint Petersburg,Russia")
to <- c("Moscow")

# Get distances on road and travel times by different modes
bycar <- mapdist(from, to, mode = "driving")
bybike <- mapdist(from, to, mode = "bicycling")
onfoot <- mapdist(from, to, mode = "walking")

# Distance on road differs by mode of transport
outcome <- data.frame(bycar[,c(1:2,4)],bybike[,4],onfoot[,4])
colnames(outcome) <- c("from","to","car_(km)","bike_(km)","on_foot_(km)")
outcome

# Travel times, of course, as well
outcomet <- data.frame(bycar[,c(1:2,8)],onfoot[,8])
colnames(outcomet) <- c("from","to","car_(h)","on_foot_(h)")
outcomet

# We can also derive road distances and travel times for addresses
from <- c("Red Square, Moscow")
to <- c("Potsdamer Platz 1, Berlin")
bycar1 <- mapdist(from, to, mode = "driving")

# Contrast the result between spherical distance, on ground city distances 
# and on ground, address distance
# Rostock and Berlin  
dist.B <- data.frame(matrix(c("Moscow","Berlin",NA,round(distances[2,4],3)
                              ,rep(NA,4)),ncol=8))
colnames(dist.B) <- colnames(bycar)
compare <- rbind(bycar[1,],bycar1,dist.B)[,-c(3,5:8)]
compare$disttype <- c("on ground - center", "on ground - address", 
                      "spherical distance - independent of roads")
compare 


################################################################################ 
#                                                                              #
# 1.3) Derive routes                                                           #
#                                                                              #
################################################################################

# Import a dataset with some addresses
rou.addr <- read.table("Addresses.csv",sep=",",head=T)

# We define a type column as we want all routes to by car.
rou.addr$type <- c("driving")
rou.addr

# Performing route calulations
# We first define a list in which we store the results
lrad <- length(rou.addr[,1])
my.routes <- list()

for (i in 1:lrad) {
# Defining start point, end point and type for route i
    from <- paste(rou.addr[i,]$from)
    to <- paste(rou.addr[i,]$to)
    type <- rou.addr[i,]$type
# Generate route i
    sinroute <- route(from, to,
    mode = type, alternatives=FALSE,
    structure = c("route"))
# Store route i
my.routes[[i]] <- sinroute
}
my.routes

# Transform the routes in a line shapefile 
lmrou <- length(my.routes)
lines <- list()
for (i in 1:lmrou) {
  toline <- my.routes[[i]]
  coordsl <- data.frame(toline$lon,toline$lat)
  line <-Line(coordsl)
  lines[[i]] <- Lines(list(line), ID=paste(i)) 
}
route.shp <- SpatialLines(lines)

# Open shapefile
shape.shp <- readShapePoly('2004_06.shp', IDvar="KREIS_KENN")
shapeogr.shp <- readOGR(".", "2004_06")
proj4string(shape.shp) <- proj4string(shapeogr.shp)

# Adjust projection
shape.ll <- spTransform(shape.shp, CRS("+proj=longlat +datum=WGS84"))

# Plot routes
plot(shape.ll,border="grey25",lwd=1)
plot(shape.ll,col=alpha("white",0.5),lty=0,add=T)
lines(route.shp[1],col=c("red"),lwd=2)
lines(route.shp[2],col=c("navy"),lwd=2)

