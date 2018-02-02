################################################################################
#                                                                              #                          
# IDEM 156 - Spatial Demography Course 2018                                    #                          
# Introduction to Mapping Functions in R                                       #                          
# Dalia A. Conde and Fernando Colchero, University of Odense                   #
# with adaptations by Sebastian Kl?sener, MPIDR                                #
#                                                                              #                          
################################################################################

# Erase all objects in workspace
rm(list=ls(all=TRUE))

# Load libraries
library(sp)
library(maptools)
library(rgdal)
library(RColorBrewer)


################################################################################
# 1) Import data and draw maps of the world                                    #
################################################################################

# If you need to change the working directory, you can use the following code:
# Set driv
# Load a shapefile type polygon
# Load a world map of 1992 (ESRI map)
world <- readShapePoly('CNTRY92', IDvar="FIPS_CODE")
worldogr <- readOGR(".", "CNTRY92")
proj4string(world) <- proj4string(worldogr)

# Print summary attributes for the 'world' layer using the 
# built-in function 'summary()':
summary(world)

# Plot the 'world' layer with the function 'plot()':
plot(world, col="grey98", border="grey25",bg="#f0f8ff")

# Load the linear shapefile 'RIVERS':
rivers <- readOGR(".", "RIVERS")
summary(rivers)

# Use the built-in function 'lines()' to plot the rivers:
lines(rivers, col = "dodgerblue2", lwd = 2)

# To improve visibility, use 'lwd' to make rivers thinner:
plot(world, col="grey98", border="grey25",bg="#f0f8ff")
lines(rivers, col = "dodgerblue2", lwd = 1)

# Load the point shapefile 'CITIES':
cities <- readOGR(".", "CITIES")
summary(cities)

# Use the built-in function 'points()' to add the cities to the main plot: 
points(cities, pch=21, col="brown", bg="brown2", cex=0.5)

# Plot cities by population size and change color of the world layer.
# From the cities layer get the POPULATION attribute
pop <- cities$POPULATION 

# Change to 'NA' all cities with population less than or equal zero:
pop[pop <= 0] = NA

# Plot the resulting layers:
plot(world, col="grey98", border="grey25",bg="#f0f8ff")
lines(rivers, col="dodgerblue2", lwd=1)

# Use cex to indicate the size of the dots according to the population size
points(cities, pch=21, col="brown", bg="brown2", 
       cex=pop / max(pop, na.rm = TRUE))

# More accurate: Use pi to create circles representing population size
plot(world, col="grey98", border="grey25",bg="#f0f8ff")
lines(rivers, col = "dodgerblue2", lwd = 1)
points(cities, pch=21, col="brown", bg="brown2", lwd=0.5,
 cex = sqrt((pop/10000000)/pi))

################################################################################
# 2) Map subregions of the world                                               #
################################################################################

# Map of a subregion of the world
# First plot command defines the subregion to be plotted.
# The plot window will be focused on that subregion.
plot(world[world$NAME=="Thailand",],bg="#f0f8ff")
plot(world, col="grey98", border="grey25", lwd=1.5,add=T) # oh shit thats tricky
plot(world[world$NAME=="Thailand",], col="lightgreen", border="grey25", 
     lwd=1.5, add=T)
lines(rivers, col = "dodgerblue2", lwd=2)
points(cities, col="brown", bg="brown2", pch = 21, lwd=2,
 cex = sqrt(pop/100000/pi))
text(coordinates(cities)[,1],
     coordinates(cities[,2])+0.35+sqrt(pop/10000000/pi),cities$NAME)
