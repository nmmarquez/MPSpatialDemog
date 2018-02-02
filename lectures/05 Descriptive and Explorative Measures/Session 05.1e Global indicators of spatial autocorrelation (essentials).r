################################################################################
#                                                                              #                         
# IDEM 156 - Spatial Demography Course 2018                                 #                          
# Calculating Global Indicators of Spatial Autocorrelation (essentials)        #
# Sebastian Klüsener, MPIDR                                                    #                         
#                                                                              #                         
################################################################################

# Erase all objects in workspace
rm(list=ls(all=TRUE))

# Load libraries
library(spdep)
library(maptools)
library(rgdal)

# Set your working directory
# If you need to change the working directory, you can use the following code:
# Set drive
main <- c("N")
# Set path to session folders
path <- ":/IMPRSD/IDEM 156 Spatial"
# Define session folder
path2 <- "/02 Sessions/05 Descriptive and Explorative Measures"
# Set working directory
setwd(paste(main,path,path2,sep=""))


################################################################################
# 1) Import and prepare data                                                   #
################################################################################ 

# Read shapefile of Germany
shape.shp <- readShapePoly('German_Districts_0406', IDvar="DISTID")
shapeogr.shp <- readOGR(".", "German_Districts_0406")
proj4string(shape.shp) <- proj4string(shapeogr.shp)
# Reprojected version of the shapefile which has a longlat projection
shape.ll <- spTransform(shape.shp, CRS("+proj=longlat +datum=WGS84"))

plot(shape.shp)                   
names(shape.shp)

# Define ID
id <- shape.shp$DISTID

# Centroids of the reprojected Shapefile with LotLang-projection (will be used
# for distance calculatons)
coords.ll <- coordinates(shape.ll)

# In order to "attach" data from a shapefile attribute table, we have to call it
# with "@data"
attach(shape.shp@data)
               
# Create a contiguity-based first order queen-neighborhood weight matrix 
nb.FOQ <- poly2nb(shape.shp, queen=TRUE)
summary(nb.FOQ)

# In the case of Germany, the poly2nb function would detect for the island-
# district of Rügen no neighbors at it is not directly bordering any other 
# district. In order to ensure that all district have neighbors, which is 
# helpful (but notnecessary) to calculate spatial autocorrelation statistics, 
# we import a gal-file which we produced with the algorithms in file 3.3 
# (Neighborhood weight files) to add neighbor-links 
# between Rügen and the districts of Stralsund und Vorpommern.
nb.FOQ.cor <- read.gal('FOQ_German_Districts_corrected.nb')

# Create a nearest neighbor weight matrix (in this exmaple based on 5 
# nearest neighbors)
l.5NN <- knearneigh(coords.ll, k=5, longlat=T)
nb.5NN <- knn2nb(l.5NN, row.names=id)

# Create a distance based matrix (distance between centroids)
# Distance of 20 km
d20km <- dnearneigh(coords.ll, 0, 20, row.names = id, longlat=T)
# Distance of 50 km
d50km <- dnearneigh(coords.ll, 0, 50, row.names = id,longlat=T)

# Turn neighborhood objects in listwise objects.
# If style is not defined, the function will use "W", which implies that the 
# weights will be row-standardized. B implies that each neighboring region gets 
# a weight of 1. For the uncorrected contiguity based matrices, in which Rügen 
# does not have any neighbors in the German file, we need to set zero.policy to 
# TRUE so that the nb2listw-command allows regions with no neighbors. The same 
# true for the distance 20 neighborhood file, which also contains for some 
# regions no neighbors
nb.FOQ.lw.W <- nb2listw(nb.FOQ, style="W",zero.policy=T)
nb.FOQ.lw.B <- nb2listw(nb.FOQ, style="B",zero.policy=T)
nb.FOQ.cor.lw.W <- nb2listw(nb.FOQ.cor, style="W")
nb.FOQ.cor.lw.B <- nb2listw(nb.FOQ.cor, style="B")
nb.d20km.lw.W <- nb2listw(d20km,zero.policy=T)
nb.d50km.lw.W <- nb2listw(d50km)
nb.5NN.lw.W <- nb2listw(nb.5NN)


################################################################################
# 2) Specify for which variable and neighborhood definition you want to        #
#    calculate the Moran's I                                                   #
################################################################################

# Define the variable for which you would like to calculate the Moran's I
# Derive variables names from attribute table and choose one (e.g. Share Non-
# Marital Births (SNM05) or Total Fertility Rate (TFR05)
names(shape.shp)

varofint <- SNM05
# varofint <- TFR05

# Define, which of the neighborhood definition you want to use
# nb.FOQ.lw.W     - Contiguity - First order queen with row-standardized weights
# nb.FOQ.lw.B     - Contiguity - First order queen 
# nb.FOQ.cor.lw.W - Contiguity - First order queen with row-standardized weights
#                   (Corrected version in which all German districts have at 
#                   least one neighbor)
# nb.FOQ.cor.lw.B - Contiguity - First order queen 
#                   (Corrected version in which all German districts have at 
#                   least one neighbor)
# nb.d20km.lw.W   - Distance between centroids less than 20km
# nb.d50km.lw.W   - Distance between centroids less than 50km
# nb.5NN.lw.W     - Five nearest neighbors (based on distance between centroids)

nb.listw <- nb.FOQ.cor.lw.W
#nb.listw <- nb.5NN.lw.W

# Define number of permutations that should be performed for the Monte Carlo-
# based permutation test.
npermutations <- 9999


################################################################################
#                                                                              #
# 3) Calculate Moran's I test                                                  #
#                                                                              #
################################################################################

# Distribution-based - Assumption of Randomization (adjusted for kurtosis)
moranI.db <- moran.test(varofint, nb.listw, alternative="greater")
moranI.db

# Permutation test for Moran's I statistic
moranI.mc <- moran.mc(varofint, nb.listw, npermutations)
moranI.mc


################################################################################
# 4) Application of Moran's I test in modelling procedures                     #
#    Simple Models testing for associations between unemployment and           #
#    demographic attributes                                                    #
################################################################################

# Define variables and models
y1 <- SNM05
y2 <- TFR05
x1 <- UER05

# Define model formula
model1 <- y1 ~ x1
model2 <- y2 ~ x1

# Calculate Moran's I test for depdendent variables
moranI.db.y1 <- moran.test(y1, nb.listw)
moranI.db.y2 <- moran.test(y2, nb.listw)

# Calculate models
model1 <- lm(model1)
model2 <- lm(model2)
summary(model1)
summary(model2)

# Run lm.morantest to test, whether residuals are autocorrelated
moranI.residuals.model1 <- lm.morantest(model1, nb.listw)
moranI.residuals.model2 <- lm.morantest(model2, nb.listw)

# Contrast Moran's I of dependent variables and models
moranI.db.y1
moranI.residuals.model1
moranI.db.y2
moranI.residuals.model2
