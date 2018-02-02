################################################################################
#                                                                              #                          
# IDEM 156 - Spatial Demography Course 2018                                    #                          
# Geographically Weighted Regression vs. Random Intercept Model                #                          
# Sebastian Klüsener, MPIDR                                                    #                          
#                                                                              #                          
################################################################################

# Erase all objects in workspace
rm(list=ls(all=TRUE))

# Load libraries
library(spgwr)
library(RColorBrewer)
library(spdep)
library(rgdal)
library(maptools)
library(INLA)

# Set your working directory
# If you need to change the working directory, you can use the following code:
# Set drive
main <- c("N")
# Set path to session folders
path <- ":/IMPRSD/IDEM 156 Spatial"
# Define session folder
path2 <- "/02 Sessions/08 Spatial Multi-level Modeling"
# Set working directory
setwd(paste(main,path,path2,sep=""))


################################################################################
# 1) Import and prepare data                                                   #
################################################################################

# Open shapefile
shape.shp <- readShapePoly('2004_06', IDvar="KREIS_KENN")
shapeogr.shp <- readOGR(".", "2004_06")
proj4string(shape.shp) <- proj4string(shapeogr.shp)
shape.shp <- spTransform(shape.shp, CRS("+proj=longlat +datum=WGS84"))

# Derive information on internal id and centroids. In this case we turn the id-
# values first into a numeric vector, which is easier to use, if we want to
# match data to the shapefile.
head(shape.shp$KREIS_KENN)
id <- as.numeric(paste(shape.shp$KREIS_KENN))/1000

# Define some standard objects
n <- length(id)
sq <- c(1:n)
sq_text <- paste(sq)
coords <- coordinates(shape.shp)

# Upload data file
mydata <- read.table(file="fertdata2005.csv",  sep=",", header=TRUE)

# Match data to shapefile
data.id <- mydata$ID

# Check, whether vector of data.ids and id of shapefile are of the same lenght
length(data.id) == length(id)

# Match by ID-column, which has to be similar in the shapefile and the datafile
o <- match(id, data.id)
mydata1 <- mydata[o,]
row.names(mydata1) <- shape.shp$KREIS_KENN
shape.shp <- spCbind(shape.shp, mydata1)

# Derive variable names
names(shape.shp)

# Create an autonumber
shape.shp$an <- c(1:length(shape.shp))

# In order to "attach" data from a shapefile attribute table, we have to call it
# with "@data"
attach(shape.shp@data)

nb.FOQ <- poly2nb(shape.shp, row.names=shape.shp$an , queen=TRUE)
lw.FOQ <- nb2INLA("FOQ_INLA",nb.FOQ)

# Defining Variables      
# Dependent Variable: Share Non-Marital Births
y <- SNM05
# Covariate 1: Female Employment Rate
x1 <- FR155005    
# Dummy East
dum_east <- DUM_E 

# Normal linear model
summary(lm(y~x1))

# Random-intercept model
formula <- y~x1+f(an,model="iid")
random_intercept_model <- inla(formula,family="gaussian",data=data.frame(an,y,x1))
summary(random_intercept_model)

# Random-slope model
formula <- y~x1+f(an,x1,model="iid")
random_slope_model <- inla(formula,family="gaussian",data=data.frame(an,y,x1))
summary(random_slope_model)

# Random-slope with spatial effects
formula <- y~x1+f(an,x1,model="besag",graph="FOQ_INLA")
random_model <- inla(formula,family="gaussian",data=data.frame(an,y,x1))
summary(random_model)

