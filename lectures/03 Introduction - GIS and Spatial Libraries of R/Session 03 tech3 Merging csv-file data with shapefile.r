################################################################################
#                                                                              #                          
# IDEM 156 - Spatial Demography Course 2018 - Technical File 3                 #                          
# Code to match data to a shapefile                                            #                          
# Sebastian Kl?sener, MPIDR                                                    #                          
#                                                                              #                          
################################################################################


# Erase all objects in memory
rm(list = ls(all = TRUE))

# Load libraries
library(maptools)
library(rgdal)

# # Please note that if you launch R-studio by double-clicking on an R-codefile, 
# # R studio sets the working directory to the folder in which the R-file that you
# # double-clicked is located. 
# # Command to check working directory 
# getwd()
# 
# # If you need to change the working directory, you can use the following code:
# # Set drive
# main <- c("N")
# # Set path to session folders
# path <- ":/IMPRSD/IDEM 156 Spatial"
# # Define session folder
# path2 <- "/02 Sessions/03 Introduction - GIS and Spatial Libraries of R"
# # Set working directory
# setwd(paste(main,path,path2,sep=""))


################################################################################
# 1) Import and prepare data                                                   #
################################################################################ 

# Import shapefile
# We are using here two different functions to open the same shapefile, which
# allows us to use different advantages of the two functions. The first
# funcion "readShapePoly" allows to define an ID-variable, which is helful
# for the linking, as the shapefile will be sorted in the order of this ID-
# variable. The "readOGR"-function offers the advantage that is automatically,
# uploads the projection inforrmation in case a prj-file is connected to the
# shapefile.
shape.shp <- readShapePoly('2008_w_fert_data', IDvar="GS")
shapeogr.shp <- readOGR(".", "2008_w_fert_data")

# Here we link the projection information to the object which we loaded with 
# the "readShapePoly"-command and continue to work with this file that also has 
# an ID-variable defined.
proj4string(shape.shp) <- proj4string(shapeogr.shp)

# Each shapefile has an attribute table in which each observation 
# (e.g. country/region/individual) has an own row. Here we display the first
# six rows of the attribute table of the shapefile using the head-command:
head(shape.shp@data)

# Import data file
mydata <- read.table("DummyEast2.csv", sep=",", dec=".",header=T)


################################################################################
# 2) Merge shapefile and dataset                                               #
################################################################################ 

# Match by ID-column, which has to be similar in the shapefile and the datafile
# Define id of shapefile (used for matching).
id <- c(shape.shp$GS)
data.id <- mydata$ID

# Are attribute table of shape file and the data frame of the data file sorted
# in the same order?
identical(id,data.id)

# Has each id in both file a match with an ID in another file?
# This seems to be the case: When we sort the id columns of the two files, the 
# "identical"-command returns a TRUE.
identical(sort(id),sort(data.id))

# So we just need to bring the observations in the data frame (by row) in the 
# same order as the observations in the shapefile. For this we use the "match"-
# command.
order <- match(id, data.id)
mydata1 <- mydata[order,]

# Now we also obtain a TRUE without sorting the data frame.
identical(id,mydata1$ID)

# Link reordered data from data file to the attribute table of the shapefile
# using the spCbind-command (spatial column bind). This command returns an
# error, if the row names (ids) of the data file object and the shapefile object
# are not identical. 
row.names(mydata1) <- mydata1$ID
shape.shp <- spCbind(shape.shp, mydata1)

# Check, whether the data from the data file is now connected to the attribute 
# table of the shapefile.
names(shape.shp)

# Check, whether file you want to export has projection information included.
# We can also export shapefiles without projection information, but it is 
# better to have it included.
proj4string(shape.shp)

# Erase area information as it is causing error messages if we save the file
erase <- which(colnames(shape.shp@data)=="SHAPE_AREA")
shape.shp@data <- shape.shp@data[-erase] 

# Saving shapefile with projection information included.
writeOGR(shape.shp,".",'shapenew',driver="ESRI Shapefile",overwrite_layer=T)


################################################################################
# 3) Merge shapefile and dataset, if dataset includes NA                       #
################################################################################ 

# The procedure works in principle also if data is missing for some
# regions/observations in the shapefile.

mydata <- read.table("DummyEast2.csv", sep=",", dec=".",header=T)

# We erase ten rows with observations from the dataset.
mydata_na <- mydata[-c(10:19),]

# Now, our datafile is ten rows shorter.
length(mydata[,1])
length(mydata_na[,1])

# Same procedure as above
id <- c(shape.shp$GS)
data_na.id <- mydata_na$ID
order <- match(id, data_na.id)
mydata_na1 <- mydata_na[order,]

# Check, whether those rows that can be matched are now in the same order.
which_na <- which(is.na(mydata_na1$ID))
identical(id[-which_na],mydata_na1$ID[-which_na])

# For the row names of the reordered data attribute table, we have now to 
# take the ids from the shapefile attribute table, as the id column of the 
# data file contains now NAs.
row.names(mydata_na1) <- id
shape_na.shp <- spCbind(shape.shp, mydata_na1)

names(shape_na.shp)

# Check, whether file you want to export has projection information included.
proj4string(shape_na.shp)

# Saving shapefile with projection information included.
writeOGR(shape_na.shp,".",'shapenew_na',
         driver="ESRI Shapefile",overwrite_layer=T)
