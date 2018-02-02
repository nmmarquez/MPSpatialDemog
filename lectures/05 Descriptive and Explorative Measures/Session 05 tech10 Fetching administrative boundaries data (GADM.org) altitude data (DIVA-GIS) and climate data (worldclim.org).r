################################################################################
#                                                                              #                          
# IDEM 156 - Spatial Demography Course 2018 - Technical File 10                #                          
# Importing shapefiles from GADM.org, DIVA-GIS and World Climate               #                          
# Sebastian Kl?sener, MPIDR                                                    #                          
#                                                                              #                          
################################################################################


# Erase all objects in memory
rm(list = ls(all = TRUE))

#install.packages(c("stringr")),

# Load libraries
library(maptools)
library(raster)
library(rgdal)
library(stringr)
library(RColorBrewer)

################################################################################
# 1) Fetch data from GADM.org                                                  #
################################################################################ 

# Download administrative boundaries data for Senegal from GADM.org
# Level defines the level of the administrative hierarchy. Download specifies
# whether also a local copy should be stored in your working directory
# Level 1, Country is chosen via ISO3 code
shape_sen1.shp <- getData("GADM", country="SEN", level=1,download=TRUE)
head(shape_sen1.shp@data)

# Level 2
shape_sen2.shp <- getData("GADM", country=c("SEN"), level=2,download=TRUE)
head(shape_sen2.shp@data)

# Plot shapefiles of first and second level administrative hierarchy
par(mfrow=c(1,2))
plot(shape_sen1.shp,col="navajowhite1", border="darkgoldenrod3")
plot(shape_sen2.shp,col="navajowhite1",border="darkgoldenrod3")


# Combine shapefiles for two countries (GADMTools library offers to download
# several countries at once, but I was not able to get it going on my 
# computer)
# Mali
shape_mli.shp <- getData("GADM", country="MLI", level=1)
# Gambia
shape_gam.shp <- getData("GADM", country="GMB", level=1)

# Combine the two shapefiles with the rbind-command
sen_mli_gam <- rbind(shape_sen1.shp, shape_mli.shp,shape_gam.shp)

par(mfrow=c(1,1))
plot(sen_mli_gam,col="navajowhite1", border="darkgoldenrod3")

# To combine shapefiles of different hierarchies, we first need to adjust the
# number of columns in the attribute tables so that they match
sen_mal <- rbind(shape_sen2.shp, shape_mli.shp)

# Level 2 file of Senegal
sen_l2 <- shape_sen2.shp
# Level 1 file of Mali
mli_l1 <- shape_mli.shp

# The unique ids of the level 2 file are in the column ID_2
# We are creating an id that includes the ISO code, and erase at the same
# time all other columns of the shapefile
sen_l2@data <- data.frame(paste(sen_l2$ISO,sen_l2$ID_2,sep="_"))

# The unique ideas of the level 1 file are in the column ID_1
mli_l1@data <- data.frame(paste(mli_l1$ISO,mli_l1$ID_1,sep="_"))

# Harmonize the name of the id columns
colnames(sen_l2@data) <- c("id")
colnames(mli_l1@data) <- c("id")

# Now we can combine the files from different hierarchies
sen2_mli1 <- rbind(sen_l2, mli_l1)

plot(sen2_mli1,col="navajowhite1", border="darkgoldenrod3")


################################################################################
# 2) Fetch elevation data from DIVA-GIS                                        #
################################################################################ 

# Download elevation data from DIVA-GIS
# mask=FALSE allows us to include or exlude elevation data of territories 
# adjacent to the country for which we are exporting the data
altitude_data_mf <- getData('alt', country='SEN', mask=FALSE,download=T)
altitude_data_mt <- getData('alt', country='SEN', mask=TRUE,download=T)

# Plot of the elevation data with administrative boundaries
par(mfrow=c(1,2))
# colNA defines the background, we use 10 terrain colors (bins are in this case
# derived manually, but we could also derive them automatically - note how 
# massively the range and categorisation varies between the two maps). 
plot(altitude_data_mf,col=terrain.colors(10),colNA="#f0f8ff")
plot(shape_sen1.shp,add=T,border="grey25")
plot(altitude_data_mt,col=terrain.colors(10))
plot(shape_sen1.shp,add=T,border="grey25")

par(mfrow=c(1,1))
altitude_data_mf <- getData('alt', country='SEN', mask=FALSE)
plot(altitude_data_mf,col=terrain.colors(10),colNA="#f0f8ff",alpha=0.9)
plot(shape_sen1.shp,add=T,border="grey25")


################################################################################
# 3) Fetch climate data from worldclimate.org                                  #
################################################################################ 

# Download climate data from worldclimate.org
clim <- getData('worldclim', var='bio', res=0.5,lon=-15,lat=15,download=TRUE)
clim1 <- crop(clim,altitude_data_mf)

# Deriving some color vectors using RColorBrewer
greys <- brewer.pal(9,"Greys")
blues <- brewer.pal(9,"Blues")

# Export map with elevation, precipitation and administrative boundaries
png(file="Map_Senegal.png",width = 2800, height = 2200, res=300)
   # Alpha allows to modify the opacity of the raster data which we would like
   # to plot. This enables us to plot several layers above each other.
   plot(altitude_data_mf,col=greys,colNA="grey99",alpha=0.9,legend=F,
        main="Senegal - Elevation and Average Annual Precipitation")
   plot(clim1$bio12_25,col=blues,alpha=0.4,add=T)
   plot(shape_sen2.shp,add=T,border="salmon3",lwd=0.7)
dev.off()

# As data on these web databases might get updated, it might be worthwhile to 
# save a local copy to make sure that your analysis results are are also in the
# future reproducible

# Adding a date stamp
date <- paste(str_sub(format(Sys.time(), "%Y"),start=3,end=4),
                     format(Sys.time(), "%m"),format(Sys.time(), "%d"),sep="")
writeOGR(shape_sen2.shp,".",paste(date,"senegal_adm2",sep="_"),
         driver="ESRI Shapefile",overwrite_layer=T)


