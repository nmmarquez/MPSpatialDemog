################################################################################
#                                                                              #                          
# IDEM 156 - Spatial Demography Course 2018 - Technical File 9                 #                          
# Opening and writing dbf files                                                #                          
# Sebastian Klüsener, MPIDR                                                    #                          
#                                                                              #                          
################################################################################

# Erase all objects in workspace
rm(list=ls(all=TRUE))

# Load library
library(foreign)


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

# Load dbf file from shapefile
my.dbf.file <- read.dbf("2008_w_fert_data.dbf")

# Load data file
mydata <- read.table("DummyEast2.csv", sep=",", head=T)


################################################################################
# 2) Merge dbf-file and dataset                                                #
################################################################################ 

# Define ids
id <- c(my.dbf.file$GS)
data.id <- mydata$ID

# Id-vectors of shapefile and datafile identical?
identical(id,data.id)
# False

# Sorted id-vectors of shapefile and datafile identical?
identical(sort(id),sort(data.id))
# True

# Sort data.table according to id-ordering in dbf-file
o <- match(id, data.id)
mydata1 <- mydata[o,]

# Vectors identical
identical(id,mydata1$ID)
# True

# Use column-bind command to combine data to your dbf-file
my.dbf.file.upd <- cbind(my.dbf.file, mydata1)

head(my.dbf.file)
head(my.dbf.file.upd)

# Here you save to a new dbf-file which is not connected to the shapefile.
# 2008_w_fert_data
write.dbf(my.dbf.file.upd,"2008_w_fert_data_upd.dbf")

# Here you save to a new dbf-file which is connected to the shapefile.
# write.dbf(my.dbf.file.upd,"2008_w_fert_data.dbf")
