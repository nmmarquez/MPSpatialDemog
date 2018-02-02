################################################################################
#                                                                              #                         
# IDEM 156 - Spatial Demography Course 2018                                    #                          
# Calculating and Mapping Local Indicators of Spatial Autocorrelation          #                                                                                                                                 #
# Sebastian Klüsener, MPIDR                                                    #                         
#                                                                              #                         
################################################################################

# Erase all objects in workspace
rm(list=ls(all=TRUE))

# Load libraries
library(spdep)
library(RColorBrewer)
library(maptools)
library(rgdal)


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

shape.shp <- readShapePoly('German_Districts_0406', IDvar="DISTID")
shapeogr.shp <- readOGR(".", "German_Districts_0406")
proj4string(shape.shp) <- proj4string(shapeogr.shp)
# Reprojected version of the shapefile which has a longlat projection
shape.ll <- spTransform(shape.shp, CRS("+proj=longlat +datum=WGS84"))

plot(shape.shp)                   
names(shape.shp)

# Define ID
id <- shape.shp$DISTID

# Number of observations
n <- length(id)

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

################################################################################            
# 2) Calculate LISA                                                            #
################################################################################

################################################################################
# 2.1) In this section we set definitions for the rest of the code             #
################################################################################

# Which variables are available in the shape.file attribute table
names(shape.shp)

# Define variable of interest
varofint <- AVAGE05

# Define name of the variable of interest 
varofint.name <- c("Average Age at Birth 2005")

# Define neighborhood matrix (type nb) 
# (choice in this case nb.FOQ.cor, d50km, nb.5NN)
nb <- nb.FOQ.cor

# Define weight style (W=row-standardized)
ws <- c("W")

# Define significance level for the cluster maps 
# (0.0001,0.001,0.01,0.05)
significance <- 0.05

# p-adjustment method (can be "none", "bonferroni", "holm",
# "hochberg","hommel","fdr")
p.ad.meth <- c("none")

# Should the cluster map show only regions belonging to significent clusters, 
# or all regions
plot.only.significant <- TRUE


################################################################################
# 2.2) In this section no changes to the code are necessary                    #
################################################################################

# Transfer nb-object in listwise object
listw <- nb2listw(nb, style=ws,zero.policy=T)

# Create lagged values of variable of interest
varlag <- lag.listw(listw, varofint)

# Calculate Lisa Test
lisa.outcomes <- localmoran(varofint,listw, alternative="two.sided", p.adjust.method=p.ad.meth)

# Get significance level
vec <- c(1:n)
vec <- ifelse(lisa.outcomes[,5] < significance, 1,0)

# Result object
lisa.outcomes
    

################################################################################            
# 3) Plot LISA-maps                                                            #
# 3.1) Preparations (no changes necessary)                                     #
################################################################################

# Calculate Mean of Variable of interest and lagged value of variable of 
# interest
m.varofint <- mean(varofint)
m.varlag <- mean(varlag)

# Derive sector
sec <- c(1:n)
for (i in 1:n) {
    if (varofint[[i]]>=m.varofint & varlag[[i]]>=m.varlag) sec[i] <- 1
    if (varofint[[i]]<m.varofint & varlag[[i]]<m.varlag) sec[i] <- 2
    if (varofint[[i]]<m.varofint & varlag[[i]]>=m.varlag) sec[i] <- 3
    if (varofint[[i]]>=m.varofint & varlag[[i]]<m.varlag) sec[i] <- 4
}

# Define colors for sectors
sec.all <- sec
colors1 <- c(1:n)
for (i in 1:n) {
    if (sec.all[i]==1) colors1[i] <- "brown2"
    if (sec.all[i]==2) colors1[i] <- "royalblue3"
    if (sec.all[i]==3) colors1[i] <- "lightblue"
    if (sec.all[i]==4) colors1[i] <- "pink"
    if (sec.all[i]==0) colors1[i] <- "white"
}

# Mark all non-significant regions white
loc.m.data <- sec*vec
colors2 <- colors1
for (i in 1:n) {
    if (loc.m.data[i]==0) colors2[i] <- "white"
}


################################################################################            
# 3.2) Plotting maps - here you might want to change the specifications of     #
#      some specifcations such as the format of the png-file which you want to #
#      export or the size of the text.                                         #                                     
################################################################################

png(file="LISA_Maps.png",width = 3200, height = 1500, res=300)
# Cluster map
par(mfrow=c(1,2),mar=c(1, 1, 3, 1))
if (plot.only.significant==TRUE) 
   {plot(shape.shp, col=colors2, border="grey25",lwd=0.7)} else 
   {plot(shape.shp, col=colors1, border="grey25",lwd=0.7)}
   legend("bottomright",fill=c("brown2","royalblue3","lightblue",
                               "pink","white"),
          legend=c("High-High","Low-Low","Low-High","High-Low"),
          border = "grey25", cex=0.8, bg="white", 
          box.col="white")
    title(paste("Significant Clusters\n",varofint.name))

# Significance map
sig.data <- data.frame(lisa.outcomes[,5])
m <- length(sig.data)
colors <- matrix(nrow=n,ncol=m, data=rep(0,n*m))
for (i in 1:m) {
    signific <- sig.data[,i]
    lb <- 5
    bins <- c(0,0.0001,0.001,0.01,0.05,1)
    colpal <- c(rev(brewer.pal(lb+1, "YlGn")[-c(1:2)]),"white")
    colors[,i] <- colpal[findInterval(signific, bins, rightmost.closed=T)]
    colors[vec==0,i] <- c("white")
    plot(shape.shp, col=colors[,i], border="grey25", lwd=0.7)
    title("Significance Level")    
    colpalad <- colpal
    colpalad[c(which(bins==significance):length(colpalad))] <- c("white")
    binsp <- c("0","0.0001","0.001","0.01","0.05",1)
    legend("bottomright",fill=colpalad,
           legend=paste(binsp[-length(bins)],"-",binsp[-1]), border="grey25", 
           cex=0.8, bg="white", box.col="white")
}
dev.off()
