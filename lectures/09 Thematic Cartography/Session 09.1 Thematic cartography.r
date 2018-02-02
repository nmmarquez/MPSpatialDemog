################################################################################
#                                                                              #
# IDEM 156 - Spatial Demography Course 2018                                    #
# Thematic Cartography                                                         #
# Sebastian Klüsener, MPIDR                                                    #
#                                                                              #
################################################################################

# Erase all objects in workspace
rm(list=ls(all=TRUE))

# Open libraries
library(spdep)
library(RColorBrewer)
library(TeachingDemos)
library(rgdal)
library(maps)
library(classInt)
library(foreign)
library(maptools)
library(scales)

# Set your working directory
# If you need to change the working directory, you can use the following code:
# Set drive
main <- c("N")
# Set path to session folders
path <- ":/IMPRSD/IDEM 156 Spatial"
# Define session folder
path2 <- "/02 Sessions/09 Thematic Cartography"
# Set working directory
setwd(paste(main,path,path2,sep=""))


################################################################################
# 1)   Import and prepare data                                                 #
################################################################################

# Open shapefile with readShapePoly, which can read .prj-files
shapeogr.shp <- readOGR(".", "Ukraine.ll")

# Open shapefile with readShapePoly, which cannot read .prj-files
shape.shp <- readShapePoly('Ukraine.ll', IDvar="ID")

proj4string(shape.shp) <- proj4string(shapeogr.shp)

# Shapefile for background
background.ll <- readOGR(".", "Background.ll")

# Add water shapefile
water.shp <- readShapePoly('UKR_water_areas_dcw')
                          
# Generate some general information
id <-  paste(shape.shp$ID)
n <- length(id)
sq <- c(1:n)
sq_text <- paste(sq)
coords <- coordinates(shape.shp)



################################################################################
# 2)   Create and export a base map of Ukraine                                 #
# 2.1) Preparations                                                            #
################################################################################

# Maps should ideally contain a north arrow and a scalebar

# North arrow function from: http://www.jstatsoft.org/v19/c01/paper
northarrow <- function(loc,size,bearing=0,cols,cex=1,...) {
                      # checking arguments
                      if(missing(loc)) stop("loc is missing")
                      if(missing(size)) stop("size is missing")
                      # default colors are white and black
                      if(missing(cols)) cols <- rep(c("white","black"),8)
                      # calculating coordinates of polygons
                      radii <- rep(size/c(1,4,2,4),4)
                      x <- radii[(0:15)+1]*cos((0:15)*pi/8+bearing)+loc[1]
                      y <- radii[(0:15)+1]*sin((0:15)*pi/8+bearing)+loc[2]
                      # drawing polygons
                      for (i in 1:15) {
                           x1 <- c(x[i],x[i+1],loc[1])
                           y1 <- c(y[i],y[i+1],loc[2])
                      polygon(x1,y1,col=cols[i])
                      }
                      # drawing the last polygon
                      polygon(c(x[16],x[1],loc[1]),c(y[16],y[1],loc[2]),
                       col=cols[16])
                      # drawing letters
                      b <- c("E","N","W","S")
                     
                      for (i in 0:3) 
                       text((size+par("cxy")[1])*cos(bearing+i*pi/2)+loc[1],
                       (size+par("cxy")[2])*sin(bearing+i*pi/2)+loc[2],b[i+1],
                       cex=cex)
}

# Position of Kiev and other big cities 

# Get coordinates from Kiev from shapefile
Kiev <- coords[which(shape.shp$NAME=="KIEW (CITY)"),]

# Import coordinates from other cities (generated from the district-level 
# dataset.
bigcitext <- t(read.table("UKR_COOR_BIG_CITIES.csv", sep=","))
bigcit <- t(data.frame(bigcitext,Kiev))
lbc <-length(bigcit[,1])

# Adjust names of cities
row.names(bigcit)[2] <- "Donets'k"
row.names(bigcit)[3] <- "L'viv"
row.names(bigcit)[5] <- "Dnipropetrovs'k"

################################################################################
# 2.2) Plot base map                                                           #
################################################################################

# Plot base map of Ukraine
# With png you can specify width and height in pixel (standard unit), 
# inches (in), cm or milimeter (mm)
png(file="Ukraine.png",width=4800,height=3600, units="px", bg = "white",res=600)
#png(file="Ukraine.png",width=20.32,height=15.24, units="cm", bg = "white",res=300)
par(mfrow=c(1,1),mar=rep(2.5,4))
plot(shape.shp, col="white", bg="#F0FFFF", lwd=1, border="grey25", axes=TRUE, 
 cex.axis=0.8, col.axis="grey25")
plot(background.ll, col="grey75",lty=0,add=T)
plot(shape.shp, col="white", lwd=1.2, border="grey25", add=TRUE)
plot(shape.shp, col=alpha("darkolivegreen2",0.8), lty=0, add=TRUE)
plot(water.shp, add=TRUE, col="#F0FFFF",border="royalblue1", lwd=0.4)
title("Ukraine")
northarrow(c(39,51.8),0.7, cex=0.4)
map.scale(23,44.5, ratio=FALSE, cex=0.4)
for (i in 1:lbc) {
      symbols(bigcit[i,1], bigcit[i,2], circles=0.2, 
       inch=FALSE, add=TRUE, bg="red")
      symbols(bigcit[i,1], bigcit[i,2], circles=0.08, 
      inch=FALSE, add=TRUE, bg="black")
      text(bigcit[i,1]-0.1, bigcit[i,2]-0.4, row.names(bigcit)[i],cex=1)
}                           
dev.off()

# One way to avoid that text is not good to read, when it intersects with 
# borders. The following function is from the following source
# http://article.gmane.org/gmane.comp.lang.r.general/147787
shadowtext <- function(x, y=NULL, labels, col='white', bg='black',
	            theta= seq(pi/4, 2*pi, length.out=8), r=0.1, ... ) {
	
             	xy <- xy.coords(x,y)
            	xo <- r*strwidth('A')
            	yo <- r*strheight('A')

            	for (i in theta) {
		              text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, 
                   labels, col=bg, ... )
            	}
	            text(xy$x, xy$y, labels, col=col, ... )
}

# Base map of Ukraine with shadowtext
png(file="Ukraine_shadowtext.png",width=4800,height=3600, units="px", bg = "white",res=600)
par(mfrow=c(1,1),mar=rep(2.5,4))
plot(shape.shp, col="white", bg="#F0FFFF", lwd=1, border="grey25", axes=TRUE, 
     cex.axis=0.8, col.axis="grey25")
plot(background.ll, col="grey75",lty=0,add=T)
plot(shape.shp, col="white", lwd=1.2, border="grey25", add=TRUE)
plot(shape.shp, col=alpha("darkolivegreen2",0.8), lty=0, add=TRUE)
plot(water.shp, add=TRUE, col="#F0FFFF",border="royalblue1", lwd=0.4)
title("Ukraine")
northarrow(c(39,51.8),0.7, cex=0.4)
map.scale(23,44.5, ratio=FALSE, cex=0.4)
for (i in 1:lbc) {
  symbols(bigcit[i,1], bigcit[i,2], circles=0.2, 
          inch=FALSE, add=TRUE, bg="red")
  symbols(bigcit[i,1], bigcit[i,2], circles=0.08, 
          inch=FALSE, add=TRUE, bg="black")
  shadowtext(bigcit[i,1]-0.1, bigcit[i,2]-0.4, label=row.names(bigcit)[i],
             col="black",bg="white",r=0.1)
  shadowtext(bigcit[i,1]-0.1, bigcit[i,2]-0.4, label=row.names(bigcit)[i],
             col="black",bg=alpha("darkolivegreen2",0.8),r=0.1)
}   
dev.off()

                                                              
# Export map of Ukraine to pdf (width and height in inches)
pdf(file="Ukraine.pdf",bg = "transparent", width=8, height=6)
par(mfrow=c(1,1),mar=rep(2.5,4))
plot(shape.shp, col="white", bg="#F0FFFF", lwd=1, border="grey25", axes=TRUE, 
     cex.axis=0.8, col.axis="grey25")
plot(background.ll, col="grey75",lty=0,add=T)
plot(shape.shp, col="white", lwd=1.2, border="grey25", add=TRUE)
plot(shape.shp, col=alpha("darkolivegreen2",0.8), lty=0, add=TRUE)
plot(water.shp, add=TRUE, col="#F0FFFF",border="royalblue1", lwd=0.4)
title("Ukraine")
northarrow(c(39,51.8),0.7, cex=0.4)
map.scale(23,44.5, ratio=FALSE, cex=0.4)
for (i in 1:lbc) {
  symbols(bigcit[i,1], bigcit[i,2], circles=0.2, 
          inch=FALSE, add=TRUE, bg="red")
  symbols(bigcit[i,1], bigcit[i,2], circles=0.08, 
          inch=FALSE, add=TRUE, bg="black")
  shadowtext(bigcit[i,1]-0.1, bigcit[i,2]-0.4, label=row.names(bigcit)[i],
             col="black",bg="white",r=0.1)
  shadowtext(bigcit[i,1]-0.1, bigcit[i,2]-0.4, label=row.names(bigcit)[i],
             col="black",bg=alpha("darkolivegreen2",0.8),r=0.1)
} 
dev.off()


################################################################################
# 3)   Avoid overlapping of text on maps - Create and export a map             #
#      of Ukraine with the names of its regions                                #
# 3.1) First draft map with name of the regions                                #
################################################################################

names <- read.dbf("UKR_Names.dbf")
idlink <- paste(shape.shp$ID)
data.id <- paste(names$ID_ISO)
identical(sort(idlink),sort(data.id))
o <- match(idlink, data.id)
names1 <- names[o,]
identical(idlink,paste(names1$ID_ISO))
Name_Reg <- paste(names1$NAME_1)

# This draft map uses the shadowtext function, but still does not look good,
# as some region names overlap.
png(file="Ukraine_region_names.png",width=4800,height=3600, bg = "white",
    res=600)
par(mfrow=c(1,1),mar=rep(2.5,4))
plot(shape.shp, col="white", bg="#F0FFFF", lwd=1, border="grey25", axes=TRUE, 
     cex.axis=0.8, col.axis="grey25")
title("Ukrainian Regions")
plot(background.ll, col="grey75",lty=0,add=T)
plot(shape.shp, col="darkseagreen2", lwd=1, border="grey25", add=TRUE)
northarrow(c(39,51.8),0.7, cex=0.8)
map.scale(23,44.5, ratio=FALSE, cex=0.7)
shadowtext(coords[,1],coords[,2], label=paste(Name_Reg), cex=0.7, 
 col="black", bg="darkseagreen2",r=0.1)
dev.off()

################################################################################
# 3.2) For this map we try a function to avoid overlapping                     #
################################################################################

# Using pointLabel-function to create new coordinates, where text 
# is not overlapping
plot(shape.shp)
newcoord <- pointLabel(coords[,1],coords[,2], label=paste(Name_Reg), cex=0.7)

# Result does not look very good
png(file="Ukraine_region_names_point_Label.png",width=4800,height=3600,
    bg = "white",res=600)
par(mfrow=c(1,1),mar=rep(2.5,4))
plot(shape.shp, col="white", bg="#F0FFFF", lwd=1, border="grey25", axes=TRUE, 
     cex.axis=0.8, col.axis="grey25")
title("Ukrainian Regions")
plot(background.ll, col="grey75",lty=0,add=T)
plot(shape.shp, col="darkseagreen2", lwd=1, border="grey25", add=TRUE)
northarrow(c(39,51.8),0.7, cex=0.8)
map.scale(23,44.5, ratio=FALSE, cex=0.7)
shadowtext(newcoord$x,newcoord$y, label=paste(Name_Reg), cex=0.7, 
 col="black", bg="darkseagreen2",r=0.1)
dev.off()

################################################################################
# 3.3) in the following map we adjust manually the position                    #
################################################################################

# Manual adjustment of the coordinates
coordsadj <- coords
coordsadj[4,2] <- coordsadj[4,2]+0.1
coordsadj[10,2] <- coordsadj[10,2]-0.1
coordsadj[10,1] <- coordsadj[10,1]+0.5
coordsadj[16,2] <- coordsadj[16,2]-0.5
coordsadj[21,1] <- coordsadj[21,1]+0.3

# Map with adjusted region names
png(file="Ukraine_region_names_man_adj.png",width=4800,height=3600,bg="white",
    res=600)
par(mfrow=c(1,1),mar=rep(2.5,4))
plot(shape.shp, col="white", bg="#F0FFFF", lwd=1, border="grey25", axes=TRUE, 
     cex.axis=0.8, col.axis="grey25")
title("Ukrainian Regions")
plot(background.ll, col="grey75",lty=0,add=T)
plot(shape.shp, col="white", lwd=2, border="grey25", add=TRUE)
plot(shape.shp, col="darkseagreen2", lwd=1, border="grey25", add=TRUE)   
northarrow(c(39,51.8),0.7, cex=0.8)
map.scale(23,44.5, ratio=FALSE, cex=0.7)
shadowtext(coordsadj[,1],coordsadj[,2], label=paste(Name_Reg), cex=0.7, 
 col="black", bg="darkseagreen2",r=0.1)
dev.off()

# Save the non-overlapping coordinates
#coordadj.df <- data.frame(ID_1,NAME_1,coordsadj)
#write.table(file="Coordsadj.csv", sep=",", coordsadj.df, row.names=T)


################################################################################
# 4)   Map with pie charts -  Nationalities by Regions - Census 2001           #
# 4.1) Union regions of shapefiles (e.g. make out of two regions one region)   #
################################################################################

# We want to merge Kiev city with Kiev oblast, and Sevastopol city with 
# the Autonomous Republic of Crimea

# The regions to be merged are identified by a vector, in which regions to be 
# merged get the same value. We use the ids of the regions and our 
# identification vector is called idunion.

NAME <- paste(shape.shp$NAME)
# Identification of the ids of the regions to be merged
kievob <- id[which(NAME=="KIEW")]
kievci <- id[which(NAME=="KIEW (CITY)")]
sevast <- id[which(NAME=="SEVASTOPOL")]
crimea <- id[which(NAME=="AUTONOME REPUNLIK KRIM")]

# Kiev City gets the id of Kiev oblast, while Sevastopol get the ID of Crimea
idunion <- id
idunion[idunion==kievci] <- kievob
idunion[idunion==sevast] <- crimea
idunion

# gpclib has to be permitted to run, as its license only allows use for non-
# commercial purposes
#gpclibPermit() 

# Function to union regions in a polygon-file
reg25shape.shp <- unionSpatialPolygons(shape.shp, idunion)

# Compare input and output
par(mfrow=c(2,1)) 
plot(shape.shp)
title("Input shapefile with Kiev city and Sevastopol city")
plot(reg25shape.shp)
title("Output shape with cities\n included in Kiev oblast and AR Crimea")

# There is no data connected to this new shapefile
reg25shape.shp@data
# But the id-columns are available as names
names(reg25shape.shp)

# Connect data from old dataframe to new dataframe
# Detect in which row of the old spatial polygon dataframe the deleted regions
# were positioned. Here you would need to maked modifications according to which
# regions you have deleted
poskc <-which(NAME=="KIEW (CITY)")
possev <- which(NAME=="SEVASTOPOL")
delreg <- c(poskc, possev)

# Delete those regions from the dataframe and store it as a seperate dataframe  
datatomatch <- shape.shp@data[-delreg,]
identical(sort(paste(datatomatch$ID)),sort(names(reg25shape.shp)))
          
# Match dataframe to shapefile as a spatial polygon dataframe
reg25shape.shp <- SpatialPolygonsDataFrame(reg25shape.shp, datatomatch)

# Now we have again all the data availabe
names(reg25shape.shp)
id <- paste(reg25shape.shp$ID)


################################################################################
# 4.2) Import data                                                             #
################################################################################

# Upload data file
mydata10 <- read.table("UKR_NAT.csv", sep=",", header=T)
mydata20 <- read.table("UKR_POP.csv", sep=",", header=T)
mydata30 <- read.table("UKR_LONG.csv", sep=",", header=T)
mydata <- data.frame(mydata10,mydata20,mydata30)
                              
# Match data to shapefile
data.id <- paste(mydata$ID_ISO)

# Check, whether vector of data.ids and id of shapefile are of the same lenght
length(data.id) == length(id)

# Match by ID-column, which has to be similar in the shapefile and the datafile
o <- match(id,data.id)
mydata1 <- mydata[o,]
row.names(mydata1) <- id
reg25shape.shp <- spCbind(reg25shape.shp, mydata1)
names(reg25shape.shp)
coordsreg25 <- coordinates(reg25shape.shp)
n <- length(reg25shape.shp)

# In order to "attach" data from a shapefile attribute table, we have to call it
# with "@data"
attach(mydata1)


################################################################################
# 4.3) If we work with choropleth maps, we would need to produce a map for     #
#      each nationality                                                        #
################################################################################


# Choropleth maps of the six biggest nationalities
var.data <- data.frame(UKR01,RUS01,BLR01,MDA01,CRTAT01,BGR01)
var.names <- c("Ukrainian","Russian","Belarussian","Moldavian","Crim Tatarian",
 "Bulgarian")
# For each map another another RColor-Brewer color scheme is used
colorsc <- c("Oranges","Blues","Reds","Purples","YlOrBr","Greens")

png(file="Ukraine_nationalities.png",width=4800,height=3600,bg = "white",
    res=600)
par(mfrow=c(2,3),mar=rep(1.2,4))
# This code produces m maps with self-defined categorization (useful, if you 
# want to contrast several maps)
m <- length(var.data)
var.color <- matrix(nrow=n,ncol=m, data=rep(0,n*m))
for (i in 1:m) {
    varofint <- var.data[,i]
    bins <- c(0,0.5,5,10,20,50,70,90,100)
    lb <- length(bins)
    colpal <- c("white",brewer.pal(lb-2, colorsc[i]))
    var.color <- colpal[findInterval(varofint, bins, rightmost.closed=T)]
    plot(reg25shape.shp, col="white", border="grey25",bg="#F0FFFF",)
    title(var.names[[i]], cex=1)
    plot(background.ll, col="grey75",lty=0,add=T)
    plot(reg25shape.shp, col=var.color, lwd=1, border="grey25", axes=TRUE, 
         cex.axis=0.8, col.axis="grey25",add=T)
    legend("bottomleft",fill=colpal,legend=paste(round(bins[-length(bins)],2),
    "-", round(bins[-1],2)),cex=0.5, bg="white")
}
dev.off()


################################################################################
# 4.4) A more effective way might be to work with pictograms. In the           #
#      following example we will use pie charts                                #
################################################################################

# Dataframe with data by nationality
var.data <- data.frame(UKR01,RUS01,BLR01,MDA01,CRTAT01,BGR01,HUN01,ROU01,POL01,
 JEW01,GRC01,TAT01,ROMA01,SVK01,GAGAU01,OTH01)

# The pie requires positive numbers. Therefore, we replace all NAs by a very 
# small positive number.
var.data[is.na(var.data)] <- 0.00000000000000000000001

# Vector of piecolors
piecolor <- c("yellow","red","orange","mediumpurple4","tan4","lightgreen",
 "forestgreen","plum3","pink","cornflowerblue","tan3","navy",
 "lightblue","coral4","darkgoldenrod","grey")

# Plot Map
png(file="Ukraine_nationalities_circles.png",width=4800,height=3600,
    bg="white",res=600)
par(mfrow=c(1,1),mar=rep(2.5,4))
plot(reg25shape.shp, col="white", bg="#F0FFFF", lwd=0.8, border="grey25", 
     axes=TRUE, cex.axis=0.8, col.axis="grey25")
title("Nationalities in Ukrainian Regions - Census 2001")
plot(background.ll, col="snow3",lty=0,add=T)
plot(reg25shape.shp, col="snow1", lwd=0.8, border="grey25", axes=TRUE, 
     cex.axis=0.8, col.axis="grey25",add=T)
northarrow(c(39,51.8),0.7, cex=0.8)
map.scale(23,44.5, ratio=FALSE, cex=0.7)
for (i in 1:n) {
    subplot(pie(c(100), labels="",radius = 0.41, lwd=2, border="grey25", 
    col="black"), coordsreg25[i,1], coordsreg25[i,2])
    subplot(pie(as.numeric(var.data[i,]), labels="",col=piecolor, 
    clockwise=TRUE, init.angle=90, radius = 0.4, lty=0), coordsreg25[i,1], 
    coordsreg25[i,2])
}                                                                          
dev.off()

################################################################################
# 4.5) This map displays the same information, but we let the pie chart vary   #
#      by population size of the region                                        #
################################################################################                                                                                                

# Calculate radius of proportional circles
TOP2001 <- POP2001/5000
rad <- sqrt(TOP2001/pi)    

# Same as in 4.4     
var.data <- data.frame(UKR01,RUS01,BLR01,MDA01,CRTAT01,BGR01,HUN01,ROU01,POL01,
 JEW01,GRC01,TAT01,ROMA01,SVK01,GAGAU01,OTH01)
var.data[is.na(var.data)] <- 0.00000000000000000000001
piecolor <- c("yellow","red","orange","mediumpurple4","tan4","lightgreen",
 "forestgreen","plum3","pink","cornflowerblue","tan3","navy",
 "lightblue","coral4","darkgoldenrod","grey")

# Plot map
png(file="Ukraine_nationalities_circles_popsize.png",width=4800,height=3600,
    bg="white",res=600)
par(mfrow=c(1,1),mar=rep(2.5,4))
plot(reg25shape.shp, col="grey25", bg="#F0FFFF", lwd=0.8, border="grey25", 
     axes=TRUE, cex.axis=0.8, col.axis="grey25")
title("Nationalities in Ukrainian Regions - Census 2001")
plot(background.ll, col="snow3",lty=0,add=T)
plot(reg25shape.shp, col="snow1", lwd=0.8, border="grey25", axes=TRUE, 
     cex.axis=0.8, col.axis="grey25",add=T)
northarrow(c(39,51.8),0.7, cex=0.8)
map.scale(23,44.5, ratio=FALSE, cex=0.7)
for (i in 1:n) {
    subplot(pie(c(100), labels="",radius = rad[i]+0.01, lwd=2, border="grey25", 
    col="black"), coordsreg25[i,1], coordsreg25[i,2])
    subplot(pie(as.numeric(var.data[i,]), labels="",col=piecolor, 
    clockwise=TRUE, init.angle=90, radius=rad[i], lty=0), coordsreg25[i,1], 
    coordsreg25[i,2])
}                                                                          
dev.off()

################################################################################
# 4.6) Here we use the layout function to produce the map with legend          #
################################################################################                                                                                                

TPOP2001 <- POP2001/5000
LegPOP2001 <- c(5000,2500,1000,500)
LegTPOP2001 <- LegPOP2001/5000 
rad <- sqrt(TPOP2001/pi)
legrad <- sqrt(LegTPOP2001/pi)    
coordsreg25adjnat <- coordsreg25 

# Plot map
png(file="Ukraine_nationalities_circles_popsize_legend.png",width=4800,
    height=3600,bg="white",res=600)
layout(matrix(c(1,2,3,3),2,2), widths=c(1,0.25), heights=c(1,0.25))
plot(reg25shape.shp, col="grey25", bg="#F0FFFF", lwd=0.8, border="grey25", 
     axes=TRUE, cex.axis=0.8, col.axis="grey25")
title("Nationalities in Ukrainian Regions - Census 2001")
plot(background.ll, col="snow3",lty=0,add=T)
plot(reg25shape.shp, col="snow1", lwd=0.8, border="grey25", axes=TRUE, 
     cex.axis=0.8, col.axis="grey25",add=T)
northarrow(c(39,51.3),0.6, cex=0.7)
map.scale(23,45, ratio=FALSE, cex=0.7)
for (i in 1:n) {
    subplot(pie(c(100), labels="",radius = rad[i], lwd=2, border="grey25", 
    col="grey25"), coordsreg25adjnat[i,1], coordsreg25adjnat[i,2])
    subplot(pie(as.numeric(var.data[i,]), labels="",col=piecolor, 
    clockwise=TRUE, init.angle=90, radius=rad[i]-0.01, lty=0), 
    coordsreg25adjnat[i,1], coordsreg25adjnat[i,2])
}
plot(c(1,0),c(1,0),col=NA,axes=F, xlab="", ylab="")
# Add legend
par(mar=c(rep(0,4)))
posx <- sort(rep(c(0,0.22,0.44,0.66),4))
posy <- rep(c(1,0.7,0.4,0.1),4)
# Nationality
points(posx,posy, bg=piecolor, pch=21, col="grey25",lwd=1.2)
text(posx+0.02,posy,adj=0, labels=c("Ukrainian","Russian","Belarussian",
 "Moldavian","Crim. Tatarian","Bulgarian", "Hungarian","Romanian","Polish",
 "Jewish","Greek","Tatarian","Roma","Slovakian","Gagausian", "Other"), cex=0.8)
text(0.95, 1,c("Population Size"), cex=0.8)
posxp <- rep(c(0.95),4)
posyp <- rep(c(0.4),4)
subplot(pie(c(100), radius=legrad[1], labels="5 Mio.", 
 lwd=2, col="white", border="grey25", lty=1, cex=0.8),0.98,0.4)
dev.off() 
 
 
################################################################################
# 5)   Choropleth map using the shadowtext function                            #
#      Male Life Expectancy in Ukrainian Regions 2003/2004                     #
################################################################################

# Plot map
var.datale <- LEML0304
var.name <- c("Life\n Expectancy\n Male\n2003/ 2004")
coordsadj25 <- coordinates(reg25shape.shp)
coordsadj25n <- coordsadj25
coordsadj25n[22,2] <- coordsadj25n[22,2]+0.3

# Here you choose your data and the number of classes
nclassint <- 6

# Here you can specify different categorizations
varofint <- var.datale
cat1 <- classIntervals(varofint, nclassint, style = "quantile")
cat2 <- classIntervals(varofint, nclassint, style = "equal")
cat3 <- classIntervals(varofint, nclassint,style = "jenks")
colpal <- brewer.pal(nclassint,"RdYlGn")

# Here you choose the category, in which the map should be printed. These
# you have defined above
categ <- cat1

# Code for map
png(file="Ukraine_life_expectancy_male.png",width=4800,height=3600,
    bg="white",res=600)
color <- findColours(categ,colpal)
bins <- categ$brks
lb <- length(bins)
layout(matrix(c(1,2),1,2), widths=c(1,0.25), heights=c(1))
par(mar=rep(2.5,4))
plot(reg25shape.shp, col="grey25", bg="#F0FFFF", lwd=0.8, border="grey25", 
     axes=TRUE, cex.axis=0.8, col.axis="grey25")
title("Life Expectancy Males in Ukrainian Regions - 2003/2004")
plot(background.ll, col="grey75",lty=0,add=T)
plot(reg25shape.shp, col=color, lwd=0.8, border="grey25", axes=TRUE, 
     cex.axis=0.8, col.axis="grey25",add=T)
shadowtext(coordsadj25n[,1],coordsadj25n[,2], 
           label=paste(reg25shape.shp$REGIONS), cex=0.6, col="black", 
           bg=color,r=0.1)
northarrow(c(39,51.8),0.7, cex=0.8)
map.scale(23,44.5, ratio=FALSE, cex=0.7) 
par(mar=rep(0,4))
plot(c(0,1),c(0,1),col=NA,axes=F, xlab="", ylab="")
points(rep(0.1,lb-1), seq(0.1,lb*0.025+0.1,length=lb-1), pch=15, cex=2.5, col=1)
points(rep(0.1,lb-1), seq(0.1,lb*0.025+0.1,length=lb-1), pch=15, cex=2,
 col=colpal)
text(rep(0.6,lb-1), seq(0.1,lb*0.025+0.1,length=lb-1),
 paste(round(bins[-lb],2),"-",round(bins[-1],2)), cex=0.9)
text(0.5, lb*0.025+0.2, var.name, cex=0.9, font=2)

op <- par(mar=c(0.5,0.5,0.5,0.5), fig=c(0.75,0.95,0.75,0.95), new = TRUE, las=1)
plot(density(varofint), xlab="", ylab="", axes=F, lwd=3, main="")
axis(1, seq(min(varofint), max(varofint), length=2),cex.axis=0.9)
abline(v = bins, col=2, lwd=1.5, lty=3)
par(op)
dev.off() 

################################################################################
# 6)   Barplot map showing deviation from national average                     #
#      Life Expectancy (Total, Females, Males)                                 #
################################################################################

# Manual Adjustment, so that graphs do not overlap
coordsreg25adjnat <- coordsadj25 
coordsreg25adjnat[1,1] <- coordsreg25adjnat[1,1]+0.2
coordsreg25adjnat[2,1] <- coordsreg25adjnat[2,1]-0.2
coordsreg25adjnat[7,] <- coordsreg25adjnat[7,]+0.4
coordsreg25adjnat[8,] <- coordsreg25adjnat[8,]-0.3
coordsreg25adjnat[10,1] <- coordsreg25adjnat[10,1]+0.5
coordsreg25adjnat[10,2] <- coordsreg25adjnat[10,2]-0.2
coordsreg25adjnat[11,1] <- coordsreg25adjnat[11,1]+0.3
coordsreg25adjnat[14,1] <- coordsreg25adjnat[14,1]-0.4
coordsreg25adjnat[22,1] <- coordsreg25adjnat[22,1]+0.8
coordsreg25adjnat[25,1] <- coordsreg25adjnat[25,1]+0.4
coordsreg25adjnat[25,2] <- coordsreg25adjnat[25,2]-0.3


# Plot map
png(file="Ukraine_life_expectancy_deviations.png",width=4800,height=3600,
    bg="white",res=600)
var.datal <- data.frame(DULT0304,DULF0304,DULM0304)
layout(matrix(c(1,2,3,3),2,2), widths=c(1,0.25), heights=c(1,0.25))
par(mar=c(rep(2.5,4)))
plot(reg25shape.shp, col="grey25", bg="#F0FFFF", lwd=0.8, border="grey25", 
     axes=TRUE, cex.axis=0.8, col.axis="grey25")
title(c("Life Expectancy in Ukrainian Regions - 2003/2004"))
plot(background.ll, col="snow3",lty=0,add=T)
plot(reg25shape.shp, col="snow1", lwd=0.8, border="grey75", axes=TRUE, 
     cex.axis=0.8, col.axis="grey25",add=T)
northarrow(c(39,51.3),0.7, cex=0.8)
map.scale(23,44.5, ratio=FALSE, cex=0.7) 

# Subplot with bars
for (i in 1:n) { 
subplot(barplot(as.numeric(var.datal[i,]),ylim=c(-3.5,3.5), 
 col=c("white","orange","cornflowerblue"),border="grey50"), 
 x=coordsreg25adjnat[i,1], 
 y=coordsreg25adjnat[i,2], size=c(.3,.3))
}
# Legend
plot(c(0,1),c(0,1),col=NA,axes=F, xlab="", ylab="")
text(0.25,0.95, c("Deviation from national average in years"), cex=0.9, font=2)                                 
subplot(barplot(c(3,-2,1),ylim=c(-3.5,3.5), 
 col=c("white","orange","cornflowerblue"),border="grey50"),
 x=0.05, y=0.6, size=c(.6,.6))
text(c(0.13),c(0.6), adj=0, 
 labels=c("White - Total \nOrange - Female\nBlue - Male"))  
dev.off()
