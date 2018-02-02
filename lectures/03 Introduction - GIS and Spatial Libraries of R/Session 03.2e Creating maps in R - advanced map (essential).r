################################################################################
#                                                                              #                          
# IDEM 156 - Spatial Demography Course 2018                                    #                          
# Introduction to Mapping Functions in R - Advanced Map                        #                          
# Dalia A. Conde and Fernando Colchero, University of Odense                   #
# with adaptations by Sebastian Klüsener, MPIDR                                #
#                                                                              #                          
################################################################################

# Erase all objects in workspace
rm(list=ls(all=TRUE))

# Load libraries
library(sp)
library(maptools)
library(rgdal)
library(RColorBrewer)
library(classInt)

# Set your working directory
# If you need to change the working directory, you can use the following code:
# Set drive
main <- c("N")
# Set path to session folders
path <- ":/IMPRSD/IDEM 156 Spatial"
# Define session folder
path2 <- "/02 Sessions/03 Introduction - GIS and Spatial Libraries of R"
# Set working directory
setwd(paste(main,path,path2,sep=""))


################################################################################
# 1. Data preparation                                                          #
################################################################################

# In this code it is assumed that you have already linked your data to
# the shapefile using the code of the tech3-helpfile.
# Load shapfile
shape.shp <- readShapePoly("2008_w_fert_data",IDvar="GS")
shapeogr.shp <- readOGR(".", "2008_w_fert_data")
proj4string(shape.shp) <- proj4string(shapeogr.shp)


################################################################################
#                                                                              #
# 2. Code for nice map with density function in right                          #
# upper corner and legend in the the right lower                               #
# corner. Intervals are defined by interval library                            #
#                                                                              #
################################################################################

################################################################################
# 2.1 Control panel for specifications                                         #
################################################################################

varofint <- shape.shp$TFR07

# Here you define the name of the variable how it should appear in the legend. 
# Wherever you put a "\n" in the next, a new  will be started. In this case,
# "Total Fertility Rate" will be printed in three lines.
var.name <- c("Total\nFertility\nRate")

# Here you choose the number of classes (colors) (bins-1)
nclassint <- 4

# Here you choose color scheme used
# To check available color schemes, run: display.brewer.all()
colorv <- c("YlOrRd")

# Here you define the border color for the regions
borcol <- c("grey25")

# Color in which potential NAs should be displayed
NAcolor <- c("grey")

# Here you can specify different categorizations (e.g. "quantile","equal",
# "jenks")
# For details on available schemes rund: ?classIntervals 
catgv <- c("quantile") 

# Here you can define the number of digits at which you want the numbers in
# legends to be rounded.
lar <- 1

# Here you can define the number of digits at which you want the axis numbers
# of the density curve plot to be rounded.
dar <- 1

# Here you define the size the text of the density curve description should have
daxsiz <- 1

# Here you define the size the text of the variable name above the legend should
# have.
vnsiz <- 1.2

# With tloc you move the variable name up and down in the legend dependent on 
# whether you change it to a positive or negative number
tloc <- 0

# Here you define the size the text of the legend is printed in map i
ltsiz <- 1 

# Put inctitle to TRUE if you want to add a title in the map window
#inctitle <- TRUE
inctitle <- FALSE

# In case you want to add a title to the plot in the map window, you would need
# to increase the upper margin of that window from 0 to a higher number. If 4 is 
# still not enough, you have to increase the number further
if (inctitle==FALSE) {martop <- 0}
if (inctitle==TRUE) {martop <- 4}

# To define text for the optional title
titletext <- c("Total Fertility Rate in Germany 2007")

# Here you define the size of the text of the title
titlesiz <- 1

# If you are printing numbers with decimal positions in the legend, R would 
# omit zeros at the end. In order to avoid this you can use the "sprintf"-
# function (?sprintf), but this function would need to be respecified dependent 
# on how many digits before and after the comma need to be displayed in order 
# to display the descriptions in the legend and the description of the density 
# curve correctly. Thus, the sprintf-function is added as an optional choice. 
# If you want to use it, you would need to set below use_sprintf  from FALSE 
# to TRUE. Afterwards you would need to specify the sprintf-function 
# specification in terms of digits that should be displayed in total and in 
# decimal positions based on the information you find under ?sprintf and the 
# number of digits you have to display in order to plot all values. In this
# example, the specification have already been adopted to the TFR-data. If
# you set use_sprintf to TRUE, they will be used. In the density, curve,
# on the other hand, 0s would not be omitted at the end
use_sprintf <- TRUE

# For the numbers displayed in the legend
spfleg <- "%3.1f" 

################################################################################
# 2.2 Code to produce the map                                                  #
################################################################################

# In the code in section 3.2 you might just need to change some layout-
# specifications in case you want to adapt it to your map. Apart from this, the 
# code should run automatically and use the specifications which were defined 
# above in the control panel. Here we derive the bins based on chosen 
# specifications.

# Export funtion
tiff(file="TFR_Germany.tif",width = 2400, height = 2400, res=300, 
     compression="lzw")
#png(file="TFR_Germany.png",width = 2400, height = 2400, res=300)
   catg <- classIntervals(varofint, nclassint, 
                       style = catgv)
   colpal <- brewer.pal(nclassint,colorv)
   color <- findColours(catg,colpal)
   color[is.na(varofint)] <- NAcolor
   bins <- catg$brks
   lb <- length(bins)
   mylayout <- layout(matrix(c(1,2),1,2), widths=c(1,0.25), heights=1)
   par(mar=c(0,0,martop,0))
   plot(shape.shp, col=color,border=borcol)
   if (inctitle==TRUE) {title(titletext, cex=titlesiz)}
   plot(c(0,1),c(0,1),col=NA,axes=F, xlab="", ylab="")
   points(rep(0.1,lb-1), seq(0.1,lb*0.025+0.1,length=lb-1), 
    pch=22, cex=2.5, col=borcol,bg=colpal)
   if (use_sprintf==FALSE) {text(rep(0.6,lb-1), seq(0.1,lb*0.025+0.1,
                                                    length=lb-1),
      paste(round(bins[-lb],lar),"-",round(bins[-1],lar)), cex=ltsiz)}
   if (use_sprintf==TRUE) {text(rep(0.6,lb-1), seq(0.1,lb*0.025+0.1,
                                                   length=lb-1),
      paste(sprintf(spfleg,round(bins[-lb],lar)),"-",
      sprintf(spfleg,round(bins[-1],lar))), cex=ltsiz)}
   text(0.5, lb*0.025+0.2+tloc, var.name, cex=vnsiz, font=2)
   op <- par(mar=c(1,1,1,1), fig=c(0.75,0.95,0.75,0.95), 
             new=TRUE, las=1)
   plot(density(varofint,na.rm=TRUE), xlab="", ylab="", axes=F, lwd=3, 
        main="")
   if (use_sprintf==FALSE) {axis(1, seq(round(min(varofint,na.rm=TRUE),dar),
       round(max(varofint,na.rm=TRUE),dar), length=2),cex.axis=daxsiz)}
   if (use_sprintf==TRUE) {axis(1, sprintf(spfleg,
                                   seq(round(min(varofint,na.rm=TRUE),dar),
                                       round(max(varofint,na.rm=TRUE),dar), 
                                       length=2)),cex.axis=daxsiz)}
   abline(v = bins, col=2, lwd=1.5, lty=3)
   par(op)
dev.off()

################################################################################
#                                                                              #
# 3. Simple code in case you want to produce fast a number of maps             #
#                                                                              #
################################################################################

# Code is adapted for four maps, but you can adapt it in case you want to 
# produce more or less maps

################################################################################
#                                                                              #
# 3.1 Control panel in which you can change the specifications such as which   #
#     variables you want to plot with which number of categories and which     #
#     categorization/color schemes. You can also change the position of the    #
#     legend, among other aspects                                              #  
#                                                                              #
################################################################################

# Here you define your data
var.data <- data.frame(shape.shp$TFR07,shape.shp$SNM07,
                       shape.shp$AVAGE07,shape.shp$NMC07)

# Here you provide the names of the variables
var.names <- c("Total Fertility Rate\n2007",
               "Share Non-marital Births\n2007",
               "Average Age at Birth\n2007",
               "Share Non-marital Conceptions\n2007")

# Here you choose the number of intervals in map i
nclassintv <- c(4,4,4,8)

# Here you choose the color schemes used in map i
colorv <- c("YlOrRd","YlOrRd","YlOrRd","Greens")

# Here you define the border color for the regions
borcol <- c("grey25")

# Color in which potential NAs should be displayed
NAcolor <- c("grey")

# Here you choose the categorization used in map i 
catgv <- c("quantile","quantile","quantile","equal") 

# Here you define the number of digits at which you want the numbers in
# legends to be rounded in map i
lar <- c(2,1,1,1)

# Here you define the size the text of the title is printed in map i
tsiz <- c(1,1,1,1) 

# Here you define the position of the legend in map i
# options are, e.g.: "bottomright","bottomleft",
#                    "topright","topleft"
posleg <- c("bottomright","bottomright","bottomright","bottomright")

# Here you define the size the text of the legend is printed in map i
ltsiz <- c(0.8,0.8,0.8,0.8) 


################################################################################
#                                                                              #
# 3.2 Plotting code to run the map. Here you would just need to make changes   #
#     to the first line in case you want to respecify your plot window.        #
#     mfrow=c(2,2) divides it in two rows and two columns, but also other      #
#     specifications are possible (e.g. mfrow=c(1,3) if you want to divide the #
#     plot window in one row and three columns. Apart from this the code       #
#     should run automatically based on the specifications set in 3.1          #
#                                                                              #
################################################################################

tiff(file="Maps_Germany.tif",width = 2800, height = 2400, res=300, 
     compression="lzw")
#png(file="Maps_Germany.png",width = 2800, height = 2400, res=300)
   par(mfrow=c(2,2), mar=rep(2,4))
   n <- length(var.data[,1])
   m <- length(var.data)
   for (i in 1:m) {
       varofint1 <- var.data[,i]
       catg <- classIntervals(varofint1, nclassintv[i], 
                              style=catgv[i])
       bins <- catg$brks
       lb <- length(bins)
       colpal <- brewer.pal(lb-1, colorv[i])
       colors <- findColours(catg,colpal)
       color[is.na(varofint)] <- NAcolor    
       plot(shape.shp, col=colors, border=borcol)
       title(var.names[i], cex=tsiz[i])
       legend(posleg[i],fill=colpal,
              legend=paste(round(bins[-lb],lar[i]),
                           "-", round(bins[-1],lar[i])),cex=ltsiz[i])
   }
dev.off()
