################################################################################
#                                                                              #                          
# IDEM 156 - Spatial Demography Course 2018                                    #                          
# Introduction to Mapping Functions in R - Advanced Map                        #                          
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
library(classInt)

# Set your working directory
# If you need to change the working directory, you can use the following code:

################################################################################
# 1. Data preparation                                                          #
################################################################################

# In this code it is assumed that you have already linked your data to
# the shapefile using the code of the tech3-helpfile.

# The following read-command to load shapefiles into R allows to define an ID-
# variable according to which the shapes and the attribute table of the 
# shapefile will be sorted. 
shape.shp <- readShapePoly("2008_w_fert_data",IDvar="GS")

# The following read-command to load shapefiles into R automatically uploads the
# projection information which would need to be included manually in the
# readShapePoly-function. We open now the same shapefile just to then transfer 
# the projection-information to the object shape.shp which we uploaded with the
# readShapePoly-function. We will then continue to work with the shape.shp-
# object as it offers the advantage that the attribute table is sorted according
# to our defined ID-variable, a property which might be of help at later stages
# of the coding.
shapeogr.shp <- readOGR(".", "2008_w_fert_data")
proj4string(shape.shp) <- proj4string(shapeogr.shp)

# Test plot of the shapefile which you loaded.
plot(shape.shp)

# To obtain a summary of the shapefile characteristics
#(Coordinates of extent, projection, Attributes.)
summary(shape.shp)

# Head of data attribute table
head(shape.shp@data)

# Names of the spatial attributes(i.e the column names of 
# the *.dbf file)
names(shape.shp)

# NOTE: Look at the attributes in Excel by opening the 
# *.dbf file in Windows 2008_w_fert_data.dbf
# If the extensions of the files are hidden you can look at
# them under the file type

# Work with one attribute: TFR07 (Total Fertility Rate) in 
# 2007. A simple way to explore the data is to look at the
# data range range(shape.shp$TFR07). 

# To obtain the minimum and maximum of the TFR07
range(shape.shp$TFR07) 

# Plot density function
plot(density(shape.shp$TFR07), lwd=7, main="Density Curve")


################################################################################
# 2. Categorization to color map                                               #
################################################################################

# If we want to color the fertility data by bins of equal 
# distance then we first create the bins. "Bins" is a 
# vector object that we will use as a reference to assign a
# range of colors to the polygons with specific TFR.

# Define number of bins to be used
lb <- 4

# Create bins of equal distance to make a histogram
bins <- seq(min(shape.shp$TFR07), max(shape.shp$TFR07), length= lb)  

# An alternative code which uses quantiles
# bins = quantile(shape.shp$TFR07,seq(0,1,length=lb))

bins

# Par allows to modify plot window. Here we used a plot 
# window with one column and one row.
par(mfrow=c(1,1)) 

# Plot histogram, breaks indicates the number of sections
# in the histogram
hist(shape.shp$TFR07, breaks=20)     

# Insert red vertical lines to show bins
abline(v=bins, col="red")

################################################################################
# Assign color to each bin that cuts the TFR data in equal                     #
# intervals. There are many color libraries in R. In our                       #
# course we use color brewer, which creates colors that are                    # 
# very similar independent on whether they are displayed on                    #
# a computer screen or as a print out.                                         #
################################################################################

# Choose a color: in this case we will use the first color 
# palette
display.brewer.all()

# When you create an object of colors named "colpal" it 
# assigns "length(bins)-1" colors to the bins. Above we 
# defined four bins, which implies we will define three 
# colors to cover the three intervals between the bins.
# We use the -1 to indicate that we only want 3 colors from
# the #"YlOrRd" pallete.
colpal <- brewer.pal(length(bins)-1,"YlOrRd")

# Note when you type colpal it gives you the codes of the 3
# colors that will be displayed
colpal

# Create your color palette with the number of bins-1
# and color each bin according to the interval to which it
# belongs in the distribution  Use the findInterval: this 
# function reads in each of the data in the initial vector 
# (in this case shape.shp$TFR07 and finds to which bin it
# corresponds (1-2, 2-3 or 3-4) and gives you the position 
# of the bin. Thus, the first bin 1-2 will assign the  
# first color of the "YlOrRd" from the brewer palette.
colpal <- brewer.pal(lb-1, "YlOrRd")
colors <- colpal[findInterval(shape.shp$TFR07, bins, 
 rightmost.closed=T)]

# To show what is happening here exactly
colpal
bins

# First entry in data attribute table of shapefile
shape.shp@data[1,]
interval <- findInterval(shape.shp$TFR07[1], bins, 
 rightmost.closed=T)
colpal[interval]

# Set the plotting space to have two parallel plotting 
# windows (two columns)
par(mfrow=c(1,2))

# On the left side we create a graph which shows the density and in the 
# background the colors. Use a for-loop to fill to assign the colors to the 
# bins.
plot(density(shape.shp$TFR07),lwd=3, 
 main="Total Fertility Rate",xlab="", ylab="")
for(i in 1:(length(bins)-1)){ 
   polygon(c(bins[i],bins[i+1],bins[i+1],bins[i]),
   c(0,0,4,4), col=colpal[i],border=NA)
}

# Draw a line that shows the distribution using the density() function 
# We do first a white line thicker lwd=6
lines(density(shape.shp$TFR07),lwd=6,col='white')
# Then we overwrite the thicker line with a thinner line in black
lines(density(shape.shp$TFR07),lwd=3)

#Make a legend from scratch
plot(c(0,1),c(0,1),col=NA,axes=F, xlab="", ylab="")
points(rep(0.1,length(bins)-1),seq(0.1,0.8,
 length=length(bins)-1),pch=15,col=colpal, cex=8)
text(0.1,0.95,"Colpal",cex=1)
text(rep(0.7,length(bins)-1), seq(0.1,0.8,
 length=length(bins)-1), paste(bins[-length(bins)],
 "-",bins[-1]),cex=1)
text(0.5,0.95, "Bins", cex=1)

# Now plot the Germany map with the three colors. You will see that it is not 
# that informative.
par(mfrow=c(1,1))
plot(shape.shp,col=colors)


################################################################################
#                                                                              #
# 3. Code for nice map with density function in right                          #
# upper corner and legend in the right lower                                   #
# corner. Intervals are defined by interval library                            #
#                                                                              #
################################################################################

################################################################################
# 3.1 Control panel for specifications                                         #
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

# To define text for the optional title (currently disactivated in the code)
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
use_sprintf <- FALSE

# For the numbers displayed in the legend
spfleg <- "%3.1f" 


################################################################################
# 3.2 Code to produce the map                                                  #
################################################################################

# In the code in section 3.2 you might just need to change some layout-
# specifications in case you want to adapt it to your map. Apart from this, the 
# code should run automatically and use the specifications which were defined 
# above in the control panel. Here we derive the bins based on chosen 
# specifications.

# Export funtion
#tiff(file="TFR_Germany.tif",width = 2400, height = 2400, res=300, 
#     compression="lzw")
#png(file="TFR_Germany.png",width = 2400, height = 2400, res=300)
   catg <- classIntervals(varofint, nclassint, 
                          style = catgv)
   # Here we derive the colors based on chosen specifications
   colpal <- brewer.pal(nclassint,colorv)
   # Here we use the find colors function, which works similar to the find interval
   # function, but requires an object defined by the classInterval function 
   color <- findColours(catg,colpal)
   color[is.na(varofint)] <- NAcolor 
   # Obtaining our bins
   bins <- catg$brks
   # Get number of bins 
   lb <- length(bins)

   # Define layout of map
   mylayout <- layout(matrix(c(1,2),1,2), widths=c(1,0.25), 
                      heights=1)
   # mylayout <- layout(matrix(c(1,2,3),1,3), 
   #                    widths=c(1,0.25,0.25), heights=1)
   # mylayout <- layout(matrix(c(1,1,0,2), 2, 2, 
   #                    byrow=TRUE), respect=TRUE)
   layout.show(mylayout)

   # Define margins of plot to be equal to "0"
   par(mar=c(0,0,martop,0))

   # Plot map in first window
   plot(shape.shp, col=color,border=borcol)

   # Optional title 
   if (inctitle==TRUE) {title(titletext, cex=titlesiz)}

   # With the next plot command we start to plot the legend 
   # in the second window. The first plot is empty
   plot(c(0,1),c(0,1),col=NA,axes=F, xlab="", ylab="")

   # Now we use points to create the legend
   # With the rep command we define the x-position of the  
   # points, with the seq command the y-position of the
   # points
   points(rep(0.1,lb-1), seq(0.1,lb*0.025+0.1,length=lb-1), 
          pch=22, cex=2.5, col=borcol,bg=colpal)

   # Here we add the lower and upper bins of the categories
   if (use_sprintf==FALSE) {text(rep(0.6,lb-1), seq(0.1,lb*0.025+0.1,
                                                    length=lb-1),
      paste(round(bins[-lb],lar),"-",round(bins[-1],lar)), cex=ltsiz)}
   # Alternative, more sophisticated way to print the numbers,
   # so that zeros are not omitted in decimals at the last position
    if (use_sprintf==TRUE) {text(rep(0.6,lb-1), seq(0.1,lb*0.025+0.1,
                                                    length=lb-1),
       paste(sprintf(spfleg,round(bins[-lb],lar)),"-",
           sprintf(spfleg,round(bins[-1],lar))), cex=ltsiz)}

   # And finally we add the name of the variable 
   text(0.5, lb*0.025+0.2+tloc, var.name, cex=vnsiz, font=2)

   # The density plot is set independent of the layout 
   # function. In a first step we define the part of 
   # the plot window where we want to plot the 
   # density plot
   op <- par(mar=c(1,1,1,1), fig=c(0.75,0.95,0.75,0.95), 
         new=TRUE, las=1)
   # With fig we define were the plot is put. Enclosed an
   # alternative specification which is actually inside
   # the map
   # op = par(mar=c(1,1,1,1), fig=c(0.25,0.45,0.25,0.45), 
   # new = TRUE, las=1)

   # Here we plot the density the function
   plot(density(varofint,na.rm=TRUE), xlab="", ylab="", axes=F, lwd=3, 
       main="")
   # In another step we add a x-axis with min and max values
   # of the variable of interest 
   if (use_sprintf==FALSE) {axis(1, seq(round(min(varofint,na.rm=TRUE),dar),
       round(max(varofint,na.rm=TRUE),dar), length=2),cex.axis=daxsiz)}
   if (use_sprintf==TRUE) {axis(1, sprintf(spfleg,
                                   seq(round(min(varofint,na.rm=TRUE),dar),
                                       round(max(varofint,na.rm=TRUE),dar), 
                                       length=2)),cex.axis=daxsiz)}
   # Now we add vertical lines showing the bins used to 
   # categorize the data
   abline(v = bins, col=2, lwd=1.5, lty=3)
   # With this command we finalize the density plot 
   par(op)
#dev.off()

################################################################################
#                                                                              #
# 4 Simple code in case you want to produce fast a number of maps              #
#                                                                              #
################################################################################

# Code is adapted for four maps, but you can adapt it in case you want to 
# produce more or less maps

################################################################################
#                                                                              #
# 4.1 Control panel in which you can change the specifications such as which   #
#     variables you want to plot with which number of categories and which     #
#     categorisation/color schemes. You can also change the position of the    #
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
# 4.2 Plotting code to run the map. Here you would just need to make changes   #
#     to the first line in case you want to respecify your plot window.        #
#     mfrow=c(2,2) divides it in two rows and two columns, but also other      #
#     specifications are possible (e.g. mfrow=c(1,3) if you want to divide the #
#     plot window in one row and three columns. Apart from this the code       #
#     should run automatically based on the specifications set in 3.1          #
#                                                                              #
################################################################################

# Export function
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
#dev.off()
