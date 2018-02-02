################################################################################
#                                                                              #                          
# IDEM 156 - Spatial Demography Course 2018                                    #                          
# Example for Shadowline Map                                                   #                          
# Dalia A. Conde and Fernando Colchero, University of Odense                   #
# with adaptations by Sebastian Klüsener, MPIDR                                #
#                                                                              #                          
################################################################################

# Erase all objects in workspace
rm(list=ls(all=TRUE))

# Load libraries
library(sp)
library(maptools)
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
#                                                                              #                          
# 1) Load files                                                                #                          
#                                                                              #                          
################################################################################

# To read the sapefile
shape.shp <- readShapePoly("2008_w_fert_data",IDvar="GS")
shapeogr.shp <- readOGR(".", "2008_w_fert_data")
proj4string(shape.shp) <- proj4string(shapeogr.shp)


################################################################################
#                                                                              #                          
# 2) Control-Panel                                                             #
#    Choose variable and set specification for shadowline map                  #                          
#                                                                              #                          
################################################################################

# Define variable you are interested in
var.data <- shape.shp$TFR07

# Define name of the variable to be printed in the map
var.name <- c("Total\nFertility\nRate")

# Here you choose the number of classes (bins-1)
# If you would like to use more than three classes, you would also need to add
# information to the angle and dens(ity) vectors below as they are currently 
# just containin information for up to three different regions.
nclassint <- 3

# Define angle of the lines in the polygons
# Here we define the angle in which the lines are plotted within each regional 
# polygon. In this case we have three categories, so we need to define for each 
# category an angle. The following values are most commonly used:
# 0  (horizontal)
# 90 (vertical)
# 45 (diagonal, upward)
# -45 (diagonal, downward)
angle1 <- c(0,-45,45)
# As we want to plot crossing lines for the third category, we need to plot 
# above the first map based on the angle-parameters a second map for which we 
# used the following parameters. In contrast to angle1, we choose for the last
# category instead of a "45" a "-45". This allows us to plot first lines along
# one diagonal, and in the second step crossing lines in the other diagonal.
angle2 <- c(0,0,-45)

# Define density of lines in the polygons
# In this case we have three categories, so we need to define for each category
# a density value. The higher the density value, the closer the distance
# between one line and the next line that are drawn within the polygon. If you 
# define a density of 0, no lines will be drawn.
dens1 <- c(0,15,30)
dens2 <- c(0,0,30)

# Define line-width of the lines in the polygons. Try alternative specifications
# to make lines bigger or smaller
lwd1 <- c(1,1,1)
lwd2 <- c(1,1,1)

# Define shape of lines
# (0=blank, 1=solid (default), 2=dashed, 3=dotted, 4=dotdash, 5=longdash, 
# 6=twodash)
lty1 <- c(1,1,1)
lty2 <- c(1,1,1)


# Here you can specify different categorizations
varofint <- var.data
catg <- classIntervals(varofint, nclassint, 
                       style = "quantile")
# Examples for alternative categorizations
# ?classIntervals 
# catg <- classIntervals(varofint, nclassint, 
#  style = "equal")
# catg <- classIntervals(varofint, nclassint,
#  style = "jenks")


# Instead of colors we are now assigning the above defined
# angle, density, linewidth and linetype information to our data values.
# If e.g. the first region has a value in the second category
# (between the second and third bin, it would get the specifications in 
# the second position of the density, angle, line width and linetype 
# vectors)

# We create a sequence of natural numbers representing the number of 
# classes 
shadepal <- c(1:nclassint)

# Now we use this sequence of numbers to derive with the findInterval
# function in which category each of our values is located.
catmem <- shadepal[findInterval(var.data, catg$brks, 
                              rightmost.closed=T)]

# Now we assign the shadowline-specifications based on the category
# each of the value is located in. 
dens1d <- dens1[catmem]
dens2d <- dens2[catmem]
angle1d <- angle1[catmem]
angle2d <- angle2[catmem]
lwd1d <- lwd1[catmem]
lwd2d <- lwd2[catmem]
lty1d <- lty1[catmem]
lty2d <- lty2[catmem]

# Obtaining our bins
bins <- catg$brks
# Get number of bins 
lb <- length(bins)

# Export funtion
tiff(file="Shadowline_map.tif",width = 2400, height = 2400, res=300, 
     compression="lzw")
#png(file="Shadowline_map.png",width = 2400, height = 2400, res=300)
    # Define layout of map
    mylayout <- layout(matrix(c(1,2),1,2), widths=c(1,0.25), 
                       heights=1)

    # Define margins of plot to be equal to "0"
    par(mar=c(0,0,0,0))

    # Plot map in first window
    plot(shape.shp, density=dens1d,angle=angle1d,lwd=lwd1d,lty=lty1d)
    
    # Optional: Do a cross-hatching for the highest category
    plot(shape.shp, density=dens2d,angle=angle2d,lwd=lwd2d,lty=lty2d,add=T)

    # With the next plot command we start to plot the legend 
    # in the second window. The first plot is empty
    plot(c(0,1),c(0,1),col=NA,axes=F, xlab="", ylab="")

    # Now we use points to create the background of the legend
    # With the rep command we define the x-position of the 
    # points, with the seq command the y-position of the
    # points

    # points(rep(0.1,lb-1), seq(0.1,lb*0.025+0.1,length=lb-1), 
    #        pch=15, cex=2.5, col="black")
    # In the next step we plot above the background the colored
    # points according to the colors used in the map
    central_positions_x <- c(rep(0.1,lb-1))
    central_positions_y <- c(seq(0.1,lb*0.025+0.1,length=lb-1))
    for (i in 1:3) {
        cpx <- central_positions_x[i]
        cpy <- central_positions_y[i]
        dfcx <- 0.05
        dfcy <- 0.01
        polygon(c(cpx-dfcx,cpx+dfcx,cpx+dfcx,cpx-dfcx),
                c(cpy+dfcy,cpy+dfcy,cpy-dfcy,cpy-dfcy),
                density=dens1[i],angle=angle1[i],lwd=lwd1[i],
                lty=lty1[i], col="black")
        polygon(c(cpx-dfcx,cpx+dfcx,cpx+dfcx,cpx-dfcx),
                c(cpy+dfcy,cpy+dfcy,cpy-dfcy,cpy-dfcy),
                density=dens2[i],angle=angle2[i],lwd=lwd2[i],
                lty=lty2[i], col="black")
    }
        
    #points(rep(0.1,lb-1), seq(0.1,lb*0.025+0.1,length=lb-1), 
    #       pch=15, cex=2, col=colpal)
    # Here we add the lower and upper bins of the categories
    text(rep(0.6,lb-1), seq(0.1,lb*0.025+0.1,length=lb-1),
    paste(round(bins[-lb],2),"-",round(bins[-1],2)), cex=1)
    # And finally we add the name of the variable 
    text(0.5, lb*0.025+0.2, var.name, cex=1.2, font=2)

    # The density plot is set independent of the layout 
    # function. In a first step we define the part of 
    # the plot window where we want to plot the 
    # density plot
    op <- par(mar=c(1,1,1,1), fig=c(0.75,0.95,0.75,0.95), 
          new = TRUE, las=1)
    # With fig we define were the plot is put. Enclosed an
    # alternative specification which is actually inside
    # the map
    #op = par(mar=c(1,1,1,1), fig=c(0.25,0.45,0.25,0.45), 
    # new = TRUE, las=1)

    # Here we plot the density the function
    plot(density(varofint,na.rm = T), xlab="", ylab="", axes=F, lwd=3, 
         main="")
    # In another step we add a x-axis with min and max values
    # of the variable of interest 
    axis(1, seq(round(min(varofint,na.rm = T),2),
                round(max(varofint,na.rm = T),2), length=2))
    # Now we add vertical lines showing the bins used to 
    # categorize the data
    abline(v = bins, col="black", lwd=1.5, lty=3)
    # With this command we finalize the density plot 
    par(op)
dev.off()
