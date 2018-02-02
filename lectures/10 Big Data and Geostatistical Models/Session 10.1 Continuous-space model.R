################################################################################
#                                                                              #
# IDEM 156 - Spatial Demography Course 2018                                    #                          
# R                                                                            #
#                  Sections 1 and 2                                            #
#     CODE: FERNANDO COLCHERO and DALIA A. CONDE                               #
#                                                                              #
################################################################################

# Erase all objects in workspace
rm(list=ls(all=TRUE))

#Load libraries
library(RColorBrewer)
library(sp)
library(maptools)
library(spdep)

################################################################################
# 1) READ IN DATA:                                                             #
#                                                                              #
# Now we are working with GRIDS (raster data) and POINT data                   #
# we are using a new function: readAsciiGrid, which deals with raster files,   #
# also known as "grids."                                                       #
#                                                                              #
################################################################################

distr      <- readAsciiGrid("SIMUL.DATA/GRID/districts")
income     <- readAsciiGrid("SIMUL.DATA/GRID/income")

#reads a shape file (point)
childc     <- readShapePoints("SIMUL.DATA/SHP/childcare") 


#reads a shape file (point)
house      <- readShapePoints("SIMUL.DATA/SHP/household") 


# Plot layers:
#Note we use "mar" for margins possition: 1= bottom, 2 =left, 3 =top & 4=right
par(mfrow=c(2,2), mar=c(4,4,4,2))


# Districts
image(distr, col=brewer.pal(11,"Spectral"))
# Main text, position 3 = top, line =heingt
mtext("Districts", 3, line=1) 

# Average income
image(income,col=brewer.pal(9,"Greens")[-1]) 
mtext("Income per distr.", 3, line=1)

# Childcare facilities
plot(childc, pch=19)
mtext("Childcare facilities", 3, line=1)

# House
plot(house, pch=24)
mtext("Households", 3, line=1)

# Extract data from GIS Layers
# Districts grid
# Number of cells
nd         <- length(income$income) 

# Range of income levels
inc        <- sort(unique(income$income)) 

# Maximum income
ninc       <- max(inc) 

dcol       <- rev(grey(1:ninc/(ninc+2)))[income$income]
sort(unique(dcol))

# Household:
babies     <- house$Children
# Existing ranks
rbabies    <- sort(unique(babies))
rbabies
# Number of existing ranks
nbabies    <- length(rbabies)
nbabies

# Plot layers:
layout(matrix(c(1,2), 1, 2), heights = 1, widths=c(1,0.25))
# Margins for the map
par(mar=c(2,2,2,2))
# Fills up the first ploting window
image(income, col=dcol) 
# Cex refers to the circle size
points(house, pch=19, col='dark blue', cex=babies+1) 
points(childc, pch=15, cex=2, col='white')
points(childc, pch=15, cex=1, col='dark red')

# Legend
par(mar=c(0,0,0,0))
# Empty plot
plot(x = c(0,1), y = c(0,1), col = NA, ann = F, axes = F) 
points(x = rep(0.2,ninc), y = inc * 0.05, col = rev(grey(1:ninc/(ninc+2))), 
       pch=15, cex=5)
text(x = rep(0.6,ninc), y = inc * 0.05, labels = 1:ninc, cex = 1.5)
text(0.5, 0.35, "Income\nscale", cex=1.5)

points(rep(0.2,nbabies), 0.5+rbabies*0.05, pch=19, col='dark blue', 
       cex=rbabies+1)
text(rep(0.6,nbabies), 0.5+rbabies*0.05, 1:nbabies-1, cex=1.5)
text(0.5, 0.775, "Children", cex=1.5)

points(0.5, 0.85, pch=15, cex=3)
points(0.5, 0.85, pch=15, cex=2.5, col='white')
points(0.5, 0.85, pch=15, cex=1.25, col='dark red')
text(0.5,0.95,"Childcare\nfacilities", cex=1.5)


################################################################################
#                                                                              #
# 2) Part II Analysis                                                          #
# 1st: Number of babies as a function of                                       #
# district income and number of childcare                                      #
# facilities per district                                                      #
#                                                                              # 
################################################################################

# Childcare facilities per district
chdis     <- over(childc,distr)
chdis     <- tapply(c(rep(1,nrow(chdis)), rep(0,max(distr$districts))),
                    c(chdis$districts, unique(distr$districts)), sum)

# Childcare facilities and income associated to household
housdis   <- over(house,distr)
housdis <- housdis$districts
housch    <- chdis[housdis]
housinc <- over(house,income)
housinc <- housinc$income

# Construct design matrix (covariates) and response variable
X1        <- cbind(housch,housinc)
y         <- babies
n         <- length(y)

dataset <- cbind(y,X1)
head(dataset)

# Run maximum-likelihood analysis
# -1 in the equation specifies that the intercept term is removed
?glm
out1      <- glm(y~X1-1, family=poisson(link='log'))
summary(out1)
# Based on aggregate-level information we get the result that availability of 
# of childcare is negative related with fertility.
idh <- house$IDh
coordinates <- coordinates(house)

# Create Nearest Neighbors Weight Matrix
l.5NN <- knearneigh(coordinates, k=5,RANN=F)
nb.5NN <- knn2nb(l.5NN, row.names=idh)
plot(nb.5NN,coordinates)
nb.5NN.lw.W <- nb2listw(nb.5NN)

# Moran's I test on residuals
mt.avag.nb.5NN <- moran.test(residuals(out1), 
                             nb.5NN.lw.W,alternative="two.sided")
mt.avag.nb.5NN

################################################################################
#                                                                              #    
# 2nd: Number of babies as a function of                                       #
# district income and distance to nearest                                      #
# childcare facility                                                           #  
#                                                                              #
################################################################################

xyh       <- coordinates(house)
xyc       <- coordinates(childc)

dch       <- rep(0,n)
for(i in 1:n){
	dis    <- sqrt((xyh[i,1]-xyc[,1])^2 + (xyh[i,2]-xyc[,2])^2)
	id     <- which(dis==min(dis))[1]
	dch[i] <- dis[id]
}

# Construct design matrix (covariates)
X2        <- cbind(dch,housinc)

# Run GLM analysis
out2      <- glm(y~X2-1, family=poisson(link='log'))
summary(out2)
# Now our childcare variable is distance to childcare and we get the expected
# outcome: The further an household is away from childcare, the less children
# it has.

# Moran's I test
mt.avag.nb.5NN <- moran.test(residuals(out2), nb.5NN.lw.W,
                             alternative="two.sided")
mt.avag.nb.5NN

# Plot results:
disch     <- 0:30

par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(c(0,30), c(0,1.75), col=NA, frame.plot=FALSE, xlab="Distance", 
     ylab="Expected num. of children")
for(i in inc){
   lines(disch, exp(out2$coef[1]*disch + out2$coef[2] * i), 
         col= brewer.pal(ninc+1, "Greens")[-1][i], lwd=4)
}
legend('topright', paste("$", rev(inc)), 
       col=rev(brewer.pal(ninc+1, "Greens")[-1]), 
       lwd=4, bty='n', title="Income")
