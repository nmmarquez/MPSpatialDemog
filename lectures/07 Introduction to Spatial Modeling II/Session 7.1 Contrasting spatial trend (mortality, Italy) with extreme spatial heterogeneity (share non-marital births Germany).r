################################################################################
#                                                                              #
# IDEM 156 - Spatial Demography Course 2018                                    #                          
# Contrasting Spatial Trend (Mortality, Italy) with                            #
# Extreme Spatial Heterogeneity (Non-Marital Births)                           #
# Sebastian Kl?sener, Francesco Lagona, MPIDR                                  #
#                                                                              #
################################################################################

# Erase all objects in workspace
rm(list=ls(all=TRUE))

# Load libraries
library(spdep)
library(RColorBrewer)
library(maptools)
library(rgdal)


################################################################################
# 1)   Italian Dataset (data and parts of the code by Francesco Lagona)        #
# 1.1) Importing and preparing Italian data                                    #
################################################################################


#loading the italian provinces
map.italy <- readShapePoly("ita_prov",IDvar="COD")
# the names of the provinces
prov.names <- map.italy$NOME
dat <- read.table('cancer.deaths.txt',header=T)   
dat$prov <- rep(rep(map.italy$COD,each=101),2)

dat.males <- subset(dat,sex=="M")
males.obs.deaths <- tapply(dat.males$n,dat.males$prov,sum)
dat.males$risks <- rep(tapply(dat.males$n,dat.males$age,sum)/
 tapply(dat.males$E,dat.males$age,sum),103)
dat.males$expected.n <- dat.males$risks*dat.males$E
males.SMR <- males.obs.deaths/tapply(dat.males$expected.n,dat.males$prov,sum)
males.dat <- as.data.frame(males.SMR,row.names=names(males.SMR))
males <- SpatialPolygonsDataFrame(map.italy, males.dat, match.ID = TRUE)
coord <- coordinates(males)
plot(males)
                      
y <- males$males.SMR
n <- length(y)

var.data <- data.frame(y)
var.names <-c("Mortality Males")
m <- length(var.data)
colors <- matrix(nrow=n,ncol=m, data=rep(0,n*m))

for (i in 1:m) {
varofint <- var.data[,i]

lb       <- 10
bins     <- seq(min(varofint), max(varofint), length=lb)
colpal   <- brewer.pal(lb-1, "YlOrRd")
colors[,i]   <- colpal[findInterval(varofint, bins, rightmost.closed=T)]
plot(map.italy, col=colors)
title(var.names[[i]], cex=1)
legend("topright",fill=colpal,
 legend=paste(round(bins[-lb],2),"-",round(bins[-1],2)),cex=0.5)
}
                             
# Plot by Latitude and Longitude
par(mfrow=c(2,2))
for (i in 1:m) {
varofint <- var.data[,i]

lb       <- 10
bins     <- seq(min(varofint), max(varofint), length=lb)
colpal   <- brewer.pal(lb-1, "YlOrRd")
colors[,i] <- colpal[findInterval(varofint, bins, rightmost.closed=T)]
plot(map.italy, col=colors)
title(var.names[[i]], cex=1)
legend("topright",fill=colpal,
 legend=paste(round(bins[-lb],2),"-",round(bins[-1],2)),cex=0.5)
}


plot(coord[,1],y, xlab="West <-> East", ylab="Mortality Males Italy",
  main="Trends by Longitude")
abline(coef(lm(y~coord[,1])))
plot(coord[,2],y, xlab="South <-> North", ylab="Mortality Males Italy", 
 main="Trends by Latitude")
abline(coef(lm(y~coord[,2])))

easting <- coord[,1]
northing <- coord[,2]
eastwest <- lm(y~easting)
northsouth <- lm(y~northing)
summary(eastwest)
summary(northsouth)
geo <- lm(y~easting+northing)
summary(geo)

# Create first order Rook contiguity weight matrix
neigh.italy.r <- poly2nb(males,queen=F) 

# Create a listwise object
listw.italy <- nb2listw(neigh.italy.r,style="W")
listw.italy.B <- nb2listw(neigh.italy.r,style="B")

# Calculate SAR-Model
sar.it <- spautolm(as.vector(y)~easting+northing, listw=listw.italy,family="SAR")

# Calculate Error-Model
error.it <- errorsarlm(as.vector(y)~easting+northing, listw=listw.italy)

# Calculate CAR-Model
sar.b.it <- spautolm(as.vector(y)~easting+northing, 
                     listw=listw.italy.B,family="SAR")
car.b.it <- spautolm(as.vector(y)~easting+northing, 
                     listw=listw.italy.B,family="CAR")

# Similiar autocomes for SAR and Error Model
summary(sar.it)
summary(error.it)

# Outcomes for SAR and CAR models differ slightly
summary(sar.b.it)
summary(car.b.it)


################################################################################
# 2) Importing German data                                                     #
################################################################################

par(mfrow=c(1,1))
shape.shp <- readShapePoly('2004_06', IDvar="KREIS_KENN")
shapeogr.shp <- readOGR(".", "2004_06")
proj4string(shape.shp) <- proj4string(shapeogr.shp)
shape.ll <- spTransform(shape.shp, CRS("+proj=longlat +datum=WGS84"))
plot(shape.shp)                   

# Derive information on internal id and centroids. In this case we turn the id-
# values first into a numeric vector, which is easier to use, if we want to
# match data to the shapefile.
head(shape.shp$KREIS_KENN)
is.vector(shape.shp$KREIS_KENN)
id <- as.numeric(paste(shape.shp$KREIS_KENN))/1000
head(id)
is.vector(id)

n <- length(id)
sq <- c(1:n)
sq_text <- paste(sq)
coord <- coordinates(shape.shp)
coordsm <- coordinates(shape.shp)/100000


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
names(shape.shp)
y1 <- shape.shp$SNM05

# Defining position of east and west districts
west.dist <- c(1:326)
east.dist <-  c(327:429)

# Sample West
y1west <- shape.shp$SNM05[shape.shp$DUM_E==0]
c1west <- coord[,1][shape.shp$DUM_E==0]
c2west <- coord[,2][shape.shp$DUM_E==0]

# Sample East
y1east <- shape.shp$SNM05[shape.shp$DUM_E==1]
c1east <- coord[,1][shape.shp$DUM_E==1]
c2east <- coord[,2][shape.shp$DUM_E==1]

#Plot by Latitude and Longitude
par(mfrow=c(2,2))
var.data <- data.frame(y1)
var.names <-c("Non-Marital Births Ratio 2005")
m <- length(var.data)
colors <- matrix(nrow=n,ncol=m, data=rep(0,n*m))

for (i in 1:m) {
varofint <- var.data[,i]

lb       <- 10
bins     <- seq(min(varofint), max(varofint), length=lb)
colpal   <- brewer.pal(lb-1, "YlOrRd")
colors   <- colpal[findInterval(varofint, bins, rightmost.closed=T)]
plot(shape.shp, col=colors)
title(var.names[[i]], cex=1)
legend("topright",fill=colpal,
 legend=paste(round(bins[-lb],2),"-",round(bins[-1],2)),cex=0.5)
}


plot(coord[,1],y1, xlab="West <-> East", ylab="Share Non-Marital Births 2005", 
 main="Trends by Longitude")
# Points for Sample West and Sample East in different colors
points(c1west,y1west, col="orange2")
points(c1east,y1east, col="navyblue")
abline(coef(lm(y1~coord[,1])), lwd="1")
# Lines for Sample West and Sample East
abline(coef(lm(y1west~c1west)), col="red", lty=2)
abline(coef(lm(y1east~c1east)), col="red", lty=2)
legend("topleft", c("East","West"), pch=1, 
 col=c("navyblue","orange2"), bg="white")
plot(coord[,2],y1, xlab="South <-> North", 
 ylab="Share Non-Marital Birth", main="Trends by Latitude")
# Points for Sample West and Sample East in different colors
points(c2west,y1west, col="orange2", lty=2)
points(c2east,y1east, col="navyblue", lty=2)
abline(coef(lm(y1~coord[,2])), lwd="1")
# Lines for Sample West and Sample East
abline(coef(lm(y1west~c2west)), col="red", lty=2)
abline(coef(lm(y1east~c2east)), col="red", lty=2)
legend("topleft", c("East","West"), pch=1, 
 col=c("navyblue","orange2"), bg="white")


################################################################################
# 2) Creating Neighborhood Weight Matrix                                       #
################################################################################

nb.FOQ <- poly2nb(shape.shp, queen=TRUE)
summary(nb.FOQ)

# Code to connect R?gen
nb.FOQ <- poly2nb(shape.shp, queen=TRUE)
summary(nb.FOQ)

# Code to connect R?gen
nb.FOQ.cor <- nb.FOQ
nb.FOQ.cor[[362]]  <- as.integer(c(350, 358))
nb.FOQ.cor[[350]]  <- as.integer(c(358, 362))
nb.FOQ.cor[[358]]  <- as.integer(c(350,352,353,354,360,362))
nb.FOQ.lw.W <- nb2listw(nb.FOQ.cor,style="W")
nb.FOQ.lw.B <- nb2listw(nb.FOQ.cor,style="B")

################################################################################
# 3) Running Model                                                             #
################################################################################

y1 <- shape.shp$SNM05
sar <- spautolm(as.vector(y1)~coordsm[,1]+coordsm[,2], listw=nb.FOQ.lw.W, 
 family="SAR")
summary(sar)

sar.b <- spautolm(as.vector(y1)~coordsm[,1]+coordsm[,2], listw=nb.FOQ.lw.B, 
                  family="SAR")
car.b <- spautolm(as.vector(y1)~coordsm[,1]+coordsm[,2], listw=nb.FOQ.lw.B, 
            family="CAR")
summary(sar.b)
summary(car.b)

################################################################################
# 3) Contrasting Germany and Italy                                             #
################################################################################

# Running Moran's I Test
morantest.dep.it <- moran.test(as.vector(males$males.SMR), listw.italy)
morantest.model.it <- moran.test(sar.it$fit$residuals, listw.italy)
morantest.dep.de <- moran.test(y1, nb.FOQ.lw.W)
morantest.model.de <- moran.test(sar$fit$residuals, nb.FOQ.lw.W)

morantest.dep.it
morantest.model.it
morantest.dep.de
morantest.model.de

# Looking at Moran's I Plot
par(mfrow=c(2,2))
moran.plot(as.vector(males$males.SMR), listw.italy, main="Male Mortality Italy")
moran.plot(sar.it$fit$residuals, listw.italy, main="Residuals - SAR Italy")
moran.plot(y1, nb.FOQ.lw.W, main="NMB Ratio Germany")
moran.plot(sar$fit$residuals, nb.FOQ.lw.W, main="Residuals - SAR Germany")

# Maps for Italy
#png(file="Italy.png",width = 2800, height = 2000, res=300)
sar.it <- spautolm(as.vector(y)~easting+northing, listw=listw.italy,family="SAR")
var.data <- data.frame(y,sar.it$fit$signal_trend, sar.it$fit$signal_stochastic, 
 sar.it$fit$residuals)
var.names <-c("observed Male Mortality", "large scale", "small scale", 
 "residuals")
m <- length(var.data)
colors <- matrix(nrow=n,ncol=m, data=rep(0,n*m))
               
par(mfrow=c(2,2))
for (i in 1:m) {
varofint <- var.data[,i]
mean <- mean(varofint)
sd <- sd(varofint)
min <- min(varofint)
max <- max(varofint)
bins     = c(min,mean-sd, mean-sd*0.5, mean, mean+sd*0.5, mean+sd, max)
lb       = c(length(bins)-1)
med      = median(lb)
colpal   = c(brewer.pal(length(bins), "PiYG"))
colpal   = colpal[-med]
colors   = colpal[findInterval(varofint, bins, rightmost.closed=T)]
plot(map.italy, col=colors)
title(var.names[[i]], cex=1)
legend("topright",fill=colpal,
 legend=paste(round(bins[-length(bins)],2),"-",round(bins[-1],2)),cex=0.8)
}
#dev.off()


# Maps for Germany
#png(file="Germany.png",width = 2800, height = 2000, res=300)
sar <- spautolm(as.vector(y1)~coordsm[,1]+coordsm[,2], listw=nb.FOQ.lw.W, 
                family="SAR")
var.data <- data.frame(y1,sar$fit$signal_trend, sar$fit$signal_stochastic, 
 sar$fit$residuals)
var.names <-c("observed NMB Ratio", "large scale", "small scale", "residuals")
m <- length(var.data)
colors <- matrix(nrow=n,ncol=m, data=rep(0,n*m))

par(mfrow=c(2,2))
for (i in 1:m) {
varofint <- var.data[,i]
mean <- mean(varofint)
sd <- sd(varofint)
min <- min(varofint)
max <- max(varofint)
bins     = c(min,mean-sd, mean-sd*0.5, mean, mean+sd*0.5, mean+sd, max)
lb       = c(length(bins)-1)
colpal   = c(brewer.pal(length(bins), "PiYG"))
colpal   = colpal[-med]
colors   = colpal[findInterval(varofint, bins, rightmost.closed=T)]
plot(shape.shp, col=colors)
title(var.names[[i]], cex=1)
legend("topright",fill=colpal,
 legend=paste(round(bins[-length(bins)],2),"-",round(bins[-1],2)),cex=0.8)
}
#dev.off()


# -> Under the condition of extreme spatial heterogeneity with clear-cut borders 
# between regions, the large-scale spatial trend might occur in discrete intervals,
# and not continuous across space (as it is the case in Germany)


y1 <- shape.shp$SNM05
sar <-spautolm(as.vector(y1)~shape.shp$DUM_E+shape.shp$AVAGE05+coordsm[,1]+coordsm[,2], listw=nb.FOQ.lw.W,
 family="SAR")
summary(sar)

var.data <- data.frame(y1,sar$fit$signal_trend, sar$fit$signal_stochastic, 
 sar$fit$residuals)
var.names <-c("SNM07", "large scale", "small scale", "residuals")
m <- length(var.data)
colors <- matrix(nrow=n,ncol=m, data=rep(0,n*m))

# Model Germany with Dummy
var.data <- data.frame(y1,sar$fit$signal_trend, sar$fit$signal_stochastic, 
 sar$fit$residuals)
var.names <-c("observed NMB Ratio", "large scale", "small scale", "residuals")
m <- length(var.data)
colors <- matrix(nrow=n,ncol=m, data=rep(0,n*m))

par(mfrow=c(2,2))
for (i in 1:m) {
varofint <- var.data[,i]
mean <- mean(varofint)
sd <- sd(varofint)
min <- min(varofint)          
max <- max(varofint)
#if (i==1) 
bins = c(min,mean-sd*0.5, mean, mean+sd*0.5, max)
#bins     = c(min,mean-sd, mean-sd*0.5, mean, mean+sd*0.5, mean+sd, max)
# if (i!=1) bins = c(min,mean-sd, mean-sd*0.5, mean, mean+sd*0.5, mean+sd, max)
lb       = c(length(bins)-1)
colpal   = c(brewer.pal(length(bins), "PiYG"))
colpal   = colpal[-med]
colors   = colpal[findInterval(varofint, bins, rightmost.closed=T)]
plot(shape.shp, col=colors)
title(var.names[[i]], cex=1)
legend("topright",fill=colpal,legend=paste(round(bins[-length(bins)],2),"-",
 round(bins[-1],2)),cex=0.5)
}

# Be very careful, if you are faced with Extreme Spatial Heterogeneity
# Examples: Highly spatially segragated cities (e.g. in North and South America)
#           West and East Germany
#           West and East Ukraine
#           North and South Korea
# In such situations you might easily treat a proxy as a covariate and
# consider the association as causal, while it is actually just an association.
