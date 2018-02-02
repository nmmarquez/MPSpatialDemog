################################################################################
#                                                                              #
# IDEM 156 Spatial Demography Course 2018                                      #                          
# Multi Level Models and Spatial Autocorrelation                               #
# Sebastian Klüsener, MPIDR                                                    #
#                                                                              #
################################################################################

# Opening libraries
library(spdep)
library(RColorBrewer)
library(lme4)
library(maptools)
library(rgdal)

# Set your working directory
# If you need to change the working directory, you can use the following code:
# Set drive
main <- c("N")
# Set path to session folders
path <- ":/IMPRSD/IDEM 156 Spatial"
# Define session folder
path2 <- "/02 Sessions/08 Spatial Multi-level Modeling"
# Set working directory
setwd(paste(main,path,path2,sep=""))


################################################################################
# 1) Import simulated data                                                     #
#    Ind.-level:                                                               #
#    Registered Severe Heart Attacks                                           #
#    - Did Person die within the first 24 hours after the attack               #
#                                                                              #
################################################################################

# Open Individual level data
Ind.data <-  read.table(file="Ind.data.csv",sep=",", header=TRUE)

colnames(Ind.data) <-c("Autonumber","vid","v1","countw")

# Open District level data
Dis.data <-  read.table(file="Dis.data.csv",sep=",", header=TRUE)

# Open Shape file
shape.shp <- readShapePoly('shape.shp', IDvar="KREIS_KENN.1")
shapeogr.shp <- readOGR(".", "shape")
proj4string(shape.shp) <- proj4string(shapeogr.shp)
shape.ll <- spTransform(shape.shp, CRS("+proj=longlat +datum=WGS84"))


id <- shape.shp$KREIS_KENN.1
n <- length(shape.shp[,1])

# Match district data to shape file
data.id <- Dis.data$ID
o <- match(id, data.id)
mydata1 <- Dis.data[o,]
row.names(mydata1) <- id
shape.shp <- spCbind(shape.shp, mydata1)
names(shape.shp)


################################################################################
# 2) Descriptives                                                              #
################################################################################

# Overall mortality rate within 24 hours after the heart attack
per.dying <- (sum(Ind.data$v1)/ length(Ind.data$v1))*100
per.dying

# Summing up the the deaths and heart attacks by district
xsum <- rowsum(Ind.data$v1, Ind.data$vid)
ysum <- rowsum(Ind.data$count, Ind.data$vid)
per.dying.d <- (xsum/ysum*100)
dframe <- data.frame(id,per.dying.d)
dframe

# Mapping Data on mortality after heart attack
n <- length(dframe$id)
par(mfrow=c(1,2))
var.data <- data.frame(per.dying.d)
var.names <-c("Mortality Rate within 24 hours after heart attack")
m <- length(var.data)
colors <- matrix(nrow=n,ncol=m, data=rep(0,n*m))
for (i in 1:m) {
    varofint <- var.data[,i]
    mean <- mean(varofint)
    sd <- sd(varofint)
    min <- min(varofint)          
    max <- max(varofint)
    bins = c(min,mean-sd*0.5, mean, mean+sd*0.5, max)
    lb       = c(length(bins)-1)
    colpal   = brewer.pal(length(bins-1), "Reds")
    colors   = colpal[findInterval(varofint, bins, rightmost.closed=T)]
    plot(shape.shp, col=colors)
    title(var.names[[i]], cex=1)
    legend("topright",fill=colpal,legend=paste(round(bins[-length(bins)],2),
     "-",round(bins[-1],2)),cex=0.5)
}

# Map data on local quality of emergency care institutions
var.data <- data.frame(shape.shp$HEALTH)
var.names <-c("Quality of Emergency Care")
m <- length(var.data)
colors <- matrix(nrow=n,ncol=m, data=rep(0,n*m))
for (i in 1:m) {
    varofint <- var.data[,i]
    mean <- mean(varofint)
    sd <- sd(varofint)
    min <- min(varofint)          
    max <- max(varofint)
    bins <- c(min,mean-sd*0.5, mean, mean+sd*0.5, max)
    lb <- c(length(bins)-1)
    med <- median(1:c(lb+1))
    colpal <- brewer.pal(length(bins), "PiYG")[-med]
    colors   = colpal[findInterval(varofint, bins, rightmost.closed=T)]
    plot(shape.shp, col=colors)
    title(var.names[[i]], cex=1)
    legend("topright",fill=colpal,legend=paste(round(bins[-length(bins)],2),
     "-",round(bins[-1],2)),cex=0.5)
}

# Creating neighborhood weight file based on First Order Queen definition
nb.FOQ <- poly2nb(shape.shp, queen=TRUE)
nb.FOQ.lw.W <- nb2listw(nb.FOQ)

# Calculating Moran's I on aggregates of available variables
mt.hall <- moran.test(per.dying.d, nb.FOQ.lw.W)
mth.health <- moran.test(shape.shp$HEALTH, nb.FOQ.lw.W)
mt.hall 
mth.health


################################################################################
# 3)   Modeling                                                               #
################################################################################

# Preparation of Dataset for Multi-Level Model in order to control for 
# the quality of the health emergency system
colnames(Dis.data) <- c("vid","D","Health")
ML.data <- merge(Ind.data, Dis.data, by="vid")
District <- ML.data$vid
Death <- ML.data$v1
Qual_Emer_Care <- ML.data$Health


################################################################################
# 3.1) Simple logistic model, in which we do not take account for the multi-   #
#      level-structure of our data                                             #
################################################################################

model <- glm(Death ~ Qual_Emer_Care, family = binomial())
summary(model)
names(model)

# Derive model residuals
res <- data.frame(Ind.data$vid, model$residuals)

# Split data into groups to derive calcullate the mean values of the 
# residuals
res.sp <- split(res, res$Ind.data.vid)
resd <- c(rep(0,n))
for (i in 1:n) {
   resd[i] <- mean(res.sp[[i]]$model.residuals)
}                
resd

# Map the mean values of the residuals
n <- length(dframe$id)
par(mfrow=c(1,2))
var.data <- data.frame(resd)
var.names <-c("Residuals - Health")
m <- length(var.data)
colors <- matrix(nrow=n,ncol=m, data=rep(0,n*m))
for (i in 1:m) {
    varofint <- var.data[,i]
    mean <- mean(varofint)
    sd <- sd(varofint)
    min <- min(varofint)          
    max <- max(varofint)
    bins = c(min,mean-sd*0.5, mean, mean+sd*0.5, max)
    lb <- c(length(bins)-1)
    med <- median(1:c(lb+1))
    colpal <- brewer.pal(length(bins), "PiYG")[-med]
    colors <- colpal[findInterval(varofint, bins, rightmost.closed=T)]
    plot(shape.shp, col=colors)
    title(var.names[[i]], cex=1)
    legend("topright",fill=colpal,
     legend=paste(round(bins[-length(bins)],2),"-",round(bins[-1],2)),cex=0.5)
}

# Moran's I test on the mean values of the residuals
mt.res <- moran.test(resd, nb.FOQ.lw.W)
mt.res


################################################################################
# 3.2) Multi-level model, in which we account for the multi-level structure of #
#      our data                                                                # 
################################################################################

ML.Model <- glmer(formula=Death~Qual_Emer_Care+(1|District), family = binomial)          
summary(ML.Model)

# Check the residuals
res2 <- data.frame(Ind.data$vid, residuals(ML.Model))
res.sp2 <- split(res2, res2$Ind.data.vid)

# Derive for each district the average value of the individual-level residuals 
res2d <- c(rep(0,n))
for (i in 1:n) {
    res2d[i] <- mean(res.sp2[[i]]$residuals.ML.Model.)
}
res2d

# Random effects of the Model
random.effects_2 <- ranef(ML.Model)

# Moran's I of the Residuals
mt.mod <- moran.test(resd, nb.FOQ.lw.W)
mt.mod2 <- moran.test(res2d, nb.FOQ.lw.W)
mt.mod_re2 <- moran.test(random.effects_2$District[,1], nb.FOQ.lw.W)

# Comparison
# mt.mod:  Moran'I of district level resiudals of glm model with multi-level 
# structure
# mt.mod2: Moran'I of district level residuals of multi-level model
mt.mod
mt.mod2
mt.mod_re2

# Map of the Residuals
n <- length(dframe$id)
par(mfrow=c(1,2))
var.data <- data.frame(res2d)
var.names <-c("Residuals of Multi-Level Model")
m <- length(var.data)
colors <- matrix(nrow=n,ncol=m, data=rep(0,n*m))
for (i in 1:m) {
    varofint <- var.data[,i]
    mean <- mean(varofint)
    sd <- sd(varofint)
    min <- min(varofint)          
    max <- max(varofint)
    bins = c(min,mean-sd, mean, mean+sd, max)
    lb <- c(length(bins)-1)
    med <- median(1:c(lb+1))
    colpal <- brewer.pal(length(bins), "PiYG")[-med]
    colors <- colpal[findInterval(varofint, bins, rightmost.closed=T)]
    plot(shape.shp, col=colors)
    title(var.names[[i]], cex=1)
    legend("topright",fill=colpal,legend=paste(round(bins[-length(bins)],2),
     "-",round(bins[-1],2)),cex=0.5)
}

var.data <- data.frame(random.effects_2$District[,1])
var.names <-c("Random Effects of Multi-Level Model")
m <- length(var.data)
colors <- matrix(nrow=n,ncol=m, data=rep(0,n*m))
for (i in 1:m) {
  varofint <- var.data[,i]
  mean <- mean(varofint)
  sd <- sd(varofint)
  min <- min(varofint)          
  max <- max(varofint)
  bins = c(min,mean-sd, mean, mean+sd, max)
  lb <- c(length(bins)-1)
  med <- median(1:c(lb+1))
  colpal <- brewer.pal(length(bins), "PiYG")[-med]
  colors <- colpal[findInterval(varofint, bins, rightmost.closed=T)]
  plot(shape.shp, col=colors)
  title(var.names[[i]], cex=1)
  legend("topright",fill=colpal,legend=paste(round(bins[-length(bins)],2),
                                             "-",round(bins[-1],2)),cex=0.5)
}

# Introduction of lagged variable - Quality of Health Care System                             
varlag1 <- lag.listw(nb.FOQ.lw.W, shape.shp$HEALTH)
varlag2 <- (shape.shp$HEALTH*0.2)+(varlag1*0.8)

# Choose between the two definition of the lagged Quality of the
# Health Care System
varlag <- varlag2

par(mfrow=c(1,2))
var.data <- data.frame(shape.shp$HEALTH,varlag)
var.names <-c("Quality of Emergency Care", "Lagged Quality Emergency Care")
m <- length(var.data)
colors <- matrix(nrow=n,ncol=m, data=rep(0,n*m))
for (i in 1:m) {
    varofint <- var.data[,i]
    mean <- mean(varofint)
    sd <- sd(varofint)
    min <- min(varofint)          
    max <- max(varofint)
    bins = c(min,mean-sd*0.5, mean, mean+sd*0.5, max)
    lb <- c(length(bins)-1)
    med <- median(1:c(lb+1))
    colpal <- brewer.pal(length(bins), "PiYG")[-med]
    colors <- colpal[findInterval(varofint, bins, rightmost.closed=T)]
    plot(shape.shp, col=colors)
    title(var.names[[i]], cex=1)
    legend("topright",fill=colpal,legend=paste(round(bins[-length(bins)],2),
     "-",round(bins[-1],2)),cex=0.5)
}

# Conntect variable to dataset
to <- data.frame(id,varlag)
colnames(to) <- c("vid","varlag")
ML.data <- merge(Ind.data,to, by="vid")
District <- ML.data$vid
Death <- ML.data$v1
Lagged_Qual_Emer_Care <- ML.data$varlag


# Multi-level model
ML.Model1 <-glmer(formula=Death~Lagged_Qual_Emer_Care+(1|District), 
 family = binomial)          
summary(ML.Model1)

# Residuals of model
res3 <- data.frame(Ind.data$vid, residuals(ML.Model1))
res.sp3 <- split(res3, res3$Ind.data.vid)
resd3 <- c(rep(0,n))
for (i in 1:n) {
    resd3[i] <- mean(res.sp3[[i]]$residuals.ML.Model1.)
}
resd3

# Random effects of the Model
random.effects_3 <- ranef(ML.Model1)

# Moran's I test on the residuals of the multi-level model with a lagged 
# variable
mt.mod3 <- moran.test(resd3, nb.FOQ.lw.W)
mt.mod_re3 <- moran.test(random.effects_3$District[,1], nb.FOQ.lw.W)

# Map of the Residuals
n <- length(dframe$id)
par(mfrow=c(1,2))
var.data <- data.frame(resd3)
var.names <-c("Residuals of Multi-Level Model\nwith Lagged Variable")
m <- length(var.data)
colors <- matrix(nrow=n,ncol=m, data=rep(0,n*m))
for (i in 1:m) {
  varofint <- var.data[,i]
  mean <- mean(varofint)
  sd <- sd(varofint)
  min <- min(varofint)          
  max <- max(varofint)
  bins = c(min,mean-sd, mean, mean+sd, max)
  lb <- c(length(bins)-1)
  med <- median(1:c(lb+1))
  colpal <- brewer.pal(length(bins), "PiYG")[-med]
  colors <- colpal[findInterval(varofint, bins, rightmost.closed=T)]
  plot(shape.shp, col=colors)
  title(var.names[[i]], cex=1)
  legend("topright",fill=colpal,legend=paste(round(bins[-length(bins)],2),
                                             "-",round(bins[-1],2)),cex=0.5)
}

var.data <- data.frame(random.effects_3$District[,1])
var.names <-c("Random Effects of Multi-Level Model\nwith Lagged Variable")
m <- length(var.data)
colors <- matrix(nrow=n,ncol=m, data=rep(0,n*m))
for (i in 1:m) {
  varofint <- var.data[,i]
  mean <- mean(varofint)
  sd <- sd(varofint)
  min <- min(varofint)          
  max <- max(varofint)
  bins = c(min,mean-sd, mean, mean+sd, max)
  lb <- c(length(bins)-1)
  med <- median(1:c(lb+1))
  colpal <- brewer.pal(length(bins), "PiYG")[-med]
  colors <- colpal[findInterval(varofint, bins, rightmost.closed=T)]
  plot(shape.shp, col=colors)
  title(var.names[[i]], cex=1)
  legend("topright",fill=colpal,legend=paste(round(bins[-length(bins)],2),
                                             "-",round(bins[-1],2)),cex=0.5)
}



# Comparison
# mt.mod:  Moran'I of district level resiudals of glm model with multi-level 
# structure
# mt.mod2: Moran'I of district level residuals of multi-level model
# mt.mod_re2: Moran'I of district level random effects of multi-level model
# mt.mod3: Moran'I of district level residuals of multi-level model with 
# lagged Quality Care variable
# mt.mod_re3: Moran'I of district level random effects of multi-level model
# with lagged Quality Care variable
mt.mod
mt.mod2
mt.mod_re2
mt.mod3
mt.mod_re3
