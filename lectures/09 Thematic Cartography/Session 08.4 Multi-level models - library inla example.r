################################################################################
#                                                                              #
# IDEM 156 Spatial Demography Course 2018                                      #                          
# Multi Level Models and Spatial Autocorrelation - INLA                        #
# Sebastian Kl?sener, MPIDR                                                    #
#                                                                              #
################################################################################

# Erase all objects in memory
rm(list = ls(all = TRUE))

# Opening libraries
library(spdep)
library(RColorBrewer)
library(INLA)
library(maptools)
library(rgdal)
library(lme4)


################################################################################
# 1) Import simulated data                                                     #
#    Ind.-level:                                                               #
#    Registered Severe Heart Attacks                                           #
#    - Did Person die within the first 24 hours after the attack               #
#                                                                              #
################################################################################

setwd("~/Documents/Classes/MPSpatialDemog/lectures/08 Spatial Multi-level Modeling/")
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


# Create a neighborhood file
nb.FOQ <- poly2nb(shape.shp, queen=TRUE)

# Create an INLA Weitght matrix
lw.FOQ <- nb2INLA("FOQ_INLA",nb.FOQ)


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

y <- Death
x1 <- Qual_Emer_Care


################################################################################
# 3.1) Simple logistic model, in which we do not take account for the multi-   #
#      level-structure of our data                                             #
################################################################################

# Simple logistic model using glm
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

an <- data.frame(1:length(shape.shp),shape.shp$KREIS_KENN.1)
o <- match(Ind.data$vid,an$shape.shp.KREIS_KENN.1)
ID <- an$X1.length.shape.shp.[o]

ID <- data.frame(District=as.character(row.names(shape.shp@data)), 
           ID=1:nrow(shape.shp@data)) %>% 
    right_join(data.frame(District=as.character(District)), "District") %>% 
    .$ID

formula=Death~Qual_Emer_Care + f(ID, model="iid")

result <- inla(formula,data=data.frame(Death,Qual_Emer_Care,ID,District),family="binomial")
summary(result)

ML.Model <- glmer(formula=Death~Qual_Emer_Care+(1|District), family = binomial)          
summary(ML.Model)


# Create a neighborhood file
nb.FOQ <- poly2nb(shape.shp, queen=TRUE)

# Create an INLA Weitght matrix
lw.FOQ <- nb2INLA("FOQ_INLA", nb.FOQ)

# Calculating as a Besag Model (Spatial Random Effects)
formula=Death~Qual_Emer_Care+f(ID,model="besag",graph="FOQ_INLA")
result1 <- inla(formula,data=data.frame(ID,Death,Qual_Emer_Care),family="binomial")
summary(result1)
