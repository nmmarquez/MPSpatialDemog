################################################################################
#                                                                              #
# IDEM 156 - Spatial Demography Course 2018                                    #
# Spatial Lag and Spatial Error Models                                         #
# Sebastian Kl?sener, MPIDR                                                    #
#                                                                              #
################################################################################

# Erase all objects in memory
rm(list = ls(all = TRUE))

# Load libraries
library(spdep)
library(RColorBrewer)
library(rgdal)
library(maptools)
library(lrmest)
library(fBasics)
library(car)
library(perturb)

# # Set your working directory
# # If you need to change the working directory, you can use the following code:
# # Set drive
# main <- c("N")
# # Set path to session folders
# path <- ":/IMPRSD/IDEM 156 Spatial"
# # Define session folder
# path2 <- "/02 Sessions/06 Introduction to Spatial Modeling I"
# # Set working directory
# setwd(paste(main,path,path2,sep=""))


################################################################################
# 1) Import and prepare data                                                   #
################################################################################

# Open shapefile
shape.shp <- readShapePoly('2004_06', IDvar="KREIS_KENN")
shapeogr.shp <- readOGR(".", "2004_06")
proj4string(shape.shp) <- proj4string(shapeogr.shp)
shape.ll <- spTransform(shape.shp, CRS("+proj=longlat +datum=WGS84"))
plot(shape.shp)                   

names(shape.shp)
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

# Coordinates of UTM-projected shapefile
coords <- coordinates(shape.shp)
# Coordinates from reprojected Shapefile with LongLat-projection 
coords.ll <- coordinates(shape.ll)

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

# In order to "attach" data from a shapefile attribute table, we have to call it
# with "@data"
attach(shape.shp@data)

# Create two random variables
ran.y <- rnorm(n)
ran.x <- rnorm(n)

# Create neighborhood weight matrices
# First order queen contiguity
nb.FOQ <- poly2nb(shape.shp, queen=TRUE)
summary(nb.FOQ)

# Code to connect R?gen
#edit.nb(nb.FOQ, polys=shape.shp)
nb.FOQ[[350]]
nb.FOQ[[358]]
nb.FOQ[[362]]

nb.FOQ.cor <- nb.FOQ
nb.FOQ.cor[[362]]  <- as.integer(c(350, 358))
nb.FOQ.cor[[350]]  <- as.integer(c(358, 362))
nb.FOQ.cor[[358]]  <- as.integer(c(350,352,353,354,360,362))
lw.FOQ <- nb2listw(nb.FOQ.cor) 

# Create a nearest neighbor weigh matrix
l.5NN <- knearneigh(coords.ll, k=5, longlat=T)
nb.5NN <- knn2nb(l.5NN, row.names=id)
lw.5NN <- nb2listw(nb.5NN)

# Create a distance based matrix (distance between centroids)
# Distance of 50 km
d50km <- dnearneigh(coords.ll, 0, 50, row.names = id, longlat=T)
lw.50km <- nb2listw(d50km)
# Distance of 100 km
d100km <- dnearneigh(coords.ll, 0, 100, row.names = id,longlat=T)
lw.100km <- nb2listw(d100km)

################################################################################
# 2)   Create Spatial Variables                                                #
# 2.1) Distance to nearest big town                                            #
################################################################################

# Vector of names of big cities, which we want to consider in our "distance to
# nearest town variable". Instead of the name, we could also use the IDs for 
# these towns, if we know them.
bigcit <- c("Berlin","Hamburg","Koeln","Muenchen","Frankfurt am Main",
 "Stuttgart","Dresden","Leipzig", "Nuernberg","Duesseldorf","Dortmund",
 "Region Hannover","Bremen")
length.bc <- length(bigcit)

# Determine the position of big towns by first determining their row-wise
# position in the dataframe. For this we use a list-element
pos <- list()
for (i in 1:length.bc) {
    pos[[i]] <- c(which(NAMES==bigcit[i]))
}
pos

# There are two districts called M?nchen (M?nchen Stadt and Land)
# We have to detect the city of M?nchen by its ID
alt <- c(pos[[4]])
num <- which(alt==which(ID==9162))
pos[[4]] <- pos[[4]][num]

# Unlist the list into a vector
cbc <- c(unlist(pos))

# Derive the coordinates of the big towns, which are stored in the 
# coords-vector in the same position (e.g. element 9) as the information for 
# this town in the shapefile data attribute table (row 9)
coorbc.ll <- as.matrix(coords.ll[cbc,])

# Create a n*cbc matrix containing the distances between all regional centroids 
# and all biggest towns.
distances <- spDists(coords.ll, coorbc.ll,longlat=T)

dist_dataf <- data.frame(distances)
colnames(dist_dataf) <- bigcit
# ids as rownames
rownames(dist_dataf) <- id
# ids as seperate column in the dataframe
dist_dataf$ID <- id

# Check derived distances
head(dist_dataf)


# Derive row-wise the minimum distance between a region and the nearest big
# town.
DISTBC <- c(1:n)
for (i in 1:n) {
    DISTBC[i] <- min(distances[i,])
}

# Plot distances
var.data <- data.frame(DISTBC)
var.names <- c("Distance to big cities")
m <- length(var.data)
par(mfrow=c(1,1))
for (i in 1:m) {
    varofint <- var.data[,i]
    bins <- c(0,10,20,30,40,50,100,200,300)
    lb <- c(length(bins)-1)
    colpal <- brewer.pal(length(bins), "RdYlGn")
    colpal <- c(colpal[1:5],colpal[7:11])
    colors <- colpal[findInterval(varofint, bins, rightmost.closed=T)]
    plot(shape.shp, col=colors, lty=0, lwd=0.0001, bg="grey5")
    title(var.names[[i]], cex=1)
    legend("bottomright",fill=colpal,legend=paste(round(bins[-length(bins)],2),
    "-",round(bins[-1],2)),cex=0.6, bg="white")
}


################################################################################
# 2.2) Easting and Northing                                                    #
################################################################################

# Create variables Easting and Northing
EAST <- coords.ll[,1]
NORTH <- coords.ll[,2]

# Plot East and North 
var.data <- data.frame(EAST, NORTH)
var.names <- c("Easting","Northing")
m <- length(var.data)
par(mfrow=c(1,2))
for (i in 1:m) {
    varofint <- var.data[,i]
    mean <- mean(varofint)
    sd <- sd(varofint)
    min <- min(varofint)
    max <- max(varofint)
    bins <- c(min,mean-sd, mean, mean+sd, max)
    lb <- length(bins)
    med <- median(1:lb)
    colpal <- brewer.pal(lb, "PRGn")[-med]
    colors <- colpal[findInterval(varofint, bins, rightmost.closed=T)]
    plot(shape.shp, col=colors, lty=0, lwd=0.0001, bg="grey5")
    title(var.names[[i]], cex=1)
    legend("bottomright",fill=colpal,legend=paste(round(bins[-length(bins)],2),
    "-",round(bins[-1],2)),cex=0.6, bg="white")
}


################################################################################
# 2.3) Create Lagged Variables                                                 #
################################################################################

# Create lagged covariates
# Define neigborhood weight matrix to create lagged variables
listw <- nb2listw(nb.FOQ.cor) 

LOGPD05 <- log(PD05)
INCOMET05 <- INCOME05/1000
LOGTOUR05 <- log(TOUR05)

# Create lagged variables
lagvar <- data.frame(SNM05,AVAGE05,TFR05,SHB152505,INCOME05,
 UER05,PD05,SHFF05, TOUR05, FR155005,CC_u3,CC_3_6)
llv <- length(lagvar[1,])

for (i in 1:llv) {
    lagvar[,i] <- lag.listw(listw,lagvar[,i])
}
# Divide income by 1000
lagvar[,5] <- (lagvar[,5]/1000)
# Create log version of variables
lagvar[,7] <- log(lagvar[,7])
lagvar[,9] <- log(lagvar[,9])

# Define names of variables 
colnames(lagvar) <- c("LAGSNM05","LAGAVAGE05","LAGTFR05","LAGSHB152505",
 "LAGINCOMET05","LAGUER05","LAGLOGPD05","LAGSHFF05", "LAGLOGTOUR05", 
 "LAGFR155005","LAG_CC_u3","LAG_CC_3_6")
attach(lagvar)

# Plot lagged variables
lagvar1 <- data.frame(SNM05,AVAGE05,TFR05,SHB152505, INCOMET05,
 UER05,LOGPD05,SHFF05, LOGTOUR05, FR155005, CC_u3, CC_3_6)

par(mfrow=c(2,2))
for (i in 1:llv) {
    plot(lagvar1[,i],lagvar[,i], main=colnames(lagvar1[i]), 
    xlab=colnames(lagvar1[i]), 
    ylab=paste("Spatially Lagged",colnames(lagvar1[i])))
    legend("bottomright",paste("Cor:",round(cor(lagvar1[,i],lagvar[,i]),2)), 
     bty="n")  
}


################################################################################
# 3)   Preparation of Model Analysis                                           #
# 3.1) Create correlation matrix                                               #
################################################################################

# Construct correlation matrix
av.var <- data.frame(SNM05,AVAGE05,TFR05,SHB152505,INCOMET05,
 UER05,log(PD05),SHFF05,log(TOUR05),FR155005,DISTBC,CC_u3,CC_3_6,EAST,NORTH)
av.var.nam <- colnames(av.var)
len.av.var <- length(av.var[1,])

cor.mat <- data.frame(matrix(nrow=len.av.var, ncol=len.av.var))
 
for (i in 1:len.av.var) {
    for (j in 1:len.av.var) {
    cor.mat[i,j] <- round(cor(av.var[,i],av.var[,j]),3) 
    }
}
colnames(cor.mat) <- av.var.nam
rownames(cor.mat) <- av.var.nam
res.cor.mat <- cor.mat
res.cor.mat[upper.tri(res.cor.mat)] <- c(" ")
res.cor.mat

#write.table(res.cor.mat,"Correlation-matrix.csv", sep=",")


################################################################################
# 3.2) Maps of the variables of interest and explanatory variables             #
################################################################################

# Create dataframe, in which CC_u3 (column 12) East (column 14) and 
# North (column 15) are excluded
var.data <- av.var[,-c(12,14:15)]
var.names <- av.var.nam[-c(12,14:15)]

par(mfrow=c(1,1))
# This code produces m maps with standard deviation categorization
m <- length(var.data)
var.color <- matrix(nrow=n,ncol=m, data=rep(0,n*m))
for (i in 1:m) {
    varofint <- var.data[,i]
    mean <- mean(varofint)
    sd <- sd(varofint)
    min <- min(varofint)
    max <- max(varofint)
    bins <- c(min,mean-sd, mean, mean+sd, max)
    lb <- length(bins)
    med <- median(1:lb)
    colpal <- brewer.pal(lb, "PRGn")[-med]
    var.color[,i] <- colpal[findInterval(varofint, bins, rightmost.closed=T)]
    plot(shape.shp, col=var.color[,i], bg="grey5",lty=0)
    title(var.names[[i]], cex=1)
    legend("bottomright",fill=colpal,legend=paste(round(bins[-length(bins)],2),
    "-", round(bins[-1],2)),cex=0.7, bg="white")
}

var.data <- data.frame(av.var[,c(12)])
var.names <- av.var.nam[c(12)]

par(mfrow=c(1,1))
# This code produces m maps with 0.5 standard deviation categorization,
# as for variable 12 mean-sd is < 0.
m <- length(var.data)
var.color <- matrix(nrow=n,ncol=m, data=rep(0,n*m))
for (i in 1:m) {
  varofint <- var.data[,i]
  mean <- mean(varofint)
  sd <- sd(varofint)
  min <- min(varofint)
  max <- max(varofint)
  bins <- c(min,mean-(sd*0.5), mean, mean+(sd*0.5), max)
  lb <- length(bins)
  med <- median(1:lb)
  colpal <- brewer.pal(lb, "PRGn")[-med]
  var.color[,i] <- colpal[findInterval(varofint, bins, rightmost.closed=T)]
  plot(shape.shp, col=var.color[,i], bg="grey5",lty=0)
  title(var.names[[i]], cex=1)
  legend("bottomright",fill=colpal,legend=paste(round(bins[-length(bins)],2),
  "-", round(bins[-1],2)),cex=0.7, bg="white")                                              
}
                                                 
################################################################################
# 3.3) Create function to show significance levels with stars                  #
################################################################################

# Here we write a function that assigns significance stars to p-values
# Whenever we will call sig.stars(x) with x being a number, it will perform 
# this code.

sig.stars <- function(x) {
            len <- length(x)
            res <- c(1:len)       
            for (i in 1:len) {
                pv <- x[i]
                if (pv <0.001) pv <- c("***")
                if (pv >=0.001 & pv <0.01) pv <- c("**")
                if (pv >=0.01 & pv <0.05) pv <- c("*")
                if (pv >=0.05 & pv <0.1) pv <- c("'")
                if (pv >=0.1) pv <- c(" ")
                res[i] <- pv
            }    
            res
}


################################################################################
# 4)   Model Data                                                              #
# 4.1) Define Linear Models - Here you can fill in as many models as you       #
#      like with models numbered in ascending order                            #
################################################################################

# Available variables (see Word-Sheet)
# SNM05, AVAGE05, TFR05, SHB152505
# INCOMET05, UER05, LOGPD05, SHFF05, LOGTOUR05, FR155005, CC_u3, CC_3_6, 
# random.y, random.x
# 
# LAGSNM05, LAGAVAGE05, LAGTFR05, LAGSHB152505
# LAGINCOMET05, LAGUER05, LAGLOGPD05, LAGSHFF05, LAGLOGTOUR05 
# LAGFR155005, LAGCC_u3, CC_3_6
# 
# DISTBC, EAST, NORTH

# Available weight matrices
# lw.FOQ, lw.5NN, lw.50km, lw.100km

# Define models and weight matrices 
# Only objects of type lm() (linear model) accepted!
our.models           <- list()
our.weight.mats      <- list()

our.models[[1]]        <- lm(ran.y~ran.x)
our.weight.mats[[1]]   <- lw.FOQ
our.models[[2]]        <- lm(AVAGE05~INCOMET05)
our.weight.mats[[2]]   <- lw.FOQ
our.models[[3]]        <- lm(AVAGE05~INCOMET05)
our.weight.mats[[3]]   <- lw.100km
our.models[[4]]        <- lm(AVAGE05~INCOMET05)
our.weight.mats[[4]]   <- lw.5NN
#our.models[[5]]        <- 
#our.weight.mats[[5]]   <- 
#our.models[[6]]       <- 
#our.weight.mats[[6]]  <- 
#our.models[[7]]       <-
#our.weight.mats[[7]]  <- 
#our.models[[8]]       <-
#our.weight.mats[[8]]  <- 
#our.models[[9]]       <-
#our.weight.mats[[9]]  <- 
#our.models[[10]]      <-
#our.weight.mats[[10]] <- 


################################################################################
# 4.2) The following automatized code calculates linear models and tests and   #
#      fills the results up into a list-object                                 #
#      Please only change it if you know what you doing!                       #
################################################################################

# Obtain information on number of models
num.models <- length(our.models)

# Derive names of variables used as covariates in the models
namvar <- list()
for (i in 1:num.models) {
    namvar[[i]] <- c(colnames(our.models[[i]]$model)[-1])
}
namvarl <- unique(unlist(namvar))
l.namvarl <- length(namvarl)

# Models and test calculations. Please notice that in each step we are further
# filling up the list "our.models" with new listwise elements. The nice thing 
# about lists is that there are no limitations with regard to the object that 
# you can put into a list (e.g. vector, matrix, data.frame, model results 
# objects, etc.). In case we have 6 models, the first six list-elements will
# contain the model call (as defined in the preceeding section). The summary of
# the model (see below) will then be stored in the list at the position 7-12, 
# the Moran's I on the model in the positions 13-18, and so on. 

# Summary of models
for (i in 1:num.models) {
    our.models[[i+num.models]] <- summary(our.models[[i]])
}

# Calculate Moran's I test for Spatial Autocorrelation in the error term
for (i in 1:num.models) {
    our.models[[i+2*num.models]] <- lm.morantest(our.models[[i]],
     our.weight.mats[[i]])
}

# Calculate Lagrange Multiplier tests (Spatial Lag vs. Spatial Error)
for (i in 1:num.models) {
    our.models[[i+3*num.models]] <- lm.LMtests(our.models[[i]], 
     our.weight.mats[[i]], test=c("LMerr","RLMerr","LMlag","RLMlag"))
}

# Calculate Moran's I test for Dependent Variables
for (i in 1:num.models) {
    dpvar <- our.models[[i]]$model[1]
    our.models[[i+4*num.models]] <- moran.test(dpvar[,1], our.weight.mats[[i]]) 
}

# Calculate AIC
for (i in 1:num.models) {
    our.models[[i+5*num.models]] <- AIC(our.models[[i]])
}

# Checking for Multicollinearity with kappa function
for (i in 1:num.models) {
  our.models[[i+6*num.models]] <- colldiag(our.models[[i]])
}

# Calculate Jarque-Bera Test on Residuals (whether skewness and kurtosis are
# matching a normal distribution)
for (i in 1:num.models) {
  our.models[[i+7*num.models]] <- jarqueberaTest(our.models[[i]]$residuals)
}

################################################################################
# 4.3) Code for table showing results of your models                           #
################################################################################

# The following code automatically creates a result table based on the model 
# specifications which you set in section 4.1
results <- data.frame(matrix(ncol=num.models*2,nrow=3+l.namvarl))

for (i in 1:num.models) {
    colnames(results)[(i*2)-1] <- paste("Model",i)
    colnames(results)[(i*2)] <- paste(" ")
    results[1,(i*2)-1] <- paste(colnames(our.models[[i]]$model)[1])
    results[1,(i*2)] <- paste(" ")
    results[2,] <- paste(" ")
    results[3,] <- rep(c("Coefficients","p-value"),num.models)
    paste("  ",our.models[[i]]$call[2],"  ")
    varinmod <- c(colnames(our.models[[i]]$model)[-1])
    ncof <- length(varinmod)
    ord <- c(1:ncof)
    for (j in 1:ncof) {
        ord[j] <- which(namvarl == varinmod[j])+1
    }
    ordi <- c(1,ord)
    l.ordi <- length(ordi)
    for (k in 1:l.ordi) {
        results[ordi[k]+3,(i*2)-1] <- 
         round(our.models[[i+num.models]]$coefficients[k,1],3)
        results[ordi[k]+3,(i*2)] <- 
         sig.stars(our.models[[i+num.models]]$coefficients[k,4])
    }
    results[5+l.namvarl,] <- rep("------",2*num.models)
    results[6+l.namvarl,(i*2)-1] <- 
     round(our.models[[i+num.models]]$adj.r.squared,2)
    results[6+l.namvarl,(i*2)] <- c("-")
    results[7+l.namvarl,(i*2)-1] <- round(our.models[[i+5*num.models]],2)
    results[7+l.namvarl,(i*2)] <- c("-")
    results[8+l.namvarl,] <- rep("------",2*num.models) 
    results[9+l.namvarl,(i*2)-1] <- 
     paste(attributes(our.weight.mats[[i]])$call[2])
    results[9+l.namvarl,(i*2)] <- c("-")
    results[10+l.namvarl,(i*2)-1] <- 
     round(our.models[[i+4*num.models]]$estimate[1],2)
    results[10+l.namvarl,(i*2)] <- 
     sig.stars(our.models[[i+4*num.models]]$p.value)
    results[11+l.namvarl,(i*2)-1] <- 
     round(our.models[[i+2*num.models]]$estimate[1],2)
    results[11+l.namvarl,(i*2)] <- 
     sig.stars(our.models[[i+2*num.models]]$p.value)
    results[12+l.namvarl,(i*2)-1] <- 
     round(our.models[[i+3*num.models]]$LMerr$statistic,3)
    results[12+l.namvarl,(i*2)] <- 
     sig.stars(our.models[[i+3*num.models]]$LMerr$p.value)
    results[13+l.namvarl,(i*2)-1] <- 
     round(our.models[[i+3*num.models]]$RLMerr$statistic,3)
    results[13+l.namvarl,(i*2)] <- 
     sig.stars(our.models[[i+3*num.models]]$RLMerr$p.value)
    results[14+l.namvarl,(i*2)-1] <- 
     round(our.models[[i+3*num.models]]$LMlag$statistic,3)
    results[14+l.namvarl,(i*2)] <- 
     sig.stars(our.models[[i+3*num.models]]$LMlag$p.value)
    results[15+l.namvarl,(i*2)-1] <- 
     round(our.models[[i+3*num.models]]$RLMlag$statistic,3)
    results[15+l.namvarl,(i*2)] <- 
     sig.stars(our.models[[i+3*num.models]]$RLMlag$p.value)
    results[16+l.namvarl,(i*2)-1] <- 
     round(sort(our.models[[i+6*num.models]]$condindx)[length(our.models[[i+6*num.models]]$condindx)],3)
    results[16+l.namvarl,(i*2)] <- c("-")
    results[17+l.namvarl,(i*2)-1] <- 
     round(as.numeric(paste(our.models[[i+7*num.models]]@test$statistic)),3)
    results[17+l.namvarl,(i*2)] <- 
     sig.stars(as.numeric(paste(our.models[[i+7*num.models]]@test$p.value)))
}

# Change row names
row.names(results) <- c("Dependent Variable"," ","   ","Intercept",namvarl,
 "-----------","Adj.r2", "AIC","+-+-+-+-+-+","NB Weight Matrix",
 "Moran's I (Dep. Var.)", "Moran's I (Model)", "LMerr", "RLMerr","LMlag",
 "RLMlag", "Multicollinearity","Jarque Bera")

results
# Write results in file
#write.table(results,"results.csv", sep=",")


################################################################################
#                                                                              #
# 4.4) Variable and residual maps                                              #
#                                                                              #
################################################################################

# Derive data
res.data <- data.frame(matrix(ncol=num.models*2,nrow=n))
for (i in 1:num.models) {
    res.data[,(i-1)*2+1] <-  our.models[[i]]$model[1]
    res.data[,(i-1)*2+2] <-  round(our.models[[i]]$residuals,4)
}
# Generate titles and subtitles
res.names <- rep("NA",num.models*2)
sub.names <- rep("NA",num.models*2)
for (i in 1:num.models) {
    res.names[(i-1)*2+1] <-  
     paste("Dep. Var. Mod.",i,":\n",colnames(our.models[[i]]$model[1]))
    res.names[(i-1)*2+2] <-  
     paste("Residuals Mod.",i,":\n",our.models[[i]]$call[2])
#     paste("Residuals Mod.",i,":\n",colnames(our.models[[i]]$model)[-1])
}

res.names <- rep("NA",num.models*2)
sub.names <- rep("NA",num.models*2)
for (i in 1:num.models) {
    res.names[(i-1)*2+1] <-  paste("Dep. Var. Mod.",i,":")
    sub.names[(i-1)*2+1] <- colnames(our.models[[i]]$model[1])
    res.names[(i-1)*2+2] <- paste("Residuals Mod.",i,":")
    sub.names[(i-1)*2+2] <- paste(our.models[[i]]$call[2])
}

# Plot maps
par(mfrow=c(2,2))
# This code produces m maps with standard deviation categorization
m <- length(res.data)
var.color <- matrix(nrow=n,ncol=m, data=rep(0,n*m))
for (i in 1:m) {
    varofint <- res.data[,i]
    mean <- mean(varofint)
    sd <- sd(varofint)
    min <- min(varofint)
    max <- max(varofint)
    bins <- c(min,mean-sd, mean, mean+sd, max)
    lb <- length(bins)
    med <- median(1:lb)
    colpal <- brewer.pal(lb, "PRGn")[-med]
    var.color[,i] <- colpal[findInterval(varofint, bins, rightmost.closed=T)]
    plot(shape.shp, col=var.color[,i], bg="grey5",lty=0)
    title(res.names[[i]],sub=sub.names[[i]], cex.main=1,cex.sub=1)
    legend("bottomright",fill=colpal,legend=paste(round(bins[-length(bins)],2),
    "-", round(bins[-1],2)),cex=0.7, bg="white")
}

# Your task: Choose one of the available demographic covariates (SNM05, AVAGE05,
# TFR05, SHB152505) and try to construct a meaningful model. Use test-statistics
# to evaluate the performance and try different neighborhood definitions. 
# You might face problem with multicollinearity (see correlation matrix), and 
# you might be at risk of running into ecological fallacies (e.g. with the 
# variable Share Foreign Females). Determine, whether it is necessary to use a 
# a spatial lag or spatial error model, or whether it can be avoided.


################################################################################
# 5)   LM vs. Spatial Lag vs. Spatial Error Model                              #
# 5.1) Define Model                                                            # 
################################################################################

# In case your test statistics indicate that substantial spatial 
# auto-correlation remains in every meaningful lm-Model specifiation, you can
# run a spatial lag and/ or spatial error model
# e.g AVAGE05~UER05+LOGPD05
modeltext <- formula(AVAGE05~UER05+LOGPD05)

# Define neighborhood weight matrix
# e.g. lw.FOQ
listwise <- lw.FOQ


################################################################################
# 5.2) Calculate Model                                                         # 
################################################################################

# Calculate three models - lm, lagmodel and errormodel
# Linear model
lmmodel <- lm(modeltext)

# Spatial lag model
lagmodel <- lagsarlm(modeltext,listw=listwise)

# Spatial Error Model
errormodel <- errorsarlm(modeltext,listw=listwise)

# Summaries of these three models
summary(lmmodel)
summary(lagmodel)
summary(errormodel)

# Spatial Durbin Model
durbinmodel <- lagsarlm(modeltext,listw=listwise,type="Durbin")
summary(durbinmodel)

################################################################################
# 5.3) Calculate Moran's I on model/ model residuals                           # 
################################################################################

# Test for spatial auto-correlation of the residuals
lm.morantest(lm(modeltext),listwise) 
summary(lagmodel)
moran.test(errormodel$residuals,listwise) 
