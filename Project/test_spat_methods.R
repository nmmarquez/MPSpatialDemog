rm(list=ls())

library(spdep)
library(rgdal)
library(maptools)
library(INLA)
library(ggplot2)
library(dplyr)
library(surveillance)
library(spgwr)
library(gridExtra)

setwd("~/Documents/Classes/MPSpatialDemog/lectures/07 Introduction to Spatial Modeling II/")

base_plot <- function(datur, var_, continuous=F, inverse=F){
    datur$QVAR <- cut_interval(datur[,var_], 8)
    datur$VAR <- datur[,var_]
    inv <- ifelse(inverse, 1, -1)
    bmap <- tidymap %>% left_join(datur) %>%
        ggplot(aes(x=long, y=lat)) +
        theme_classic()
    if(!continuous){
        bmap <- bmap + geom_polygon(aes(group=group, fill = QVAR)) +
            scale_fill_brewer(palette = "Spectral", direction=inv) +
            scale_color_brewer(palette = "Spectral", direction=inv) + 
            theme(axis.line = element_blank(),
                  legend.title=element_blank(), axis.text=element_blank(),
                  axis.ticks=element_blank(), axis.title=element_blank())
    }
    else{
        bmap <- bmap + 
            geom_polygon(aes(group=group, fill = VAR)) +
            scale_fill_distiller(palette = "Spectral", direction=inv) + 
            theme(axis.line = element_blank(),
                  legend.title=element_blank(), axis.text=element_blank(),
                  axis.ticks=element_blank(), axis.title=element_blank())
    }
    return(bmap + geom_path(aes(group=group), size=.1))
}

# Open shapefile
shape.shp <- readShapePoly('2004_06', IDvar="KREIS_KENN")
shape.shp@data$id <- row.names(shape.shp@data)
tidymap <- fortify(shape.shp)
shapeogr.shp <- readOGR(".", "2004_06")
proj4string(shape.shp) <- proj4string(shapeogr.shp)
shape.shp <- spTransform(shape.shp, CRS("+proj=longlat +datum=WGS84"))

# Derive information on internal id and centroids. In this case we turn the id-
# values first into a numeric vector, which is easier to use, if we want to
# match data to the shapefile.
id <- as.numeric(paste(shape.shp$KREIS_KENN))/1000

# Define some standard objects
n <- length(id)
sq <- c(1:n)
sq_text <- paste(sq)
coords <- coordinates(shape.shp)

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


# Create neighborhood weight matrices
# First order queen contiguity
nb.FOQ <- poly2nb(shape.shp, queen=TRUE)
# Create Adjacency map
smat <- poly2adjmat(shape.shp)
# Create spatial variable indicators
shape.shp@data$ID <- 1:nrow(shape.shp@data)
shape.shp@data$ID2 <- 1:nrow(shape.shp@data)

# run the geography weighted model with no dummy variable
gwmodel <- with(
    shape.shp@data, 
    gwr.sel(SNM05 ~ FR155005, coords=coords, longlat=T) %>%
        gwr(SNM05 ~ FR155005, coords=coords, bandwidth=., longlat=T))

# run the in lamodel with no dummy variable
f_ <- SNM05 ~ FR155005 + 
    f(ID, model="besag", graph=inla.read.graph(smat), 
      hyper=list(prec=list(prior="loggamma",param=c(1, 10)))) +
    f(ID2, FR155005, model="besag", graph=inla.read.graph(smat) , 
      hyper=list(prec=list(prior="loggamma",param=c(1, 10))))  
    
remodel <- inla(f_, data=shape.shp@data, family="Gaussian") 

# plot the results
grid.arrange(
    gwmodel$SDF@data %>% mutate(id=shape.shp@data$id) %>% 
        base_plot("(Intercept)") + labs(title="GWR Intercept"),
    
    remodel$summary.random$ID %>% mutate(id=shape.shp@data$id) %>% 
        base_plot("mean") + labs(title="RE Intercept"),
    
    gwmodel$SDF@data %>% mutate(id=shape.shp@data$id) %>% 
        base_plot("FR155005") + labs(title="GWR Employment Effects"),
    
    remodel$summary.random$ID2 %>% mutate(id=shape.shp@data$id) %>% 
        base_plot("mean") + labs(title="RE Employment Effects"),
    nrow=2
)

# run the geography weighted and INLA models with the new east west dummy variable
gwmodelDum <- with(
    shape.shp@data, 
    gwr.sel(SNM05 ~ FR155005 + DUM_E, coords=coords, longlat=T) %>%
        gwr(SNM05 ~ FR155005 + DUM_E, coords=coords, bandwidth=., longlat=T))


f_dum <- SNM05 ~ FR155005 + DUM_E +
    f(ID, model="besag", graph=smat, 
      hyper=list(prec=list(prior="loggamma",param=c(1, 10)))) +
    f(ID2, FR155005, model="besag", graph=smat, 
      hyper=list(prec=list(prior="loggamma",param=c(1, 10))))

remodel_dum <- inla(f_dum, data=shape.shp@data, family="Gaussian") 

grid.arrange(
    gwmodelDum$SDF@data %>% mutate(id=shape.shp@data$id) %>% 
        base_plot("(Intercept)") + labs(title="GWRDUM Intercept"),
    
    remodel_dum$summary.random$ID %>% mutate(id=shape.shp@data$id) %>% 
        base_plot("mean") + labs(title="REDUM Intercept"),
    
    gwmodelDum$SDF@data %>% mutate(id=shape.shp@data$id) %>% 
        base_plot("FR155005") + labs(title="GWRDUM Employment Effects"),
    
    remodel_dum$summary.random$ID2 %>% mutate(id=shape.shp@data$id) %>% 
        base_plot("mean") + labs(title="REDUM Employment Effects"),
    nrow=2
)

# plot the correlation between betas in the GWR model 
gwmodel$SDF@data %>% mutate(id=shape.shp@data$id) %>% 
    left_join(shape.shp@data %>% mutate(FR155005=NULL)) %>% 
    ggplot(aes(x=`(Intercept)`,y=FR155005, color=as.factor(DUM_E))) + 
    geom_point()
