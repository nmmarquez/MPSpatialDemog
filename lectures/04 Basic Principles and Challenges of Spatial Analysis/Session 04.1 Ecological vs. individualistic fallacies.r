################################################################################
#                                                                              #                          
# IDEM 156 - Spatial Demography Course 2018                                    #                          
# Ecological vs. individualistic fallacies                                     #                          
# Sebastian Kl?sener, MPIDR                                                    #                          
#                                                                              #                          
################################################################################


# Erase all objects in workspace
rm(list=ls(all=TRUE))

# Load libraries
library(RColorBrewer)
library(lme4)

# Set your working directory

################################################################################
# 1.)    Assumption: Awareness on how to live a healthy life should            #  
#        be higher among pupils, whose parents have a higher income            #
# 1.1)   Importing and preparing data                                          #
################################################################################

# This data is simulated data
mydata <- read.table(file="fallacies.csv", sep=",", header=T)

# Check, which variables are available
names(mydata)
# SecLev <- School ID
# HealthSc <- Health awareness of an individual pupil
# IncPar <- Income of parents of an individual pupil
# TeachSc: Score of school in teaching health awareness    


# Check first six lines of dataset
head(mydata)

# Attach dataset, so that we can call variables directly without mentioning the 
# object: E.g. instead of mydata$IncPar
attach(mydata)

# Here we group the individuals according to the school which they are 
# attending. In total there are 20 schools and in each school 20 pupils
# have been interviewed

# Define grouping variable (School-ID)
group <- SecLev

# Creation of a listwise object, which lists the data the pupils separately,
# depending on the school they are attending
SL <- split(mydata, group)
is.list(SL)

# Derive number of schools
lg <- length(SL)

# Create aggregated dataset (ecological analysis)
# Option a) Using the aggregate function
coltoagg <- which(colnames(mydata)%in%c("HealthSc","IncPar","TeachSc"))
Eco <- aggregate(mydata[coltoagg],by=list(SecLev),FUN="mean")
colnames(Eco) <- c("School.ID","E.HealthSc","E.IncPar","E.TeachSc")

# Option b) Using a loop
Eco1 <- matrix(nrow=lg,ncol=4, data=rep(0,lg*4))
Eco1[,1] <- as.numeric(names(SL))

for (i in 1:lg) {
    Eco1[i,2] <- mean(SL[[i]]$HealthSc)
    Eco1[i,3] <- mean(SL[[i]]$IncPar)
    Eco1[i,4] <- mean(SL[[i]]$TeachSc)
}
Eco1 <- data.frame(Eco1)
colnames(Eco1) <- c("School.ID","E.HealthSc","E.IncPar","E.TeachSc")

# Another example of a false negative
identical(Eco,Eco1)
Eco-Eco1

attach(Eco)


################################################################################
# 2)     Plot individual and aggregated data                                   #
#        Awareness of a school kid on health-related issues                    #
################################################################################

# Check dependent variables
par(mfrow=c(1,2),mar=rep(4,4))
hist(HealthSc) 
qqnorm(HealthSc)
qqline(HealthSc)

par(mfrow=c(1,2),mar=rep(4,4))
# Plot individual level data
plot(IncPar, HealthSc, xlab="Income", ylab="Health Awareness Score")
abline(lm(HealthSc~IncPar))
summary(lm(HealthSc~IncPar))

# Plot aggregated data
plot(Eco$E.IncPar,Eco$E.HealthSc, xlab="Income", ylab="Health Awareness Score")
abline(lm(Eco$E.HealthSc~Eco$E.IncPar))
summary(lm(Eco$E.HealthSc~Eco$E.IncPar))

# Plot multi-level
mypalette  <- rep(brewer.pal(5,"Set1"),4)
par(mfrow=c(1,1))
plot(IncPar, HealthSc, xlab="Income", ylab="Health Awareness Score", col="grey")

# Plot values for four schools seperatly
for(i in c(1,4,8,10)) {
points(SL[[i]]$IncPar,SL[[i]]$HealthSc,col=mypalette[i])
abline(lm(SL[[i]]$HealthSc~SL[[i]]$IncPar),col=mypalette[i])
}

################################################################################
# 3)     Models                                                                #
################################################################################

# Individual level model controlling for income of parents
summary(lm(HealthSc~IncPar))

# Ecological model controlling for income of parents
summary(lm(E.HealthSc~E.IncPar))

# Ecological model performs worse, as long as we do not control for the school-
# effect

# Individual level controlling for income and level of health education
# at school
summary(lm(HealthSc~IncPar+TeachSc))

# As a multi-level model
summary(lmer(formula=HealthSc~IncPar+(1|TeachSc)))

# Ecological model controlling for income and level of health education
# at school
summary(lm(E.HealthSc~E.IncPar+E.TeachSc))

# In this example, both the individual-level as well as the ecological model 
# with aggregate-level data perform quite well.
