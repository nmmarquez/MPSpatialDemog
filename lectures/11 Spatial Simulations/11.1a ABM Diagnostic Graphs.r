################################################################################
#                                                                              #  
# ABM model of the Swedish fertility decline                                   #
# Diagnostics                                                                  #
# Social adaptation algorithm                                                  #
# Sebastian Klüsener, MPIDR                                                    #
#                                                                              #
################################################################################

# Erase all objects in memory
rm(list = ls(all = TRUE))

# Load libraries
library(data.table)
library(spacom)
library(spdep)
library(maptools)
library(rgdal)
library(rgeos)
library(RColorBrewer)

# Set drive
main <- c("C")
# Set path to session folders
path <- ":/ownCloud/01 MPI/120 Spatial Demography 2018"
# Define session folder
path2 <- "/02 Sessions/11 Spatial Simulations"
# Set working directory
setwd(paste(main,path,path2,sep=""))


################################################################################
#                                                                              #                         
# 1. Load and prepare data                                                     #                         
#                                                                              #                         
################################################################################

# Swedish county shapefile
shapemer.shp <- readShapePoly('Sweden_tc_laen_1880_1900',IDvar="SP_ID")
# Projection information has to be added manually
proj4string(shapemer.shp) <- CRS("+proj=lcc +towgs84=0,0,0 +x_0=4000000 +y_0=2800000 +lon_0=10 +lat_0=52 +lat_1=35 +lat_2=65 +units=m")

shapemer.proj <- spTransform(shapemer.shp, CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"))
shapemer.ll <-spTransform(shapemer.proj, CRS("+proj=longlat"))
coords.ll <- coordinates(shapemer.ll)

# k-nearest neighbors
l.4NN <- knearneigh(coords.ll, k=4,longlat=T)
nb.4NN <- knn2nb(l.4NN)
listw <- nb2listw(nb.4NN)


################################################################################
#                                                                              #                         
# 1. Load and prepare data                                                     #                         
#                                                                              #                         
################################################################################

# Scenario
scen <- c("BC_100")

# Adaptation mechanism
adapt <- c("soc_adapt")
adapt_desc <- c("Social adaptation")

# Diffusion pathways
adapt_from <- c("same_ses_ror_rob_vgses_ror")


# Additional specifications
add <- c("")

# Parameter settings
wr_all <- c(0.7)
ws_all <- c(0.05)
x_all <- c(10)

parameter_grid <- expand.grid(wr_all,ws_all,x_all)

nsim_all <- 75

agg_list <- list()


################################################################################
#                                                                              #                         
# 2. Load files with the results and derive mean values from the different     #
#    iterations of the model                                                   #                         
#                                                                              #                         
################################################################################

for (i in 1:length(parameter_grid[,1])) {
  # Influence of share adopted in place of residence and place of birth
  wr <- parameter_grid[i,1]
  
  # Influence of vanguard group in place of residence
  ws <- parameter_grid[i,2]
  
  # Maximum adaptation rate
  x <- parameter_grid[i,3]
  
  # Set directory
  path3 <- c("/Model results")

  setwd(paste(main,path,path2,path3,sep=""))
  
  nsim <- nsim_all 
  res <- list()
  for (j in 1:nsim) {
    res[[j]] <- read.table(paste("ABM_",scen,"_",adapt,"_",adapt_from,add,"_wr_",wr,"_ws_",ws,
                                 "_x_",x,"_sim_",j,".csv",sep=""),sep=",",head=T)   
  }
  res_to_agg <- data.frame(rbindlist(res))
  
  mean_count_reg <- aggregate(res_to_agg[,c(12:length(res_to_agg))],
                              by=list(res_to_agg$regkey,res_to_agg$ses),FUN=mean)
  
  agg_list[[i]] <- data.frame(res[[1]][,c(1:12)],mean_count_reg[,-c(1:3)])
}    

toexport <- rbindlist(agg_list)         


write.table(toexport,paste("Agg_",scen,"_",adapt,"_",
                           adapt_from,add,".csv",sep=""),sep=",",row.names=F)  


toplot <- read.table(paste("Agg_",scen,"_",adapt,"_",
                           adapt_from,add,".csv",sep=""),sep=",",head=T)    


################################################################################
#                                                                              #                         
# 3. Derive diagnosticts                                                       #                         
#                                                                              #                         
################################################################################

splitvec <- paste(toplot$wr,toplot$ws,toplot$x)    

res_split <- split(toplot,splitvec)
l_res_sp <- length(res_split)

# Preparations
t0 <- which(colnames(res_split[[1]])=="X0")
tmax <- length(colnames(res[[1]]))
n_tim <- tmax-t0
identical(paste(res[[1]]$scenario[1]),scen)
if (scen=="BC_100") scenname <- "Diffusion from big cities"

setwd(paste(main,path,path2,sep=""))

obs <- read.table("151012_observed_by_ses_and_region.csv",sep=",",head=T) 
regnames <- read.table("151012_region_names.csv", sep=",",head=T)
identical(rep(regnames$id,3),res[[1]]$regkey)
agses <- sort(rep(1:3,25))
agreg <- rep(regnames$reg_agg4,3)

regagg <- paste(agses,agreg,sep="_")

# List to store the results
perc_reg <- list()
ses_mean_of_reg <- list()
ses_mean_SWE <- list()
ses_sd <- list()
ses_mean_reg <- list() 
moran_test <- list()    
cor <- list()
cortot <- list()

# Iterations for the Moran's I
nsim_m <- 99

for (i in 1:l_res_sp) {
  for_diag <- res_split[[i]]
  perc_reg[[i]] <- data.frame(for_diag[,c(1:(t0-1))],(for_diag[,c(t0:tmax)]/for_diag$n*100))
  
  # Mean of regions
  mean_reg <- aggregate(perc_reg[[i]][,c(t0:tmax)],by=list(for_diag$ses,for_diag$ses_name),FUN=mean)
  ses_mean_of_reg[[i]] <- data.frame(for_diag[c(1:3),c(1:6)],mean_reg)
  
  # Swedish mean
  ses_count_SWE <- aggregate(for_diag[,c(t0:tmax)],by=list(for_diag$ses,for_diag$ses_name),FUN=sum)
  tot_count_SWE <- aggregate(for_diag$n,by=list(for_diag$ses,for_diag$ses_name),FUN=sum)
  ses_count_SWE[3:length(ses_count_SWE)] <- ses_count_SWE[3:length(ses_count_SWE)]/tot_count_SWE$x*100
  ses_mean_SWE[[i]] <- data.frame(for_diag[c(1:3),c(1:6)],ses_count_SWE)
  
  # Standard deviation of regions
  sd_reg <- aggregate(perc_reg[[i]][,c(t0:tmax)],by=list(for_diag$ses,for_diag$ses_name),FUN=sd)
  ses_sd[[i]] <- data.frame(for_diag[c(1:3),c(1:6)],sd_reg) 
 
  # Mean of regions
  adcount_reg <- aggregate(for_diag[,c(t0:tmax)],by=list(regagg),FUN=sum)
  totcount_reg <- aggregate(for_diag$n,by=list(regagg),FUN=sum)
  adcount_reg[2:length(adcount_reg)] <- adcount_reg[2:length(adcount_reg)]/totcount_reg$x*100
  ses_mean_reg[[i]] <- data.frame(for_diag[c(1:9),c(1:6)],adcount_reg) 
  moran_test[[i]] <- ses_sd[[i]]
  moran_test[[i]][,c(9:length(moran_test[[i]]))] <-0 
  for (j in (t0-4):(tmax-4)) {
    moran_test[[i]][1,j] <- moran.mc(perc_reg[[i]][1:25,(j+4)],listw,nsim_m)$statistic[1]
    #if(mean(perc_reg[[i]][1:25,(j+4)])>99) moran_test[[i]][1,j] <- NA
    moran_test[[i]][2,j] <- moran.mc(perc_reg[[i]][26:50,(j+4)],listw,nsim_m)$statistic[1]
    #if(mean(perc_reg[[i]][26:50,(j+4)])>99) moran_test[[i]][2,j] <- NA
    moran_test[[i]][3,j] <- moran.mc(perc_reg[[i]][51:75,(j+4)],listw,nsim_m)$statistic[1]
    #if(mean(perc_reg[[i]][51:75,(j+4)])>99) moran_test[[i]][3,j] <- NA
  }
  cor[[i]] <- ses_sd[[i]]
  cortot[[i]] <- ses_sd[[i]][1,]
  cortot[[i]][,8] <- paste(cortot[[i]][,8])
  cortot[[i]][1,8] <- c("Total")
  cor[[i]][,c(9:length(cor[[i]]))] <-0
  cortot[[i]][,c(9:length(cor[[i]]))] <-0
  for (j in (t0-4):(tmax-4)) {
    cor[[i]][1,j] <- cor(perc_reg[[i]][1:25,(j+4)],obs$dec9000p[1:25]*-1)
    # Stockholm city excluded
    cor[[i]][2,j] <- cor(perc_reg[[i]][27:50,(j+4)],obs$dec9000p[27:50]*-1)
    cor[[i]][3,j] <- cor(perc_reg[[i]][51:75,(j+4)],obs$dec9000p[51:75]*-1)
    cortot[[i]][1,j] <- cor(perc_reg[[i]][c(1:25,27:75),(j+4)],obs$dec9000p[c(1:25,27:75)]*-1)
  }
}    

################################################################################
#                                                                              #                         
# 4. Plot model outcomes                                                       #                         
#                                                                              #                         
################################################################################

res_ses_mean_SWE <- data.frame(rbindlist(ses_mean_SWE))
res_ses_sd <- data.frame(rbindlist(ses_sd))
res_moran <- data.frame(rbindlist(moran_test))
res_ses_mean_reg <- data.frame(rbindlist(ses_mean_reg))
res_cor <- data.frame(rbind(rbindlist(cor),rbindlist(cortot)))

plot_wr <- 0.7
plot_ws <- 0.05
plot_x <- 10
l_plo <- length(plot_x)
matchvec_pl <- paste(plot_wr,plot_ws,plot_x,sep="_")
matchvec_dat <- paste(res_ses_mean_SWE$wr,res_ses_mean_SWE$ws,res_ses_mean_SWE$x,sep="_")
matchvec_dat1 <- paste(res_cor$wr,res_cor$ws,res_cor$x,sep="_")
matchvec_dat2 <- paste(res_ses_mean_reg$wr,res_ses_mean_reg$ws,res_ses_mean_reg$x,sep="_")

cho <- match(matchvec_dat,matchvec_pl) 
cho1 <- match(matchvec_dat1,matchvec_pl) 
cho2 <- match(matchvec_dat2,matchvec_pl) 

tsi <- 1
lwdc <- 1.2
firpos <- which(colnames(res_ses_mean_SWE)=="X0")
endpos <- which(colnames(res_ses_mean_SWE)=="X200")

col_elite <- brewer.pal(9, "YlOrRd")[6]
col_worker <- brewer.pal(9, "GnBu")[8]
col_farmer <- brewer.pal(9, "YlGn")[6]
col_total <- c("grey25")

# Mean
#tiff(file=(paste("SES_mean_of_SWE_",scen,"_",plot_x[1],"_bp.tif",sep="")), width=4800, height=6400,compression="lzw",res=600)
par(mfrow=c(1,1), mar=c(4,4.5,3,1),mgp=c(2.4,1,0))
for (i in 1:l_plo) {
  dat_to_plot <- res_ses_mean_SWE[which(cho==i),]
  plot(as.numeric(paste(res_ses_mean_SWE[1,firpos:endpos])),type="l", col="transparent",lwd=lwdc,las=1,
       ylab="Share adopted\nSweden", xlab="Time period",cex.axis=tsi,cex.lab=tsi,
       ylim=c(0,100),xlim=c(0,150),xaxp=c(0,200,4),
       main=paste("wr:",plot_wr[[i]],"ws:",plot_ws[[i]],"x:",plot_x[[i]]))  
  abline(v=50,col="grey75")
  abline(v=100,col="grey75")
  abline(v=150,col="grey75")
  abline(v=200,col="grey75")
  lines(as.numeric(paste(dat_to_plot[3,firpos:endpos])),col=col_worker,lty=4,lwd=lwdc)
  lines(as.numeric(paste(dat_to_plot[2,firpos:endpos])),col=col_farmer,lty=3,lwd=lwdc)  
  lines(as.numeric(paste(dat_to_plot[1,firpos:endpos])),col=col_elite,lty=1,lwd=lwdc)
  legend("bottomright",c("Elite","Workers","Farmers"),col=c(col_elite,col_worker,col_farmer),lty=c(1,4,3),lwd=lwdc,cex=tsi,bg="white")    
}
#dev.off()

# Standard deviation

#tiff(file=(paste("SES_sd_SWE_",scen,"_",plot_x[1],"_bp.tif",sep="")), width=4800, height=6400,compression="lzw",res=600)
par(mfrow=c(1,1), mar=c(4,4.5,3,1),mgp=c(2.4,1,0))
for (i in 1:l_plo) {
  dat_to_plot <- res_ses_sd[which(cho==i),]
  plot(as.numeric(paste(res_ses_sd[1,firpos:endpos])),type="l", col="transparent",lwd=lwdc,las=1,
       ylab="Standard deviation of share adopted\n25 Swedish regions", xlab="Time period",cex.axis=tsi,cex.lab=tsi,
       ylim=c(0,40),xlim=c(0,150),xaxp=c(0,200,4),
       main=paste("wr:",plot_wr[[i]],"ws:",plot_ws[[i]],"x:",plot_x[[i]]))  
  abline(v=50,col="grey75")
  abline(v=100,col="grey75")
  abline(v=150,col="grey75")
  abline(v=200,col="grey75")
  lines(as.numeric(paste(dat_to_plot[3,firpos:endpos])),col=col_worker,lty=4,lwd=lwdc)
  lines(as.numeric(paste(dat_to_plot[2,firpos:endpos])),col=col_farmer,lty=3,lwd=lwdc)  
  lines(as.numeric(paste(dat_to_plot[1,firpos:endpos])),col=col_elite,lty=1,lwd=lwdc)  
  legend("topright",c("Elite","Workers","Farmers"),col=c(col_elite,col_worker,col_farmer),lty=c(1,4,3),lwd=lwdc,cex=tsi,bg="white")
  
}
#dev.off()

# Correlation coefficient

#tiff(file=(paste("09_SES_cor_of_SWE_",scen,"_",plot_x[1],"_bp_alt.tif",sep="")), width=4800, height=6400,compression="lzw",res=600)
par(mfrow=c(1,1), mar=c(4,4.5,3,1),mgp=c(2.4,1,0))
for (i in 1:l_plo) {
  dat_to_plot <- res_cor[which(cho1==i),]
  plot(as.numeric(paste(res_cor[1,firpos:endpos])),type="l", col="transparent",lwd=lwdc,las=1,
       ylab="Share adopted\nSweden", xlab="Time period",cex.axis=tsi,cex.lab=tsi,ylim=c(0,1),xlim=c(0,200),xaxp=c(0,200,4),
       main=paste("wr:",plot_wr[[i]],"ws:",plot_ws[[i]],"x:",plot_x[[i]]))  
  abline(v=50,col="grey75")
  abline(v=100,col="grey75")
  abline(v=150,col="grey75")
  abline(v=200,col="grey75")
  lines(as.numeric(paste(dat_to_plot[3,firpos:endpos])),col=col_worker,lty=4,lwd=lwdc)
  lines(as.numeric(paste(dat_to_plot[2,firpos:endpos])),col=col_farmer,lty=3,lwd=lwdc)  
  lines(as.numeric(paste(dat_to_plot[1,firpos:endpos])),col=col_elite,lty=1,lwd=lwdc)
  lines(as.numeric(paste(dat_to_plot[4,firpos:endpos])),col=col_total,lty=2,lwd=lwdc)
  legend("topright",c("Elite","Workers","Farmers","Total"),col=c(col_elite,col_worker,col_farmer,col_total),
         lty=c(1,4,3,2),lwd=lwdc,cex=tsi,bg="white",ncol=2)    
}
#dev.off()

# Moran's I

#tiff(file=(paste("SES_Moran_SWE_",scen,"_",plot_x[1],"_4NN.tif",sep="")), width=4800, height=6400,compression="lzw",res=600)
par(mfrow=c(1,1), mar=c(4,4.5,3,1),mgp=c(2.4,1,0))
for (i in 1:l_plo) {
  dat_to_plot <- res_moran[which(cho==i),]
  plot(as.numeric(paste(res_moran[1,firpos:endpos])),type="l", col="transparent",lwd=lwdc,las=1,
       ylab="Moran's (Spatial Clustering of Decline)\n25 Swedish regions", xlab="Time period",cex.axis=tsi,cex.lab=tsi,ylim=c(-0.2,0.6),
       main=paste("wr:",plot_wr[[i]],"ws:",plot_ws[[i]],"x:",plot_x[[i]]))  
  abline(v=50,col="grey75")
  abline(v=100,col="grey75")
  abline(v=150,col="grey75")
  abline(v=200,col="grey75")
  lines(as.numeric(paste(dat_to_plot[1,firpos:endpos])),col=col_elite,lty=2,lwd=lwdc)
  lines(as.numeric(paste(dat_to_plot[3,firpos:endpos])),col=col_worker,lty=4,lwd=lwdc)
  lines(as.numeric(paste(dat_to_plot[2,firpos:endpos])),col=col_farmer,lty=1,lwd=lwdc)  
  legend("topleft",c("Elite","Workers","Farmers"),col=c(col_elite,col_worker,col_farmer),lty=c(2,4,1),lwd=lwdc,cex=tsi,bg="white") 
}
#dev.off()
