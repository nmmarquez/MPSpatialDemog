################################################################################
#                                                                              #  
# ABM model of the Swedish fertility decline                                   #
# Diffusion from big cities scenario                                           #
# Social adaptation algorithm                                                  #
# Sebastian Klüsener, MPIDR                                                    #
#                                                                              #
################################################################################

# Erase all objects in memory
rm(list = ls(all = TRUE))

# Set drive
main <- c("N")
# Set path to session folders
path <- ":/IMPRSD/IDEM 156 Spatial"
# Define session folder
path2 <- "/02 Sessions/11 Spatial Simulations"
# Set working directory
setwd(paste(main,path,path2,sep=""))


################################################################################
#                                                                              #                         
# 1. Load and prepare data                                                     #                         
#                                                                              #                         
################################################################################

# Load agentset
women1880 <-  read.table("151012_women_1880_sorted_100.csv",
                         sep=",",head=T)
sample_sparse <- 
  read.table("151012_big_city_sample_sparse_strat_reg_ses_100.csv",
             sep=",",head=T)
  
# Total number of women
n <- length(women1880[,1])

# Prepare dataset
# Code SES-variable
ses_code <- rep(0,length(women1880[,1]))
ses_code[women1880$ses=="Elite"] <- 1
ses_code[women1880$ses=="Farmer"] <- 2
ses_code[women1880$ses=="Worker"] <- 3
ses_name <- c("Elite","Farmer","Worker")

# Code region of residence and region of birth
ror_code <- sprintf("%02.0f",women1880$laen)
rob_code <- sprintf("%02.0f",women1880$bpllaen)

# Adapted
adapt <- women1880$ad

# Age group
age_code <- women1880$agegr 

# Dataset for simulation
women1880sim <- data.frame(ses_code,ror_code,rob_code,adapt)

# Prepare vectors to derive share of adapters by SES-groups, region of
# residence (ror) and/or age group
agg_ses_ror     <- paste(ses_code,ror_code,sep="_")
agg_ses         <- paste(ses_code,sep="_")
agg_ror         <- paste(ror_code,sep="_")
agg_ses_rob     <- paste(ses_code,rob_code,sep="_")
agg_ses         <- paste(ses_code,sep="_")
agg_rob         <- paste(rob_code,sep="_")

# Load regionnames
regnames <- read.table("151012_region_names.csv", sep=",",head=T)
regcode <- sprintf("%02.0f",regnames$id)
ord <- match(sort(unique(ror_code)),regcode)
regnames_ord <- regnames[ord,]


################################################################################
#                                                                              #                         
# 2. Set specifications of ABM                                                 #                         
#                                                                              #                         
################################################################################

# Scenario
scen <- c("BC_100")

# Adaptation algorithm
ad_alg <-     c("soc_adapt")
adapt_from <- c("same_ses_ror_rob_vgses_ror")

# Number of simulations 
#n_sim <- 100
n_sim <- 3

# Number of time periods in each simulation
n_time <- 200

# Specify parameter combinations to be considered
#wr_all <- c(0.5,0.7,0.9)
#ws_all <- c(0.05,0.10)
#x_all <- c(10,15,20)

wr_all <- c(0.7)
ws_all <- c(0.05)
x_all <- c(10)

parameter_grid <- expand.grid(wr_all,ws_all,x_all)

add <- c("")

for (h in 1:length(parameter_grid[,1])) {
# Influence of share adopted in place of residence and place of birth
    wr <- parameter_grid[h,1]

# Influence of vanguard group in place of residence
    ws <- parameter_grid[h,2]

# Maximum adaptation rate
    x <- parameter_grid[h,3]

# Barriers of diffusion (e.g. ses, region, age)
    group <- agg_ses_ror
    group.b <- agg_ses_rob
    group.reg <- agg_ror
    group.ses <- agg_ses
    #group <- agg_ses
    #group <- agg_ror

################################################################################
#                                                                              #                         
# 3. Run agent-based model                                                     #                         
#                                                                              #                         
################################################################################

# Derive total number of women by SES and/or region
    headcount  <- rep(1,length(women1880sim[,1]))
    n_by_group <- aggregate(headcount,list(group),FUN="sum")
    n_by_group.reg <- aggregate(headcount,list(group.reg),FUN="sum")
    n_by_group.ses <- aggregate(headcount,list(group.ses),FUN="sum")
    n_group    <- length(n_by_group[,1])

# Prepare dataframe with results
    res_dat <- data.frame(matrix(nrow=n_group,ncol=n_time+13))
    colnames(res_dat) <- c("scenario","adapt_alg","adapt_from","wr","ws","x",
                           "ses","ses_name","regkey","regname","age","n",
                           paste(0:n_time))

    res_dat[,1] <- scen
    res_dat[,2] <- ad_alg
    res_dat[,3] <- adapt_from
    res_dat[,4] <- wr
    res_dat[,5] <- ws
    res_dat[,6] <- x

    if (group[1]==agg_ses_ror[1]) {
       res_dat[,7] <- sort(rep(unique(ses_code),25))
       res_dat[,8] <- sort(rep(ses_name,25))
       res_dat[,9] <- rep(sort(unique(ror_code)),3)
       res_dat[,10] <- rep(regnames_ord$name,3)
       res_dat[,11] <- rep("20-49",n_group)
       res_dat[,12] <- n_by_group$x
    }   
    if (group[1]==agg_ses[1]) {
       res_dat[,7] <- sort(unique(ses_code))
       res_dat[,8] <- ses_name
       res_dat[,8] <- rep("SWE",n_group)
       res_dat[,10] <- rep("Sweden",n_group)
       res_dat[,11] <- rep("20-49",n_group)
       res_dat[,12] <- n_by_group$x
    } 

    if (group[1]==agg_ror[1]) {
       res_dat[,7] <- rep("ALL",n_group)
       res_dat[,8] <- rep("All SES",n_group)
       res_dat[,9] <- sort(unique(ror_code))
       res_dat[,10] <- regnames_ord$name
       res_dat[,11] <- rep("20-49",n_group)
       res_dat[,12] <- n_by_group$x
    }

# Save Location
#path3 <- paste("/Results/",scen,"_",ad_alg,
#               "_",adapt_from,"_wr_",wr,"_ws_",ws,"_x_",x,add,sep="")
    path3 <- paste("/Model results")
    dir.create(file.path(paste(main,path,path2,path3,sep="")))
    setwd(paste(main,path,path2,path3,sep=""))

    regs <- rep(1:25,3)
    ror <- sprintf("%02.0f",1:25)

# ABM Model
    for (i in 12:n_sim) {
        women <- women1880sim
        women$adapt[sample_sparse[,i]] <- 1
        for (j in 1:n_time) {
            adapters.ror <- aggregate(women$adapt,list(group),FUN="sum")
            sh.ad.r.ses <- adapters.ror$x/n_by_group$x*100
            sp <- split(sh.ad.r.ses,regs)
            ses.max.reg <- unlist(c(lapply(sp,max)))
# We consider our 25 regions and three groups for the outcome
            res_dat[,j+12] <- adapters.ror$x
            ord_r <- match(group.reg,ror)
            ord_r.ses <- match(group,adapters.ror$Group.1)
            ord_b.ses <- match(group.b,adapters.ror$Group.1)
            SR.van <- ses.max.reg[ord_r]
            SR.ses <- sh.ad.r.ses[ord_r.ses]
            SB.ses <- sh.ad.r.ses[ord_b.ses]
            random_num <- runif(n,0,100)
            #RA <- ((SR*w+SB*(1-w))/100)*x
            RA <- (((SR.ses*(1-ws)+SR.van*ws)*wr+SB.ses*(1-wr))/100)*x
            a <- which(RA>random_num)  
            women$adapt[a] <- 1 
        }
        adapters.ror <- aggregate(women$adapt,list(group),FUN="sum")
        res_dat[,n_time+13] <- adapters.ror$x
        write.table(res_dat,paste("ABM_",scen,"_",ad_alg,"_",adapt_from,"_wr_",wr,
                                  "_ws_",ws,"_x_",x,"_sim_",i,".csv",sep=""),
                                  sep=",",row.names=F)
    }
}