################################################################################
#                                                                              #
# IDEM 156 - Spatial Demography Course 2018                                    #
# Theil's Index of Inequality (essentials)                                     #
# Sebastian Kl?sener, Giancarlo Camarda, MPIDR                                 #
#                                                                              #
################################################################################

################################################################################
# 1)   Load and data and define upper hierarchy                                #
################################################################################

# Data should be organized the following way: 
# first column: in this case "ID" (e.g. of individual/region); 
# second column: group variable: in this case "g" (e.g. region/country)
# third-n column: data by time, second-n row: regions)

# Theil index can only be calculated with data having values >0
# It can take ordinal and continuous data

# Load data
ineq_data <-  read.table(file="tfr_sel_european_states_1991-2008.csv", 
                         sep=",", header=TRUE)
 
# This command defines the upper hierarchy for the decomposed Theil measure
group <- ineq_data$g


################################################################################
# 2)   Calculate Theil (if you followed the instructions above correctly, it   #
#      is not necessary to adopt the code in section 2                         #
#                                                                              #
# 2.1) Preparations                                                            #
################################################################################

# Split data by groups
data_groups <- split(ineq_data, group)

# Drop id and group column 
data_idp <- data.frame(ineq_data[,-c(1:2)])

# Derive information on number of t and n
years <- length(data_idp)
groups <- length(data_groups)

# Create empty data matrix for results
Theil.data <- matrix(nrow=5+(2*groups),ncol=years, 
                     data=rep(0,(5+(2*groups))*years))


################################################################################
# 2.2) Calculations                                                            #
################################################################################

for (i in 1:years) {
    yi <- data_idp[,i]

    # Calculating Global Theil
    ni <- length(yi)
    nil <- rep(ni, ni)
    sum_yi <- rep(sum(yi), ni)
    si <- yi/sum_yi
    Theil <- sum(si*log(nil*si))

    # Calculating Decomposed Theil
    TBgj <- as.vector(rep(0, groups))
    TWgj <- as.vector(rep(0, groups))
    TBgjt <- as.vector(rep(0, groups))
    TWgjt <- as.vector(rep(0, groups))

    for (j in 1:groups) {
        y_gj <- data.frame(data_groups[[j]][,-c(1:2)])
        yi_gj <- y_gj[,i]
        n_gj <- length(yi_gj)
        nil_gj <- rep(n_gj, n_gj)
        sum_yi_gj <- rep(sum(yi), n_gj)
        sum_yigj_gj <- rep(sum(yi_gj), n_gj)

        s_gj <- sum(yi_gj/sum_yi_gj)
        si_gj <- yi_gj/sum_yigj_gj
        TBgj[j]  <- s_gj*log(ni/n_gj*s_gj)
        TWgj[j] <-  s_gj*sum(si_gj*log(nil_gj*si_gj))
        TBgjt[j] <- paste("TB.s",j)
        TWgjt[j] <- paste("TW.s",j)
    }
    # "between group" component of inequality
    TB <- sum(TBgj)
    # "within group" component of inequality
    TW <- sum(TWgj)
    T1 <- TB+TW
    # Polarization index
    P <- TB/TW
    # Store results in data frame
    Theil.data[,i] <- c(Theil, TB, TW, T1, P, TBgj, TWgj)
}


################################################################################
# 2.3) Preparation of results table                                            #
################################################################################

# Generate row names with Country ID
cnames <- names(data_groups)
cnamesB <- paste("TB.s",cnames)
cnamesW <- paste("TW.s",cnames)

# Rename row names and column names 
rownames(Theil.data) <- c("Theil.s","TB.s","TW.s","T1.s","P",cnamesB,cnamesW)
colnames(Theil.data) <- c(colnames(ineq_data[,-c(1:2)])) 

# Results
Theil.data.r <- t(round(Theil.data,6))
Theil.data.r


################################################################################
# 3) Export of table                                                           #
################################################################################

# Theil.s: Overall Theil index
# TB.s: Overall Theil Between index
# TW.s: Overall Theil Within index
# T1.s: Overall Theil index derived by adding Theil Between and Theil within
# P: Polarization Index (Theil Between divided by Theil Within)
# Export data                                     
Theil.data.exp <- data.frame(row.names(Theil.data),Theil.data)
colnames(Theil.data.exp)[1] <- c("Name")
write.table(file="Theil.data.Europe.csv", sep=",", Theil.data.exp, row.names=F)


################################################################################
# 4) In case you want to generate a pseudo p-value for the Polarization index  #
#    based on a permutation prodedure, you would run the following code        #
################################################################################

# Load dataset
ineq_datac <- read.table(file="tfr_sel_european_states_1991-2008.csv", 
                         sep=",", header=TRUE)

# In this example we just focus on regional data for Austria and Spain
# But the code could also be run over the complete dataset.
#ineq_data <- ineq_datac
ineq_data <- ineq_datac[ineq_datac$g=="AUT"|ineq_datac$g=="ESP",]

# Define number of premutations
nper <- 99

# Define group level variable
group <- paste(ineq_data$g)


################################################################################
# 5)   Calculate Theil (if you followed the instructions under 1 correctly, it #
#      is not necessary to adopt the code in section 5                         #
# 5.1) Preparations                                                            #
################################################################################

# Define matrix for results
res.mat <- matrix(nrow=ncol(ineq_data)-2,ncol=nper+1)

# Create a list to store the random samples
perdata <- list()                               

# Sample dataset columnwise  
for (k in 1:nper) {

  perdata[[k]] <- data.frame(ineq_data[,c(1,2)],
                             apply(data.frame(ineq_data[,-c(1:2)]), 2, sample)) 
}

# The last list element is the observed dataset
perdata[[nper+1]] <- ineq_data 


################################################################################
# 5.2) Permutation Procedure                                                   #
################################################################################

# Calculate Theil Index for number of permutations
for (k in 1:(nper+1)) {
    # Split data by groups
    data_groups <- split(perdata[[k]], group)
  
    # Drop id and group column 
    data_idp <- data.frame(perdata[[k]][,-c(1:2)])
  
    # Derive information on number of t and n
    years <- length(data_idp)
    groups <- length(data_groups)
  
    for (i in 1:years) {
        yi <- data_idp[,i]
    
        # Calculating Global Theil
        ni <- length(yi)
        nil <- rep(ni, ni)
        sum_yi <- rep(sum(yi), ni)
        si <- yi/sum_yi
        Theil <- sum(si*log(nil*si))
    
        # Calculating Decomposed Theil
        TBgj <- as.vector(rep(0, groups))
        TWgj <- as.vector(rep(0, groups))
        TBgjt <- as.vector(rep(0, groups))
        TWgjt <- as.vector(rep(0, groups))
    
        for (j in 1:groups){
            y_gj <- data.frame(data_groups[[j]][,-c(1:2)])
            yi_gj <- y_gj[,i]
            n_gj <- length(yi_gj)
            nil_gj <- rep(n_gj, n_gj)
            sum_yi_gj <- rep(sum(yi), n_gj)
            sum_yigj_gj <- rep(sum(yi_gj), n_gj)
      
            s_gj <- sum(yi_gj/sum_yi_gj)
            si_gj <- yi_gj/sum_yigj_gj
            TBgj[j]  <- s_gj*log(ni/n_gj*s_gj)
            TWgj[j] <-  s_gj*sum(si_gj*log(nil_gj*si_gj))
            TBgjt[j] <- paste("TB.s",j)
            TWgjt[j] <- paste("TW.s",j)
        }
        # "between group" component of inequality
        TB <- sum(TBgj)
        # "within group" component of inequality
        TW <- sum(TWgj)
        T1 <- TB+TW
        # Polarization index
        P <- TB/TW
        # Store results
        res.mat[i,k] <- c(P)
    }
}

# Derive observed Polarization index
Pobs <- res.mat[,nper+1]

# Sort the rows of the results matrix in ascending order
Psort <- t(apply(res.mat, 1, sort))

# Detect, at which position is the observed P value in the P value derived by
# randomly permuting the values across the regions
pseudo <- c(rep(0,years))
for (i in 1:years) {
    pseudo[i] <- which(Psort[i,] == Pobs[i])
}

# Calculate pseudo-p-value
pseudo.p <- 1-(pseudo/(nper+1))

# Prepare results in dataframe
pseudo.dat <- data.frame(colnames(data.frame(ineq_data[,-c(1:2)])),pseudo.p)
colnames(pseudo.dat) <- c("Year","pseudo.p")


################################################################################
# 6) Export of table with pseudo-values of Polarization index                  #
################################################################################

write.table(file="Polarization_pseudo_p.csv", sep=",", pseudo.dat, row.names=T)
