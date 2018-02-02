################################################################################
#                                                                              #
# IDEM 156 - Spatial Demography Course 2018                                    #
# Theil's Index of Inequality                                                  #
# Sebastian Kl?sener, Giancarlo Camarda, MPIDR                                 #
#                                                                              #
################################################################################


################################################################################
# 1)   Load and data and define upper hierarchy                                #
################################################################################

# Data should be organized the following way: 
# first column: ID (e.g. of individual/region); 
# second column: group variable (e.g. region/country); 
# third-n column: data by time, second-n row: regions)
# Theil index can only be calculated with data having values >0
# It can take ordinal and continuous data
ineq_data <-  read.table(file="tfr_sel_european_states_1991-2008.csv", 
                         sep=",", header=TRUE)
 
# This commands the upper hierarchy g for the decomposed Theil measure
group <- ineq_data$g


################################################################################
# 2)   Calculate Theil                                                         #
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
    # Store results in data frame
    Theil.data[,i] <- c(Theil, TB, TW, T1, P, TBgj, TWgj)
}


################################################################################
# 3) Preparation of results table                                              #
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

# Export data                                     
Theil.data.exp <- data.frame(row.names(Theil.data),Theil.data)
colnames(Theil.data.exp)[1] <- c("Name")
write.table(file="Theil.data.Europe.csv", sep=",", Theil.data.exp, row.names=F)


################################################################################
# 4) Plot Theil - Specifically adapted to the dataset                          #
################################################################################

# Plot information
colnames(Theil.data) <- c(1991:2008) 
Theil.datat <- data.frame(t(Theil.data))

# Plot Theil, Theil Between, Theil Within and Polarization
par(mfrow=c(2,2))
# Plot Theil
plot(Theil.datat$Theil.s, type="l", lwd=2, col="black", xaxt="n", 
ylab="Theil Index", ylim=c(0,0.03),xlab="Year", main="Theil Index")
axis(1, at=c(1,6,11,16), labels=c(1991,1996,2001,2006))
legend("topright", c("Theil"), lwd=2, col=c("black"))

# Plot Theil and Within
plot(Theil.datat$Theil.s, type="l", lwd=2, col="black", xaxt="n", 
ylab="Theil Index", ylim=c(0,0.03),xlab="Year", main="Theil and Theil Within")
axis(1, at=c(1,6,11,16), labels=c(1991,1996,2001,2006))
lines(Theil.datat$TW.s, col="orange", lwd=2)
legend("topright", c("Theil","Theil Within"), lwd=2, col=c("black","orange"))

# Plot Theil and Theil Between
plot(Theil.datat$Theil.s, type="l", lwd=2, col="black", xaxt="n", 
ylab="Theil Index", ylim=c(0,0.03),xlab="Year", main="Theil and Theil Between")
axis(1, at=c(1,6,11,16), labels=c(1991,1996,2001,2006))
lines(Theil.datat$TB.s, col="blue3", lwd=2)
legend("topright", c("Theil","Theil Between"), lwd=2, col=c("black","blue3"))

# Plot Polarization Index
plot(Theil.datat$P, type="l", lwd=2, col="darkred", xaxt="n", ylim=c(0,6), 
ylab="Degree of Polarization", xlab="Year", main="Polarization Index")
axis(1, at=c(1,6,11,16), labels=c(1991,1996,2001,2006))
legend("bottomright", c("Polarization Index"), lwd=2, col=c("darkred"))


# Plot Theil within by country 
lenTh <- length(Theil.datat[1,])
lw <- length(cnamesW)
TW <- Theil.datat[,c((lenTh-lw+1):lenTh)]

# Prepare plots
color <- c("red", "darkred", "darkblue", "darkgrey", "green3",
 "orange", "yellow2")
countrynames <- c("Austria","Spain","France","Germany",
 "Iceland","Netherlands","Sweden")


# Plot Theil within by country (with for-loop)
par(mfrow=c(1,2))
plot(Theil.datat$TW.s.AUT, type="l", lwd=2, col="red", xaxt="n", 
 ylim=c(0,0.008), ylab="Theil Index", xlab="Year", main="Theil Within")
axis(1, at=c(1,6,11,16), labels=c(1991, 1996,2001,2006))
for (i in 2:lw) {
    lines(TW[,i], col=color[i], lwd=2)
}
legend("topright", countrynames, lwd=2, col=color)

# Plot Theil between by country (without for-loop)
plot(Theil.datat$TB.s.AUT, type="l", lwd=2, col="red", xaxt="n", 
 ylim=c(-0.06,0.18), ylab="Theil Index", xlab="Year", main="Theil Between")
axis(1, at=c(1,6,11,16), labels=c(1991, 1996,2001,2006))
lines(Theil.datat$TB.s.ESP, col="darkred", lwd=2)
lines(Theil.datat$TB.s.FRA, col="darkblue", lwd=2)
lines(Theil.datat$TB.s.GER, col="darkgrey", lwd=2)
lines(Theil.datat$TB.s.ICE, col="green3", lwd=2)                                                 
lines(Theil.datat$TB.s.NLD, col="orange", lwd=2)                                                    
lines(Theil.datat$TB.s.SWE, col="yellow2", lwd=2)  
legend("topright", countrynames, lwd=2, col=color)


################################################################################
# 5) Randomization Procedure                                                   #
################################################################################

# Create-random permutation
# Define number of premutations
nper <- 99

# Load dataset
ineq_datac <- read.table(file="tfr_sel_european_states_1991-2008.csv", 
                         sep=",", header=TRUE)
ineq_data <- ineq_datac
ineq_data <- ineq_datac[ineq_datac$g=="AUT"|ineq_datac$g=="ESP",]

# Define group level variable
group <- paste(ineq_data$g)

# Define matrix for results
res.mat <- matrix(nrow=ncol(ineq_data)-2,ncol=nper+1)

# Create a list to store the random samples
perdata <- list()                               

# Sample dataset columnwise  
for (k in 1:nper) {

  perdata[[k]] <- data.frame(ineq_data[,c(1,2)],apply(data.frame(ineq_data[,-c(1:2)]), 2, sample)) 
}

# The last list element is the observed dataset
perdata[[nper+1]] <- ineq_data 

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

# Plot development of pseudo.p of Polarization Index over time
# If pseudo.p is close to zero, this indicates that between-group variation
# dominates over within-group variation
par(mfrow=c(1,2))
plot(Pobs, type="l", lwd=2, col="darkred", xaxt="n", ylim=c(0,6), 
 ylab="Polarization Index", xlab="Years", main="Polarization Index")
axis(1, at=c(1,6,11,16), labels=c(1991,1996,2001,2006))
legend("topright", c("Polarization Index"), lwd=2, col=c("darkred"))
plot(pseudo.dat$pseudo.p,type="l", lwd=2,ylim=c(0,1), xaxt="n",
 xlab="Years", ylab="pseudo.p of Polarization Index", main="Pseudo-p")
axis(1, at=c(1,6,11,16), labels=c(1991,1996,2001,2006))
legend("topleft", c("Pseudo p-value"), lwd=2, col=c("black"))
