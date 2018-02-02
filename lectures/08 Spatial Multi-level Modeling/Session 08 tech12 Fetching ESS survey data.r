################################################################################
#                                                                              #                          
# IDEM 156 - Spatial Demography Course 2018 - Technical File 12                #                          
# Extracting European Social Survey data                                       #                          
# Sebastian Klüsener, MPIDR                                                    #                          
#                                                                              #                          
################################################################################


# Erase all objects in memory
rm(list = ls(all = TRUE))

# Load libraries
library(ess)

# This library has been written by Jorge Cimentada. He actually does not use
# an application programming interface for this, but basically fetches the data
# directly from the ESS webpage. Thus, if the ESS webpage would change, this
# library would need to be adapted. 

# Available ESS countries
show_countries()

# Available rounds for a specific country 
# (Belgium participated oin round 2 and 4)
show_country_rounds("Belgium")

# Fetch round 2 ESS data for Turkey. If you are a registered ESS users,
# you would just need to add your email address here which you used to 
# register.
ess_country <-
  ess_country(
    country = "Belgium",
    rounds = c(8),
    your_email = "your_email@random.com"
  )

# Data is loaded as a tibble object (first introduced in the dplyer library),
# which is an alternative way to store a data frame.
head(ess_country)

# To turn it into a data.frame
mydata <- as.data.frame(ess_country)

# Number of observations by region
table(mydata$region)



# Other functions
# Fetch all rounds for one country
austria_all_rounds <- ess_all_cntrounds("Austria", "your_email@random.com")

# Fetch all countries in specific rounds
rounds_1_3_6 <-  ess_rounds(c(1, 3, 6), "your_email@random.com

# Fetch all rounds
ess_all_rounds("your_email@random.com")

# The library also allows to download the data as Stata-files (see ??ess)
