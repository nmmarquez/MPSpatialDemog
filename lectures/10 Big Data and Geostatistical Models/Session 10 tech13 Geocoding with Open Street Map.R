################################################################################
#                                                                              #                          
# IDEM 156 - Spatial Demography Course 2018 - Technical File 13                #                          
# Geocoding with Open Street Map                                               #                                                                             
# Sebastian Kl?sener, MPIDR                                                    #                          
#                                                                              #                          
################################################################################

# Erase all objects in memory
rm(list = ls(all = TRUE))

# Install packages
#install.packages(c("RJSONIO"))

# Load Library
library(RJSONIO)


################################################################################
#                                                                              #                          
# 1. Example - City Locations                                                  #                          
#                                                                              #
################################################################################

# Here we are defining a dataframe with placenames and ISO-countrycodes in which
# the cities are situated. Be aware that names are not necessarily unique and 
# that there might be several cities with the same name in a country.
# In this example, the location of "New York" will not give you the location of 
# New York City, but of New York, Florida. If you write New York City, you get
# the right location. How place names are stored, you can find out by running
# the search function on http://www.openstreetmap.org/, where a search after 
# "New York" will give you all the places stored under that name.
# You can derive with function also the location of settlements which do not 
# have city status.
cityname <- c("New York", "New York City","Hamburg","Amsterdam")
countrycode <- c("US","US","DE", "NL")
mydata <- data.frame(cityname,countrycode)
colnames(mydata) <- c("CityLong","CountryCode")

# Extracting geocodes (latlong-projected)
# Based on a code suggested by Jochem
# http://stackoverflow.com/questions/13905098/how-to-get-the-longitude-and-latitude-coordinates-from-a-city-name-and-country-i
# Cities
nrow <- nrow(mydata)
counter <- 1
mydata$lon[counter] <- 0
mydata$lat[counter] <- 0
while (counter <= nrow){
  CityName <- gsub(' ','%20',mydata$CityLong[counter]) #remove space for URLs
  CountryCode <- mydata$Country[counter]
  url <- paste(
    "http://nominatim.openstreetmap.org/search?city="
    , CityName
    , "&countrycodes="
    , CountryCode
    , "&limit=9&format=json"
    , sep="")
  x <- fromJSON(url)
  if(is.vector(x)){
    mydata$lon[counter] <- x[[1]]$lon
    mydata$lat[counter] <- x[[1]]$lat    
  }
  counter <- counter + 1
}
mydata

# As a next step you could save the object mydata in a csv-file or turn it in a
# point-shapefle using the code from 
# "Session 03 tech8 Merging csv-file data with shapefile"


################################################################################
#                                                                              #                          
# 2. Example - Address Locations                                               #                          
#                                                                              #
################################################################################

# This time we also define a column with address-information
# Now I get the right location in New York City, although I define the city just
# as "New York".
address <- c("1 Times Square","Jungfernstieg 10","Rembrandtplein 1")
cityname <- c("New York", "Hamburg", "Amsterdam")
countrycode <- c("US", "DE", "NL")
mydata1 <- data.frame(address,cityname,countrycode)
colnames(mydata1) <- c("AddressLong","CityLong","CountryCode")

# Extracting geocodes (latlong-projected)
nrow <- nrow(mydata1)
counter <- 1
mydata1$lon[counter] <- 0
mydata1$lat[counter] <- 0
while (counter <= nrow){
  Address <- gsub(' ','%20',mydata1$AddressLong[counter]) #remove space for URLs
  CityName <- gsub(' ','%20',mydata1$CityLong[counter]) #remove space for URLs
  CountryCode <- mydata1$Country[counter]
  url <- paste(
    "http://nominatim.openstreetmap.org/search?street="
      , Address
      , "&city="
      , CityName
      , "&countrycodes="
      , CountryCode
      , "&limit=9&format=json"
      , sep="")
  x <- fromJSON(url)
  if(is.vector(x)){
    mydata1$lon[counter] <- x[[1]]$lon
    mydata1$lat[counter] <- x[[1]]$lat    
  }
  counter <- counter + 1
}

mydata1

# As a next step you could save the object mydata in a csv-file or turn it in a
# point-shapefle using the code from 
# "Session 03 tech8 Merging csv-file data with shapefile"

