rm(list=ls())
library(foreign)
library(dplyr)
library(ggplot2)
library(ipumsr)
library(data.table)
library(dtplyr)

setwd("~/Documents/Classes/MPSpatialDemog/IPUMS/")

cps_data <- fread("./ipumsi_00001.csv")

cps_ddi <- read_ipums_ddi("ipumsi_00001.xml")
cps_ddi
cps_data <- read_ipums_micro(cps_ddi, verbose=TRUE)

fwrite(cps_data, "./ipumsi_00001.csv", nThread=6)
