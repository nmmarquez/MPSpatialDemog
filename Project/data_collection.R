rm(list=ls())
setwd("~/Documents/Classes/datascipop/Lab_Week2_Twitter_Data_Canvas/")
# set wd to this files location

library(ROAuth)
library(streamR)
library(twitteR)
library(stringr)
library(plyr)
library(dplyr)
library(rjson)
library(pander)
source('functions_by_pablobarbera.R')
source('./sentiment.r')

# run only on initial run through

load(file="../credentials/secretkeys.Rdata")
load("../credentials/twitCred.Rdata")

setup_twitter_oauth(api_key, api_secret, access_token, access_token_secret)

tweets_DACA <- searchTwitter("DACA", n=500)
