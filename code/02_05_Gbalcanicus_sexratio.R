#### ------------------- Script for calculating sex ratio of Gammarus balcanicus --------------------- ##
## Does leafbreakdown correlate with Gammarid abundance?
## Gamma biomass?



###### Section 1 Loading packages and file paths ######
library(data.table)
library(dplyr)
library(readxl)
library(RPostgreSQL)

## load connections
wd <- setwd("/home/moritz/Nextcloud/PHD/Projects/Romania/UKL_RO_Macroinv")
## check the wd
getwd()

## create path to data folders
datadir <- '/home/moritz/Nextcloud/PHD/Projects/Romania/UKL_RO_Macroinv/data'


###### Section 2 Loading data from server ######
## load gammarus data
gbal <- read.csv("./data/invertebrate_data/g_balcanicus.csv", stringsAsFactors = FALSE)
## remove all entries that are not from June
gbal_june <- gbal[gbal$date >= "2016-05-24",]
## remove all juveniles
gbal_june <- gbal_june[!gbal_june$life_stage %in% c("j", "sum"),]
## remove all rows but site, life_stage and sum
gbal_june <- gbal_june[, c(1,4,32)]
## remove O.DOWN, rename O.UP
gbal_june[gbal_june$site == "O.UP", 1] <- c(rep("O", 3))
gbal_june <- gbal_june[gbal_june$site != "O.DOWN",]
## sum up females
gbal_female <- gbal_june %>%
  group_by(site) %>%
  filter(life_stage != "m") %>%
  summarize(f = sum(sum))

gbal_ratio <- merge(gbal_female, gbal_june[gbal_june$life_stage == "m", c(1,3)], by = "site")
names(gbal_ratio)[3] <- "m"
###### Section 3 Calculate sex ratio ######

gbal_ratio$ratio <- gbal_ratio$f / gbal_ratio$m

## safe data frame
write.csv(gbal_ratio, "./data/parasite_data/sex_ratio.csv", row.names = FALSE)
