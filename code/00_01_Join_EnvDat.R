##// Script for selection of Environmental Variables //##


########### Section 1 ###########
## -------------------------------------- enable connections ------------------------------ ##
## packages
library(dplyr)
require(RPostgreSQL)
library(data.table)
library(ggplot2)



setwd("/home/moritz/Nextcloud/PHD/Projects/Romania/UKL_RO_Macroinv")
## check the wd
getwd()

## create path to data folders
datadir <- '/home/moritz/Nextcloud/PHD/Projects/Romania/UKL_RO_Macroinv/data'

## enable connection
source("/home/moritz/Nextcloud/PHD/Projects/General_Scripts/Connection_Data.R")
drv = dbDriver("PostgreSQL")
con = dbConnect(drv, dbname = NEW_DB, user = DBuser_m, host = DBhost, port = DBport, password = DBpassword_m)


########### Section 2 ###########
## -------------------------------------- Data Preparation ------------------------------ ##

## ////////////////////// Load PHCH DATA //////////////////// ##
phch <- dbGetQuery(con, "SELECT * FROM phch.reg_phch")
## use first values around first mzb sampling, leaf decomp 
## and create mean value for variables
phch_sub <- phch[phch$variable %in% c("water_temp", "ph", "el_cond", "o2_mg", "cl", "no3", "no2", "no4",
                                      "nh4", "po4", "flow_1", "flow_2", "flow_3"), ]

phch_sub <- phch_sub[phch_sub$date >= "2016-06-01" | phch_sub$date <= "2016-04-30",]
## create a month variable
phch_sub$month <- ifelse(phch_sub$date <= "2016-04-30",
                         paste0("April"),
                         ifelse(phch_sub$date >= "2016-06-01",
                                paste0("June"),
                                NA))
## cast into wide format
phch_wide <- reshape2::dcast(phch_sub, site + month ~ variable,  value.var = "value")
## average flow measurements
phch_wide$flow <- rowMeans(phch_wide[,5:7])
phch_wide <- phch_wide[,-c(5:7)]



## ////////////////////// Load PESTICIDE DATA //////////////////// ##
pesticides <- dbGetQuery(con, "SELECT * FROM pesticide_data.pest_gradient")
## only keep mean, meadian and max TU values
clmns <- c( "site", "max_sumtu_iv")
tu_ms <- pesticides[clmns]




## ////////////////////// Load HABITAT DATA //////////////////// ##
habitat <- dbGetQuery(con, "SELECT * FROM habitat.sites_habitat")
## use second sampling here
habitat_2 <- habitat[habitat$date >= "2016-06-01"| habitat$date <= "2016-04-30", ]
habitat_2$month <- ifelse(habitat_2$date <= "2016-04-30",
                          paste0("April"),
                          ifelse(habitat_2$date >= "2016-06-01",
                                 paste0("June"),
                                 NA))
## variable for bank cover
habitat_bc <- habitat_2[habitat_2$variable %in% c("bank_cover_meadow", "bank_w/o_veg"), ]
habitat_bc$value <- as.numeric(habitat_bc$value)
## sum up values per site and take mean
setDT(habitat_bc)
rip_veg <- habitat_bc[, sum(value)/2, keyby = .(site, month)]
names(rip_veg)[3] <- "low_rip_veg"

## variable for agriculture on banks and nearby fields
habitat_agri <- habitat_2[habitat_2$variable %in% c("bank_cover_agri", "land_use_agri"),]
habitat_agri$value <- as.numeric(habitat_agri$value)
## sum up values per site and take mean
setDT(habitat_agri)
rip_agri <- habitat_agri[, sum(value)/2, keyby = .(variable, site, month)]
## change into wide format
rip_agri <- reshape2::dcast(rip_agri, site + month ~ variable, value.var = "V1")

## get the field distances left and right and keep only the minimum
## attach results to rip_veg
field_dist <- habitat[habitat$date >= "2016-06-01"| habitat$date <= "2016-04-30", ]
field_dist$month <- ifelse(field_dist$date <= "2016-04-30",
                           paste0("April"),
                           ifelse(field_dist$date >= "2016-06-01",
                                  paste0("June"),
                                  NA))

field_dist <- field_dist[field_dist$variable %in% c("dist_field"),]
field_dist$value <- as.numeric(field_dist$value)
setDT(field_dist)
min_field <- field_dist[, min(value), keyby = .(site, month)]
names(min_field)[3] <- "min_dist_field"
## merge rip_veg with min_field
rip_hab <- merge(rip_veg, min_field, by = c("site", "month"), all = TRUE)
rip_hab <- merge(rip_hab, rip_agri[,c(1,2,4)], by = c("site", "month"), all = TRUE)


## extract variable "shading", "riffels"
shading <- habitat[habitat$date >= "2016-06-01"| habitat$date <= "2016-04-30", ]
shading$month <- ifelse(shading$date <= "2016-04-30",
                        paste0("April"),
                        ifelse(shading$date >= "2016-06-01",
                               paste0("June"),
                               NA))
shading <- shading[shading$variable %in% c("shading_from_south", "riffels", "curvature"), c(1,10,4,7)]
shading <- reshape2::dcast(shading, site + month ~ variable, value.var = "value")
names(shading)[c(4,5)] <- c("riffles", "shading")
## merge with rip_hab
rip_hab <- data.frame(merge(rip_hab, shading, by  = c("site", "month"), all = TRUE))
num_col <- c("riffles", "shading", "curvature")
rip_hab[num_col] <- sapply(rip_hab[num_col], as.numeric)




## ////////////////////// Load SUBSTRATE DATA //////////////////// ##
substrate <- dbGetQuery(con, "SELECT * FROM habitat.sites_substrate")
## second sampling here
substrate_2 <- substrate[substrate$date >= "2016-06-01" | substrate$date <= "2016-04-30",]
substrate_2$month <- ifelse(substrate_2$date <= "2016-04-30",
                            paste0("April"),
                            ifelse(substrate_2$date >= "2016-06-01",
                                   paste0("June"),
                                   NA))
substrate_2 <- substrate_2[substrate_2$variable %in% c("Psammal", "Argyllal"),]
substrate_2$value <- as.numeric(substrate_2$value)
## sum up values for both variables
setDT(substrate_2)
stream_sub <- substrate_2[, sum(value), by = .(site,month)]
names(stream_sub)[3] <- "fine_sed"
## divide by 100 to get a range from 0 to 1
stream_sub$fine_sed <- stream_sub$fine_sed/100


## ////////////////////// Load LAND USE DATA //////////////////// ##
ext_agri <- dbGetQuery(con, "SELECT * FROM spatial_derived.land_use")
agri <- ext_agri[,c(1,22:27)]




## ////////////////////// Load RIPARIAN LAND USE DATA //////////////////// ##
rip_lu <- dbGetQuery(con, "SELECT * FROM spatial_derived.riparian_land_use_all")
rip_lu <- rip_lu[,c(1,20:23)]
names(rip_lu)[2:5] <- c("rip_urban", "rip_forest", "rip_past", "rip_agri")


## ////////////////////// Load UPSTREAM REFUGE DATA //////////////////// ##
refuge_dist <- read.csv("/home/moritz/Nextcloud/PHD/Projects/Romania/Daten/refuge_dist.csv", stringsAsFactors = FALSE)
refuge_dist$refuge_forest <- as.integer(ifelse(refuge_dist$type == "forest",
                                               1,
                                               0))



########### Section 3 ###########
## -------------------------------------- merge all variables into one table ------------------------------ ##
data <- merge(phch_wide, tu_ms, by = "site", all = TRUE)
data <- merge(data, rip_hab, by = c("site", "month"), all=TRUE)
data <- merge(data, stream_sub, by = c("site", "month"), all=TRUE)
data <- merge(data, agri, by = "site", all = TRUE)
data <- merge(data, rip_lu, by = "site", all = TRUE)
data <- merge(data, refuge_dist[,c(1,4)], by = "site", all = TRUE)

## combine O and Oup, remove Odown
data[grep("^O$", data$site),]
data[grep("O", data$site),]
# data[data$site == "O",2] <- 15.32515
data[data$site == "O" & data$month == "June", "min_dist_field"] <- 4
data[data$site == "O" & data$month == "June", "curvature"] <- 7
data[data$site == "O" & data$month == "June", "riffles"] <- 60
data[data$site == "O" & data$month == "June", "fine_sed"] <- 0.1
## replace NA with 0 in refuge forest
data[is.na(data$refuge_forest),"refuge_forest"] <- c(rep(0, 13))

## drop all other columns with Oup and Odown column
data <- data[-c(29:30),]

## change value range of columns in percentage to 0 - 1
## low_rip_veg, land_use_agri, riffles, shading, urban, ext_agri, int_agri, nat_nonforest, forest, sum_agri
percent_data <- c("low_rip_veg", "land_use_agri", "riffles", "shading", "urban", "pastures", 
                  "agriculture", "nat_nonforest", "forest")
data[percent_data] <- data[percent_data] / 100
## divide total area by 10^6
data$totalarea <- data$totalarea / 10^6


## round values to 3 digits
data <- data %>%                  
  mutate_if(is.numeric,
            round,
            digits = 3)

## change column name refuge_forest to refugium
names(data)[31] <- "refugium"

########### Section 3 ###########

## write to file
write.csv(data, file.path(datadir, "environmental_data/env_predictors.csv"), row.names = FALSE)
