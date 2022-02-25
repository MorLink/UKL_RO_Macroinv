## -------------------------------- Gammarid and Shredder preparation ------------------------------ ##
## This script prepares the dataframe for the shredder community
## 
###### Section 1 Loading packages and file paths ######
library(data.table)
library(dplyr)
library(readxl)
library(RPostgreSQL)

## load connections
os_sys <- switch(Sys.info()[['sysname']],
                 Windows= {print("Windows")},
                 Linux  = {print("Linux")})



wd <- if (os_sys == "Windows") {
  setwd("C:/Users/Moritz/Documents/PHD/Projects/Romania/UKL-Romania-SEM")
} else {
  setwd("/home/moritz/Nextcloud/PHD/Projects/Romania/UKL-Romania-SEM")
}

## check the wd
getwd()

## create path to data folders
datadir <- if (os_sys == "Windows") {
  'C:/Users/Moritz/Nextcloud/PHD/Projects/Romania/Daten'
} else {
  '/home/moritz/Nextcloud/PHD/Projects/Romania/Daten'
}

datadir2 <- if (os_sys == "Windows") {
  'C:/Users/Moritz/Documents/PHD/Projects/Romania/Daten'
} else {
  '/home/moritz/Documents/PHD_Spatial_Files/Projects/Romania/Daten'
}

## enable connection
source(file.path(datadir, "Scripts/00_Connection_Data.R"))
drv = dbDriver("PostgreSQL")
con = dbConnect(drv, dbname = NEW_DB, user = DBuser_m, host = DBhost, port = DBport, password = DBpassword_m)


###### Section 2 Load Gammarid data  and environmental data ######
## size classes and abundances
gam <- dbGetQuery(con, "SELECT * FROM inv_aquatic.gammarus_balcanicus")
## save an offline copy
# write.csv(gam, "/home/moritz/Nextcloud/PHD/Projects/Romania/Daten/Offline_Backups/gam.csv", row.names = FALSE)
## subset gam to samples of second sampling date
gam <- gam[grepl("2016-06", gam$date),]

## load table with size-to-weight ratio
sw_ratio <- dbGetQuery(con, "SELECT * FROM inv_aquatic.g_balcanicus_size_to_weight")

## overwrite the size classes 1.5 to 3.0 in order to get them all in the same group
sw_ratio[1:18,1] <- c(rep(3.5, 18))
## calculate mean weights for the size classes in sw_ratio
group_means <- data.frame(sw_ratio %>%
                            group_by(length_mm) %>%
                            summarise_at(vars(mass_mg), list(mean_weight_g = ~ mean(., na.rm=TRUE), 
                                                             sd_weight_g = ~ sd(., na.rm=TRUE))), 
                          stringsAsFactors = FALSE)

## taxa data with actual number of gammarids
all_jn <- read.csv(file.path(datadir, "MZB/alltaxajune_wide.csv"), stringsAsFactors = FALSE)
all_jn$site <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "K", 
                 "L", "M", "N", "O", "P", "Q", "R", "S", "T")
## keep only Gammarids
all_jn <- all_jn[,c("site", "Gammarus.balcanicus")]

## -------- Get Biomass of all size groups ------##
## transpose gam_sum so sites are columns and size classes are rows
gam_sum <- gam[gam$life_stage == "sum", -c(2,3,4,32)]
nam <- gam_sum$site
rows <- names(gam_sum)
gam_sum <- data.frame(t(gam_sum))
gam_sum <- gam_sum[-1,]
rownames(gam_sum) <- NULL
names(gam_sum) <- nam
gam_sum$size_class <- rows[-1]
gam_sum <- gam_sum[,c(21,1:20)]
## turn all site columns to numeric
gam_sum[nam] <- sapply(gam_sum[nam],as.numeric)
## delete sizeclass 26 since empty and also missing in group_means
gam_sum <- gam_sum[-26,]
## get columns sums
sum_all <- data.frame(colSums(gam_sum[,-1]))
sum_all$site <- rownames(sum_all)
names(sum_all)[1] <- "subset_sum"
## drop O.DOWN and rename O.UP
sum_all <- sum_all[sum_all$site != "O.DOWN",]
sum_all[14,]$site <- "O"

## multiply weights with sum of size classes
biomass <- round(gam_sum[,2:21] * group_means$mean_weight_g, 3)
## get column sums
biomass_all <- data.frame(colSums(biomass))
biomass_all$site <- rownames(biomass_all)
names(biomass_all)[1] <- "subset_biomass"
## drop O.UP and rename O.DOWN
biomass_all <- biomass_all[-15,]
biomass_all[biomass_all$site == "O.UP", 2] <- "O"

## merge sum_all and biomass_all
biomass_fin <- merge(sum_all, biomass_all, by = "site")
names(biomass_fin)
## add column with average body weight
biomass_fin$gam_bodyweight <- biomass_fin$subset_biomass / biomass_fin$subset_sum


## merge biomass_fin with all_jn
biomass_fin <- merge(biomass_fin, all_jn, by = "site")
## calculate biomass of Gammarus balcanicus at each site
biomass_fin$total_biomass <- biomass_fin$gam_bodyweight * biomass_fin$Gammarus.balcanicus
biomass_fin[is.na(biomass_fin$gam_bodyweight),]$gam_bodyweight <- c(0,0)
biomass_fin[is.na(biomass_fin$total_biomass),]$total_biomass <- c(0,0)
names(biomass_fin)[5] <- "total_abund"
# write.csv(biomass_fin, file.path(datadir, "MZB/gammarid_biomass.csv"), row.names = FALSE)


###### Section 3 Generate the Shredder data ######
## ////////////////////// Load MZB DATA //////////////////// ##
## load invertebrate data
inv <- read.csv(file.path(datadir, "/MZB/MZB_Rom_Leon_merged.csv"), stringsAsFactors = FALSE, na.strings=c("","NA"))

setDT(inv)

## keep only data from June and April
mzb_mi <- inv[date >= "2016-06-01" | date <= "2016-04-30",,]
mzb_mi$month <- ifelse(mzb_mi$date <= "2016-04-30",
                       paste0("April"),
                       ifelse(mzb_mi$date >= "2016-06-01",
                              paste0("June"),
                              NA))

mzb_mi[mzb_mi$site == "O.UP", ]$site <- c(rep("O",19))
mzb_mi <- mzb_mi[site != "O.DOWN",]

#####################################################
## ----- continue with traits on taxa ----- ##
#####################################################
## merge trait data from server with mzb_mi
traits_spec <- dbGetQuery(con, "SELECT * FROM inv_aquatic.traits_species")
traits_gen <- dbGetQuery(con, "SELECT * FROM inv_aquatic.traits_genus")
traits_fam <- dbGetQuery(con, "SELECT * FROM inv_aquatic.traits_family")
traits_oth <- dbGetQuery(con, "SELECT * FROM inv_aquatic.traits_class")
mzb_spec <- mzb_mi[!is.na(species),]
mzb_spec <- merge.data.table(mzb_spec, traits_spec, by.x = "species", by.y = "taxon", all.x = TRUE)
mzb_gen <- mzb_mi[is.na(species) & !is.na(genus),]
mzb_gen <- merge.data.table(mzb_gen, traits_gen, by = "genus", all.x = TRUE)
mzb_fam <- mzb_mi[is.na(genus) & !is.na(family),]
mzb_fam <- merge.data.table(mzb_fam, traits_fam, by = "family", all.x = TRUE)
mzb_class <- mzb_mi[is.na(order) & is.na(family) & class %in% c("Clitellata", "Gastropoda", "Bivalvia"),]
mzb_class <- merge.data.table(mzb_class, traits_oth, by = "class", all.x = TRUE)

mzb_rest <- mzb_mi[(is.na(family) & phylum %like% "Arthropoda") | (is.na(family) & phylum %like% "Mollusca" & is.na(class)),]

## add empty trait columns for rbind later
trait_cols <- data.frame(feed_gra = c(rep(NA, nrow(mzb_rest))),
                         feed_min = c(rep(NA, nrow(mzb_rest))),
                         feed_xyl = c(rep(NA, nrow(mzb_rest))),
                         feed_shr = c(rep(NA, nrow(mzb_rest))), 
                         feed_gat = c(rep(NA, nrow(mzb_rest))),
                         feed_aff = c(rep(NA, nrow(mzb_rest))), 
                         feed_pff = c(rep(NA, nrow(mzb_rest))), 
                         feed_pre = c(rep(NA, nrow(mzb_rest))), 
                         feed_par = c(rep(NA, nrow(mzb_rest))), 
                         feed_oth = c(rep(NA, nrow(mzb_rest))), 
                         loco_sws = c(rep(NA, nrow(mzb_rest))), 
                         loco_swd = c(rep(NA, nrow(mzb_rest))),
                         loco_bub = c(rep(NA, nrow(mzb_rest))), 
                         loco_spw = c(rep(NA, nrow(mzb_rest))), 
                         loco_ses = c(rep(NA, nrow(mzb_rest))), 
                         loco_oth = c(rep(NA, nrow(mzb_rest))), 
                         resp_teg = c(rep(NA, nrow(mzb_rest))), 
                         resp_gil = c(rep(NA, nrow(mzb_rest))), 
                         resp_pls = c(rep(NA, nrow(mzb_rest))),
                         resp_spi = c(rep(NA, nrow(mzb_rest))), 
                         resp_ves = c(rep(NA, nrow(mzb_rest))), 
                         resp_tap = c(rep(NA, nrow(mzb_rest))), 
                         resp_sur = c(rep(NA, nrow(mzb_rest))),
                         stringsAsFactors = FALSE)
mzb_rest <- cbind(mzb_rest, trait_cols)

## bring all files back together
mzb_trait <- rbind(mzb_spec, mzb_gen,mzb_fam,mzb_class,mzb_rest)

#####################################################
##  -- Shredder taxa diversity --- ##
#####################################################
## create column with lowest identification level
shr_trait_1 <- mzb_trait[!is.na(feed_shr) & feed_shr > 0,]
shr_trait_2 <- mzb_trait[!is.na(feed_shr) & feed_shr > 3,]
shr_trait_3 <- mzb_trait[!is.na(feed_shr) & feed_shr > 5,]

all_shr <- shr_trait_1[month == "June", .(all_shr = sum(abundance)), by = .(site)]
medium_shr <- shr_trait_2[month == "June", .(med_shr = sum(abundance)), by = .(site)]
strong_shr <- shr_trait_3[month == "June", .(str_shr = sum(abundance)), by = .(site)]

## merge all three shredder tables
shr_abun <- merge(all_shr, medium_shr, by = "site", all = TRUE)
shr_abun <- merge(shr_abun, strong_shr, by = "site", all = TRUE)
shr_abun[is.na(med_shr), ]$med_shr <- c(0,0)
shr_abun[is.na(str_shr), ]$str_shr <- c(0,0)
# write.csv(shr_abun, file.path(datadir, "MZB/shredder_abundance.csv"), row.names = FALSE)


###### Section 4 Merge Gamma data and Shredder data ######
shr_com <- merge(biomass_fin[,c(1,4,5,6)], shr_abun[,c(1,2)], by = "site", all = TRUE)
shr_com$nongam_shr <- shr_com$all_shr - shr_com$total_abund
shr_com <- shr_com[,-4]
names(shr_com)[3] <- "gam_shr"

## write as csv file
write.csv(shr_com, file.path(datadir, "/MZB/shredder_community.csv"), row.names = FALSE)
