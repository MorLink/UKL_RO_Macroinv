## ------------------------------------ Script for creating Abundance data frame ----------------- ##



###### Section 1 Load packages and set paths ######
library(data.table)
library(dplyr)
library(vegan)
library(readxl)

## connection to server
os_sys <- switch(Sys.info()[['sysname']],
                 Windows= {print("Windows")},
                 Linux  = {print("Linux")})



wd <- setwd("/home/moritz/Nextcloud/PHD/Projects/Romania/UKL_RO_Macroinv")
## check the wd
getwd()

## create path to data folders
datadir <- '/home/moritz/Nextcloud/PHD/Projects/Romania/UKL_RO_Macroinv/data'


###### Section 2 Load and clean taxa data ######
## ------------------------ Taxa Data ---------------------------- ##
## load the taxa data
mzb_dat <- read.csv(file.path(datadir, "invertebrate_data/macroinv_full_ap_jn.csv"), stringsAsFactors = FALSE, na.strings=c("","NA"))

## rename O.UP and remove O.DOWN
mzb_dat[mzb_dat$site %like% "O.UP",]$site <- c(rep("O", 19))
mzb_dat <- mzb_dat[!mzb_dat$site %like% "O.DOWN",]
setDT(mzb_dat)

## drop sites X, Y, Z
mzb_dat <- mzb_dat[!site %in% c("X", "Y", "Z"),]

## keep only data from June and April
mzb_dat <- mzb_dat[date >= "2016-06-01" | date <= "2016-04-30",,]
mzb_dat$month <- ifelse(mzb_dat$date <= "2016-04-30",
                        paste0("April"),
                        ifelse(mzb_dat$date >= "2016-06-01",
                               paste0("June"),
                               NA))

## create column with lowest identification level
mzb_dat$taxa <- ifelse(!is.na(mzb_dat$species),
                       paste0(mzb_dat$species),
                       ifelse(!is.na(mzb_dat$genus),
                              paste0(mzb_dat$genus, ".sp."),
                              ifelse(!is.na(mzb_dat$family),
                                     paste0(mzb_dat$family, ".Gen.sp."),
                                     ifelse(!is.na(mzb_dat$order),
                                            paste0(mzb_dat$order),
                                            ifelse(!is.na(mzb_dat$sub_class),
                                                   paste0(mzb_dat$sub_class),
                                                   ifelse(!is.na(mzb_dat$class),
                                                          paste0(mzb_dat$class),
                                                          ifelse(!is.na(mzb_dat$phylum),
                                                                 paste0(mzb_dat$phylum),
                                                                 NA)))))))

## remove Copepoda, Ostracoda and Trombidiformes
mzb_dat <- mzb_dat[!taxa %in% c("Copepoda", "Ostracoda", "Trombidiformes"),]
mzb_june <- mzb_dat[month == "June",]



###### Section 3 Create species tables in wide format ######
## Aggregate species on family level
## Start with June
fam_june <- mzb_june[!is.na(family),]
ad_dat <- mzb_june[is.na(family), c("taxa", "abundance", "site")]

## overwrite entry in taxa
fam_june$taxa <- paste0(fam_june$family, ".Gen.sp.")

## combine the duplicates that are now in the fam_june dt
## check for duplicates per site taxa, sum abundances for these
sites_do <- unique(fam_june$site)
# sites_do <- sites_do[1:2]

fam_jfin <- NULL
for(i in seq_along(sites_do))
{
  s <- sites_do[i]
  site_sub <- fam_june[site == s,,]
  
  ## use dplyr to sum abundances of duplicated Taxa entries
  site_clean <- data.frame(site_sub %>%
                             group_by(taxa) %>%
                             summarize(abundance = sum(abundance)),
                           stringsAsFactors = FALSE)
  ## merge site_clean and site_sub
  site_clean <- merge(site_clean, site_sub[,c("site", "taxa")],
                      by = "taxa")
  site_clean <- site_clean[!duplicated(site_clean$taxa),]
  fam_jfin <- rbind(fam_jfin, site_clean)
}

## rbind fam_jfin and ad_dat to include Annelida and Mollusca
fam_jfin <- rbind(fam_jfin, ad_dat)

## cast into wide format, column names are the species names,
fam_wide_jn <- reshape2::dcast(fam_jfin, site ~ taxa, value.var = "abundance")
## turn all NA into 0
fam_wide_jn[is.na(fam_wide_jn)] <- 0
## drop site-column
fam_wide_jn <- fam_wide_jn[,-1]

## ----- check how many taxa have very few occurences ---- ##
## turn into presence absence data
pa_wide <- decostand(fam_wide_jn, "pa")
## get all taxa with presence on less than 1 sites
pa_few <- names(pa_wide)[colSums(pa_wide) == 1]

## drop all taxa which occur at only one site
fjune_wide <- fam_wide_jn[, -which(names(fam_wide_jn) %in% pa_few)]

## write june_wide to csv_file
write.csv(fjune_wide, file.path(datadir, "invertebrate_data/familyjune_wide.csv"), row.names = FALSE)
