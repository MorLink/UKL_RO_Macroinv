## ------------------------------------------- Parasite Data ------------------------------------- ##
## How do Parasites react to UFRs?
## Other environmental variables?
## Body size of Gammarids and Parasites?
## Fertility of Gammarids and Parasites?
## Average bodyweight and Parasites?


###### Section 1 Loading packages and file paths ######
library(data.table)
library(dplyr)
library(readxl)

## load connections
wd <- setwd("/home/moritz/Nextcloud/PHD/Projects/Romania/UKL_RO_Macroinv")
## check the wd
getwd()

## create path to data folders
datadir <- '/home/moritz/Nextcloud/PHD/Projects/Romania/UKL_RO_Macroinv/data'

###### Section 2 Load the parasite data ######
## parasite data for Baetis
bae_par <- read.csv(file.path(datadir, "parasite_data/baetis_parasites.csv"), stringsAsFactors = FALSE)
## rename O1 und O2
bae_par[bae_par$site == "O1", "site"] <- "O"
bae_par <- bae_par[-13,]

## Parasite data for Gammarids
gam_par <- read.csv(file.path(datadir, "parasite_data/gammarid_parasites.csv"), stringsAsFactors = FALSE)

## load bodyweight of gammarids
shr_com <- read.csv(file.path(datadir, "invertebrate_data/shredder_community.csv"), stringsAsFactors = FALSE)

## load environmental data
env_jn <- read.csv(file.path(datadir, "environmental_data/env_predictors_june.csv"), stringsAsFactors = FALSE)


###### Section 3 merge data frames ######
gam_dat <- merge(shr_com[,-c(4:8)], gam_par, by = "site")
## calculate average prevalence in all gammarids
gam_dat$av_prev <- round((gam_dat$microsporidia + gam_dat$nematoda + gam_dat$trematoda)/3, 1)

## merge environmental data with gam_dat
gpara_dat <- merge(gam_dat[,c(1,2,3,8)], env_jn, by = "site")

## merge environmental data with bea_par
bae_par$av_prev <- round((bae_par$microsporidia + bae_par$nematoda + bae_par$trematoda)/3, 1)
bpara_dat <- merge(bae_par[,-c(2:5)], env_jn, by = "site")

###### Section 4 model gammarid data with environmental parameters ######
## only 12 observations, too few observations for glmnet and cross validation
## check with basic correlation

## Test variables pastures, agriculture and max_sumTUiv
gp_past <- cor.test(gpara_dat$av_prev, gpara_dat$pastures)
gp_past$p.value * 3 ## Bonferoni Correction
gp_agri <- cor.test(gpara_dat$av_prev, gpara_dat$arableLand)
gp_agri$p.value * 3 ## Bonferoni Correction
gp_tu <- cor.test(gpara_dat$av_prev, gpara_dat$max_sumTUiv)
gp_tu$p.value * 3 ## Bonferoni Correction

###### Section 5 Checking correlations of baetis parasites ######
## Test variables pastures, agriculture and max_sumTUiv
bp_past <- cor.test(bpara_dat$av_prev, bpara_dat$pastures)
bp_past$p.value * 3
bp_agri <- cor.test(bpara_dat$av_prev, bpara_dat$arableLand)
bp_agri$p.value * 3
bp_tu <- cor.test(bpara_dat$av_prev, bpara_dat$max_sumTUiv)
bp_tu$p.value * 3

