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



## pairs plot function
panel.cor <- function(x, y, digits=2, prefix="", cex.cor) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, use = "complete.obs"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste(prefix, txt, sep = "")
  if (missing(cex.cor)) cex <- 1 / strwidth(txt)
  
  test <- cor.test(x, y)
  # borrowed from printCoefmat
  Signif <- symnum(
    test$p.value, corr = FALSE, na = FALSE,
    cutpoints = c(0, 0.05, 0.1, 1),
    symbols = c("*", ".", " ")
  )
  
  text(0.5, 0.5, txt, cex = cex * r)
  text(5, 5, Signif, cex = cex, col = 2)
}

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

pairs(gam_dat[,-1], lower.panel = panel.smooth, upper.panel = panel.cor)
## gammarid body weight decreases with average prevalence

## merge environmental data with gam_dat
gpara_dat <- merge(gam_dat[,c(1,2,3,8)], env_jn, by = "site")
pairs(gpara_dat[,-c(1,2,3)], lower.panel = panel.smooth, upper.panel = panel.cor)

## merge environmental data with bea_par
bae_par$av_prev <- round((bae_par$microsporidia + bae_par$nematoda + bae_par$trematoda)/3, 1)
bpara_dat <- merge(bae_par[,-c(2:5)], env_jn, by = "site")
pairs(bpara_dat[,-1], lower.panel = panel.smooth, upper.panel = panel.cor)
## microsporidia seem to be affected by po4

###### Section 4 model gammarid data with environmental parameters ######
## only 12 observations, too few observations for glmnet and cross validation
## check with basic correlation

## correlation of average prevalence with environmental variables
cor(gpara_dat[,-c(1:3)])[1,]

## Test variables pastures, agriculture and max_sumTUiv
gp_past <- cor.test(gpara_dat$av_prev, gpara_dat$pastures)
gp_past$p.value * 3 ## Bonferoni Correction
# gp_fus <- cor.test(gpara_dat$av_prev, gpara_dat$FUS)
# gp_fus$p.value * 4 ## Bonferoni Correction
gp_agri <- cor.test(gpara_dat$av_prev, gpara_dat$arableLand)
gp_agri$p.value * 3 ## Bonferoni Correction
gp_tu <- cor.test(gpara_dat$av_prev, gpara_dat$max_sumTUiv)
gp_tu$p.value * 3 ## Bonferoni Correction

###### Section 5 Checking correlations of baetis parasites ######
cor(bpara_dat[,-1])[1,]

## Test variables pastures, agriculture and max_sumTUiv
bp_past <- cor.test(bpara_dat$av_prev, bpara_dat$pastures)
bp_past$p.value * 3
bp_agri <- cor.test(bpara_dat$av_prev, bpara_dat$arableLand)
bp_agri$p.value * 3
bp_tu <- cor.test(bpara_dat$av_prev, bpara_dat$max_sumTUiv)
bp_tu$p.value * 3

