## ------------------------- Calculate macroinvertebrate community metrics ----------------------------- ##
## This script calculates:
## - total number of taxa Done
## - number of insect taxa Done
## - % EPT taxa Done
## - Shannon taxa diversity Done
## - Functional richness, divergence and evenness, based on traits body size, feeding, locomotion and aquatic stages


########### Section 1 Load packages and set paths ###########
## ----------------------------------- Load packages ----------------------------------- ##
## packages
library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)
library(vegan)
library(readxl)

## -------------------------------------- load connections and paths ------------------------------ ##
wd <- setwd("/home/moritz/Nextcloud/PHD/Projects/Romania/UKL_RO_Macroinv")
## check the wd
getwd()

## create path to data folders
datadir <- '/home/moritz/Nextcloud/PHD/Projects/Romania/UKL_RO_Macroinv/data'



########## Section 2 Load invertebrate data and compute metrics###########

## //////////////////////////////////////////////////
## ---------- Load macroinvertebrate data -----------
## //////////////////////////////////////////////////
## load invertebrate data
inv <- read.csv(file.path(datadir, "invertebrate_data/macroinv_full_ap_jn.csv"), stringsAsFactors = FALSE, na.strings=c("","NA"))

setDT(inv)

## keep only data from June and April
mzb_mi <- inv[date >= "2016-06-01" | date <= "2016-04-30",,]
mzb_mi$month <- ifelse(mzb_mi$date <= "2016-04-30",
                       paste0("April"),
                       ifelse(mzb_mi$date >= "2016-06-01",
                              paste0("June"),
                              NA))

## add column with lowest identification level
mzb_mi$classification <- ifelse(!is.na(mzb_mi$species),
                                      paste0(mzb_mi$species),
                                      ifelse(!is.na(mzb_mi$genus),
                                             paste0(mzb_mi$genus),
                                             ifelse(!is.na(mzb_mi$family),
                                                    paste0(mzb_mi$family),
                                                    ifelse(!is.na(mzb_mi$order),
                                                           paste0(mzb_mi$order),
                                                           ifelse(!is.na(mzb_mi$sub_class),
                                                                  paste0(mzb_mi$sub_class),
                                                                  ifelse(!is.na(mzb_mi$class),
                                                                         paste0(mzb_mi$class),
                                                                         ifelse(!is.na(mzb_mi$phylum),
                                                                                paste0(mzb_mi$phylum),
                                                                                NA)))))))

## remove micro crustaceans and water mites 
mzb_mi <- mzb_mi[!classification %in% c("Copepoda", "Ostracoda", "Trombidiformes"),]
## remove O.DOWN sites from datset
mzb_mi <- mzb_mi[!site == "O.DOWN",]
## rename site O.UP
mzb_mi[site == "O.UP",]$site <- c(rep("O", 18))

## //////////////////////////////////////////////////
## ---------remove pseudo taxa from mzb data -------
## //////////////////////////////////////////////////
## remove pseudo taxa like Baetis sp. from dataset
## most community metrics are based on taxa number,
## so loss of abundance data is neglectable

sites_set <- unique(mzb_mi$site)
mzb_taxanmbr <- NULL
# sites_set <- sites_set[18]

## loop on sites
for(i in seq_along(sites_set))
{
  x <- sites_set[i]
  single_site <- mzb_mi[site == x,,]
  
  ## loop on months
  month_set <- unique(single_site$month)
  # month_set <- month_set[1]
  mzb_month <- NULL
  for (j in seq_along(month_set))
  {
    pset <- NULL
    y <- month_set[j]
    site_month <- single_site[month == y,,]
    
    ## genus level;
    ## loop through genus level, check if species columns contains NAs AND species;
    ## compares nrow of genera WITH NA entries in species to nrow with actual species;
    ## if there are species AND NA entries, remove the NA entries
    genus_set <- unique(site_month[!is.na(genus),]$genus)
    if(length(genus_set) > 0)
    {
      mzb_genlev <- NULL
      for(k in seq_along(genus_set))
      {
        gen<- genus_set[k]
        gen_group <- site_month[genus == gen,,]
        na_gen <- gen_group[is.na(species),,]
        if(nrow(na_gen) > 0 & nrow(na_gen) < nrow(gen_group))
        {gen_group <- gen_group[!is.na(species),,]}
        mzb_genlev <- rbind(mzb_genlev, gen_group)
      }
      ## write rows to combined file
      pset <- rbind(pset, mzb_genlev)
    }
    
    ## family level
    ## drop rows from mzb_genlev from site_month
    mzb_sub1 <- site_month[!family %in% unique(pset[!is.na(family),,]$family),,]
    family_set <- unique(mzb_sub1[!is.na(family),]$family)
    if(length(family_set) > 0)
    {
      mzb_famlev <- NULL
      for(l in seq_along(family_set))
      {
        fam<- family_set[l]
        fam_group <- mzb_sub1[family == fam,,]
        na_fam <- fam_group[is.na(genus),,]
        if(nrow(na_fam) > 0 & nrow(na_fam) < nrow(fam_group))
        {fam_group <- fam_group[!is.na(genus),,]}
        mzb_famlev <- rbind(mzb_famlev, fam_group)
      }
      ## write rows to combined file
      pset <- rbind(pset, mzb_famlev)
    }
    
    ## order level
    mzb_sub2 <- mzb_sub1[!order %in% unique(pset[!is.na(order),,]$order),,]
    order_set <- unique(mzb_sub2[!is.na(order),]$order)
    if(length(order_set) > 0)
    {
      mzb_ordlev <- NULL
      for(m in seq_along(order_set))
      {
        ord<- order_set[m]
        ord_group <- mzb_sub2[order == ord,,]
        na_ord <- ord_group[is.na(family),,]
        if(nrow(na_ord) > 0 & nrow(na_ord) < nrow(ord_group))
        {ord_group <- ord_group[!is.na(family),,]}
        mzb_ordlev <- rbind(mzb_ordlev, ord_group)
      }
      ## write rows to combined file
      pset <- rbind(pset, mzb_ordlev)
    }
    
    ## sub_class level
    mzb_sub3 <- mzb_sub2[!sub_class %in% unique(pset[!is.na(sub_class),,]$sub_class)]
    sub_class_set <- unique(mzb_sub3[!is.na(sub_class),]$sub_class)
    if(length(sub_class_set) > 0)
    {
      mzb_subclalev <- NULL
      for(n in seq_along(sub_class_set))
      {
        subcla<- sub_class_set[n]
        subcla_group <- mzb_sub3[sub_class == subcla,,]
        na_subcla <- subcla_group[is.na(order),,]
        if(nrow(na_subcla) > 0 & nrow(na_subcla) < nrow(subcla_group))
        {subcla_group <- subcla_group[!is.na(order),,]}
        mzb_subclalev <- rbind(mzb_subclalev, subcla_group)
      }
      ## write rows to combined file
      pset <- rbind(pset, mzb_subclalev)
    }
    
    ## class level
    mzb_sub4 <- mzb_sub3[!class %in% unique(pset[!is.na(class),,]$class)]
    class_set <- unique(mzb_sub4[!is.na(class),]$class)
    if(length(class_set) > 0)
    {
      mzb_clalev <- NULL
      for(o in seq_along(class_set))
      {
        cla<- class_set[o]
        cla_group <- mzb_sub4[class == cla,,]
        na_cla <- cla_group[is.na(sub_class) & is.na(order),,]
        if(nrow(na_cla) > 0 & nrow(na_cla) < nrow(cla_group))
        {cla_group <- cla_group[!is.na(sub_class) & !is.na(order),,]}
        mzb_clalev <- rbind(mzb_clalev, cla_group)
      }
      ## write rows to combined file
      pset <- rbind(pset, mzb_clalev)
    }
    
    ## phylum level 
    mzb_sub5 <- mzb_sub4[!phylum %in% unique(pset$phylum)]
    phylum_set <- unique(mzb_sub5[!is.na(phylum),]$phylum)
    if(length(phylum_set) > 0)
    {
      mzb_phylev <- NULL
      for(p in seq_along(phylum_set))
      {
        phy<- phylum_set[p]
        phy_group <- site_month[phylum == phy,,]
        na_phy <- phy_group[is.na(class)]
        if(nrow(na_phy) > 0 & nrow(na_phy) < nrow(phy_group))
        {phy_group <- phy_group[!is.na(class),,]}
        mzb_phylev <- rbind(mzb_phylev, phy_group)
      }
      ## write rows to combined file
      pset <- rbind(pset, mzb_phylev)
    }
    mzb_month <- rbind(mzb_month, pset)
  }
  mzb_taxanmbr <- rbind(mzb_taxanmbr, mzb_month)
}


##////////////////////////////////////////////////////
## ----------- Number of Taxa ----------------------- 
##////////////////////////////////////////////////////
## richness with mzb_taxrich
## number of species per sqrt(total number of individual)
taxanmbr <- mzb_taxanmbr[, .(tax_nmbr = .N), by = .(site, month)]

##////////////////////////////////////////////////////
## ----------- Number of insect taxa  --------------- 
##////////////////////////////////////////////////////
## subset mzb_taxanmbr to insects
insect_taxanmbr <- mzb_taxanmbr[class == "Insecta", .(insect_taxanmbr = .N), by = .(site, month)]

## merge taxanmbr and insect_taxanmbr
metrics_dt <- merge.data.table(taxanmbr, insect_taxanmbr, by = c("site", "month"), all = TRUE)

##////////////////////////////////////////////////////
## ----- Shannon and Simpson Taxa Diversity  -------
##////////////////////////////////////////////////////
## reduce data table to get an abundance_dt
mzb_tn_ap <- mzb_taxanmbr[month == "April", c(1,13,14,11)]
mzb_tn_jn <- mzb_taxanmbr[month == "June", c(1,13,14,11)]

## cast into wide format, column names are the species names,
abun_wide_ap <- dcast.data.table(mzb_tn_ap, site ~ classification, value.var = "abundance")
abun_wide_jn <- dcast.data.table(mzb_tn_jn, site ~ classification, value.var = "abundance")
## turn all NA into 0
abun_wide_ap[is.na(abun_wide_ap)] <- 0
abun_wide_jn[is.na(abun_wide_jn)] <- 0
## drop site-column
abun_wide_ap <- abun_wide_ap[,-1]
abun_wide_jn <- abun_wide_jn[,-1]
names(abun_wide_ap)[1:91] <- c(paste0("sp_",rep(1:length(abun_wide_ap))))
names(abun_wide_jn)[1:66] <- c(paste0("sp_",rep(1:length(abun_wide_jn))))

## calculate shannon and simpson diversities
swd_ap_cl <- diversity(abun_wide_ap, index = "shannon", MARGIN = 1)
swd_jn_cl <- diversity(abun_wide_jn, index = "shannon", MARGIN = 1)
sd_ap_cl <- diversity(abun_wide_ap, index = "invsimpson", MARGIN = 1)
sd_jn_cl <- diversity(abun_wide_jn, index = "invsimpson", MARGIN = 1)
# relationship between taxonomic diversities (Total taxonomic richness, Shannon-Wiener Diversity and Simpson Diversity) and gradient
tax_div_ap <- data.frame(swd_ap_cl = swd_ap_cl,
                         sd_ap_cl = sd_ap_cl,
                         stringsAsFactors = FALSE)
tax_div_jn <- data.frame(swd_jn_cl = swd_jn_cl,
                         sd_jn_cl = sd_jn_cl,
                         stringsAsFactors = FALSE)
tax_div_ap$site <- sort(unique(mzb_taxanmbr[month == "April",,]$site))
tax_div_jn$site <- sort(unique(mzb_taxanmbr[month == "June",,]$site))
tax_div_ap$month <- c(rep("April", nrow(tax_div_ap)))
tax_div_jn$month <- c(rep("June", nrow(tax_div_jn)))
names(tax_div_ap) <- c("swd", "sd", "site", "month")
names(tax_div_jn) <- c("swd", "sd", "site", "month")
## rbind tax_div_ap and tax_div_jn
tax_div <- rbind(tax_div_ap, tax_div_jn)
tax_div$swd <- round(tax_div$swd, 3)
tax_div$sd <- round(tax_div$sd, 3)

## merge tax_div with metrics_dt
metrics_dt <- merge.data.table(metrics_dt, tax_div, by = c("site", "month"), all.x = TRUE)

##////////////////////////////////////////////////////
## ----------  % EPT-taxa --------------
##////////////////////////////////////////////////////
ept_taxanmbr <- mzb_taxanmbr[order %in% c("Trichoptera", "Plecoptera", "Ephemeroptera"), 
                             .(ept_taxanmbr = .N), by = .(site, month)]
## merge with taxanmbr
metrics_dt <- merge.data.table(metrics_dt, ept_taxanmbr, by = c("site", "month"), all.x = TRUE)
metrics_dt[site == "R" & month == "April",]$ept_taxanmbr <- 0
metrics_dt$epttax_percent <- round((metrics_dt$ept_taxanmbr / metrics_dt$tax_nmbr)*100, 2)


##////////////////////////////////////////////////////
## ----------  % EPT-Abundance --------------
##////////////////////////////////////////////////////
ept_abun <- mzb_mi[order %in% c("Trichoptera", "Plecoptera", "Ephemeroptera"), 
                             .(ept_abun = sum(abundance)), by = .(site, month)]
tot_abun <- mzb_mi[, .(tot_abun = sum(abundance)), by = .(site,month)]
ept_abun <- merge.data.table(ept_abun, tot_abun, by = c("site", "month"))
ept_abun$ept_percent <- round((ept_abun$ept_abun / ept_abun$tot_abun)*100, 2)

## merge with metrics_dt
metrics_dt <- merge.data.table(metrics_dt, ept_abun[,c(1,2,5)], by = c("site", "month"), all.x = TRUE)
metrics_dt[site == "R" & month == "April",]$ept_percent <- 0

##////////////////////////////////////////////////////
## ----- continue with % Log10(SEL_EPTD+1) ---------- 
##////////////////////////////////////////////////////
## log10(sum of individuals within familiy Heptageniidae, //Ephemeridae,
## Leptophlebiidae, //Brachycentridae, //Georidae, Polycentropodidae, Limnephilidae,
## //Odontoceridae, Dolichopodidae, Stratiomyidae, Dixidae, Empididae, 
## //Athericidae, //Nemouridae + 1)
sel <- c("Heptageniidae", "Ephemeridae", "Leptophlebiidae", "Brachycentridae", "Georidae",
         "Polycentropodidae", "Limnephilidae", "Odontoceridae", "Dolichopodidae", "Stratiomyidae",
         "Dixidae", "Empidiae", "Athericidae", "Nemouridae")
unique(mzb_mi[family %in% sel,][,8]) ## nine of the families are in our data
sel_eptd <- mzb_mi[family %in% sel,
                    .(log_septd = log10(sum(abundance)+1)),by=.(site, month)]
sel_eptd$log_septd <- round(sel_eptd$log_septd, 3)
## merge with metrics_dt
metrics_dt <- merge.data.table(metrics_dt, sel_eptd, by = c("site", "month"), all.x = TRUE)
metrics_dt[is.na(log_septd),]$log_septd <- c(0,0)

##////////////////////////////////////////////////////
## ------------ continue with 1-GOLD ----------------- 
##////////////////////////////////////////////////////
## 1 - ratio of Gastropoda + Oligochaeta + Diptera
gold <- mzb_mi[class %in% c("Gastropoda", "Clitellata") | order %like% "Diptera",
                    sum(abundance),by= .(site, month)]
tot_abun <- mzb_mi[!class %in% c("Gastropoda", "Clitellata") | !order %like% "Diptera",
                   sum(abundance), by= .(site, month)]
gold <- merge(gold, tot_abun, by = c("site", "month"), all = TRUE)
gold$gold <- round(gold$V1.x / gold$V1.y, 3)
## merge with mmi_dt
metrics_dt <- merge.data.table(metrics_dt, gold[,c(1,2,5)], by = c("site", "month"), all.x = TRUE)


###### Section 3 load SPEAR-values and saprobic values ######
##////////////////////////////////////////////////////
## ------------- calculate SPEAR values -------------
##////////////////////////////////////////////////////

## load taxa/spear data
spear_dat <- read.csv(file.path(datadir, "invertebrate_data/spear_taxa.csv"), stringsAsFactors = FALSE)
## change to data.table
setDT(spear_dat)

## numerator
spear_dat$sp_up <- log10(spear_dat$abundance + 1) * spear_dat$spear_class
## denominator
spear_dat$sp_down <- log10(spear_dat$abundance + 1)

## spear value
spear_fin <- spear_dat[, .(spear_pest = sum(sp_up)/sum(sp_down)), by = .(site, month)]
spear_fin$spear_pest <- round(spear_fin$spear_pest*100, 2)

## drop april
sp_jm <- spear_fin[spear_fin$month != "April",]
sp_jm[order(sp_jm$spear_pest),]


## merge with metrics_dt
metrics_dt <- merge.data.table(metrics_dt, spear_fin, by = c("site", "month"), all.x = TRUE)

##////////////////////////////////////////////////////
## ------------- export table -------------
##////////////////////////////////////////////////////
write.csv(metrics_dt, file.path(datadir, "invertebrate_data/macroinv_metrics.csv"), row.names = FALSE)
