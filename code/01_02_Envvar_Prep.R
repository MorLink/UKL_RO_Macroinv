## ------------------------------ Script for environmental variable preparation prior to regression analyses and RDA ------------------------ ##

###### Section 1: Load directories and packages ######
library(data.table)
library(dplyr)
library(ggplot2)
require(vegan)

setwd("/home/moritz/Nextcloud/PHD/Projects/Romania/UKL_RO_Macroinv")
## check the wd
getwd()

## create path to data folders
datadir <- '/home/moritz/Nextcloud/PHD/Projects/Romania/UKL_RO_Macroinv/data'


## load functions
## for correlation plot
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


###### Section 2: Load Data and split into April and June data #######
## load environmental data with single measures for all variables
env <- read.csv(file.path(datadir, "environmental_data/env_predictors.csv"), stringsAsFactors = FALSE)

## drop columns land_use_agri (better use corine data) and nat_nonforest
env <- env[ , -which(names(env) %in% c("land_use_agri", "nat_nonforest"))]

## cl = chloride in mg/L
## el_cond = electrical conductivity
## nh4 = ammonium in mg/L
## no2 = Nitrite in mg/L
## no3 = Nitrate in mg/L
## o2_mg = dissolved Oxygen in %
## ph = pH
## po4 = Phosphate in mg/L
## water_temp = water temperature
## flow = flow velocity at sampling site
## max_sumtu_iv = maximum sumTU for most sensitive invertebrate over all sampling events
## low_rip_veg = % area with low riparian vegetation (including bare soil and vegetation below 50cm heigt)
## min_dist_field = minimum distance to field on both sites of stream at sampling site in m
## curvature = degree of stream meandering, increasing values indicate increasing meandering
## riffles = ratio of riffles in stream
## shading = ratio of shaded area at sampling site
## fine_sed = ratio of fine sediment in the stream substrate at the sampling site
## totalarea = area of upstream catchment
## urban = ratio of urban area in catchment, based on corine land cover 2018
## pastures = ratio of pastures in catchment, based on corine land cover 2018
## agriculture = ratio of agriculture in catchment, based on corine land cover 2018
## forest = ratio of forest area in cathcment, based on corine land cover 2018
## rip_urban = ratio of urban area in riparian corridor, based on corine land cover 2018
## rip_forest = ratio of forest area in riparian corridor, based on corine land cover 2018
## rip_past = ratio of pasture area in riparian corridor, based on corine land cover 2018
## rip_agri = ratio of agricultural area in riparian corridor, based on corine land cover 2018
## refugium = presence/absence of upstream forest reach, acting as a refuge for sensitive taxa

## load landscape heterogeneity data
lsm <- read.csv("./data/raster_data/lsm_metrics_lu.csv", stringsAsFactors = FALSE)
## keep only variables that are not redundant
pairs(lsm[,-1], lower.panel = panel.smooth, upper.panel = panel.cor)
## most reasonable choices are pd, prd and mesh(just out of interest)
## lsm <- lsm[,c("site", "mesh", "pd", "prd")]
## merge to env
env <- merge(env, lsm, by = "site", all.x = TRUE)

## Drop April-Data, remove samplinge sites X, Y, Z
env_jn <- env[env$month == "June" & !env$site %in% c("X", "Y", "Z"),]
## add column with actual agricultural area in km^2 per catchment
env_jn$log_agriarea <- log10(env_jn$totalarea * env_jn$agriculture)



###### Section 3: Check distribution of environmental predictors of June ######
## distributions for June data
## we check whether sqrt or log transformation improves the distribution
## for all variables
par(mfrow = c(4, 4))
for (nam in names(env_jn)[c(3:10)]) {
  y <- env_jn[, nam]
  boxplot(y, main = "", xlab = nam)
  qqnorm(y)
  qqline(y)
}

for (nam in names(env_jn)[11:18]) {
  y <- env_jn[, nam]
  boxplot(y, main = "", xlab = nam)
  qqnorm(y)
  qqline(y)
}

for (nam in names(env_jn)[19:26]) {
  y <- env_jn[, nam]
  boxplot(y, main = "", xlab = nam)
  qqnorm(y)
  qqline(y)
}

par(mfrow=c(4,4))
for (nam in names(env_jn)[27:33]) {
  y <- env_jn[, nam]
  boxplot(y, main = "", xlab = nam)
  qqnorm(y)
  qqline(y)
}


par(mfrow = c(4, 4))
for (nam in names(env_jn)[c(3,5)]) {
  y <- env_jn[, nam]
  boxplot((y), main = "", xlab = nam)
  boxplot(log10(y), main = "", xlab = nam)
  boxplot(sqrt(y), main = "", xlab = nam)
  boxplot((y)^0.25, main = "", xlab = nam)
  qqnorm(y)
  qqline(y)
  qqnorm(log10(y))
  qqline(log10(y))
  qqnorm(sqrt(y))
  qqline(sqrt(y))
  qqnorm((y)^0.25)
  qqline((y)^0.25)
}
## no improvement for nh4

for (nam in names(env_jn)[c(6,7)]) {
  y <- env_jn[, nam]
  boxplot((y), main = "", xlab = nam)
  boxplot(log10(y), main = "", xlab = nam)
  boxplot(sqrt(y), main = "", xlab = nam)
  boxplot((y)^0.25, main = "", xlab = nam)
  qqnorm(y)
  qqline(y)
  qqnorm(log10(y))
  qqline(log10(y))
  qqnorm(sqrt(y))
  qqline(sqrt(y))
  qqnorm((y)^0.25)
  qqline((y)^0.25)
}
## no improvement for no2 and no3

for (nam in names(env_jn)[c(11,15)]) {
  y <- env_jn[, nam]
  boxplot((y), main = "", xlab = nam)
  boxplot(log10(y), main = "", xlab = nam)
  boxplot(sqrt(y), main = "", xlab = nam)
  boxplot((y)^0.25, main = "", xlab = nam)
  qqnorm(y)
  qqline(y)
  qqnorm(log10(y))
  qqline(log10(y))
  qqnorm(sqrt(y))
  qqline(sqrt(y))
  qqnorm((y)^0.25)
  qqline((y)^0.25)
}
## log transform min_dist_field


for (nam in names(env_jn)[c(20,21)]) {
  y <- env_jn[, nam]
  boxplot((y), main = "", xlab = nam)
  boxplot(log10(y), main = "", xlab = nam)
  boxplot(sqrt(y), main = "", xlab = nam)
  boxplot((y)^0.25, main = "", xlab = nam)
  qqnorm(y)
  qqline(y)
  qqnorm(log10(y))
  qqline(log10(y))
  qqnorm(sqrt(y))
  qqline(sqrt(y))
  qqnorm((y)^0.25)
  qqline((y)^0.25)
}
## no improvement for totalarea and pastures

par(mfrow=c(2,4))
for (nam in names(env_jn)[c(22)]) {
  y <- env_jn[, nam]
  boxplot((y), main = "", xlab = nam)
  boxplot(log10(y), main = "", xlab = nam)
  boxplot(sqrt(y), main = "", xlab = nam)
  boxplot((y)^0.25, main = "", xlab = nam)
  qqnorm(y)
  qqline(y)
  qqnorm(log10(y))
  qqline(log10(y))
  qqnorm(sqrt(y))
  qqline(sqrt(y))
  qqnorm((y)^0.25)
  qqline((y)^0.25)
}
## no improvement for agriculture

par(mfrow=c(4,3))
for (nam in names(env_jn)[c(14,18)]) {
  y <- env_jn[, nam]
  boxplot((y), main = "", xlab = nam)
  boxplot(sqrt(y), main = "", xlab = nam)
  boxplot((y)^0.25, main = "", xlab = nam)
  qqnorm(y)
  qqline(y)
  qqnorm(sqrt(y))
  qqline(sqrt(y))
  qqnorm((y)^0.25)
  qqline((y)^0.25)
}

par(mfrow=c(2,3))
for (nam in names(env_jn)[c(19)]) {
  y <- env_jn[, nam]
  boxplot((y), main = "", xlab = nam)
  boxplot(sqrt(y), main = "", xlab = nam)
  boxplot((y)^0.25, main = "", xlab = nam)
  qqnorm(y)
  qqline(y)
  qqnorm(sqrt(y))
  qqline(sqrt(y))
  qqnorm((y)^0.25)
  qqline((y)^0.25)
}


## transformation of variables
env_jn$log_mindistf <-log10(env_jn$min_dist_field)

#### --------- subset to variables of interest ------------##
## cl, el_cond, nh4, no3, pH, po4, water_temp, flow, max_sumtu_iv
## riffles, curvature, fine_sed, pastures, agriculture, forest, refugium, log_mindistf
col.names <- c("site", "cl", "el_cond", "po4", "ph", "o2_mg", "flow", "max_sumtu_iv",
               "curvature", "riffles", "fine_sed", "pastures", "water_temp", "agriculture",
               "forest", "refugium", "log_agriarea", "log_mindistf", "lsi", "frac_cv")
env_jn <- env_jn[col.names]

## check correlation plots
pairs(env_jn[,c(2:20)], lower.panel = panel.smooth, upper.panel = panel.cor)
## r > 0.7 for agriculture and forest; forest and refugium; max_sumtu_iv and log_agriarea
## pH and o2_mg
## remove cl, forest, log_agriarea, ph, curvature, riffles
env_jn <- env_jn[ , -which(names(env_jn) %in% c("cl", "forest", "log_agriarea", "ph", 
                                                "curvature", "riffles", "log_mindistf"))]
pairs(env_jn[,c(2:13)], lower.panel = panel.smooth, upper.panel = panel.cor)
## keep only lsi of lsm variables
env_jn <- env_jn[, -which(names(env_jn) %in% c("frac_cv"))]

## check vif
library(faraway)
sort(vif(env_jn[,-c(1)]))
## remove water_temp
env_jn <- env_jn[, -which(names(env_jn) %in% c("water_temp"))]
sort(vif(env_jn[,-c(1)]))
## all vif under 7
names(env_jn)[2:11] <- c("elCond", "orthoP", "mgO2", "flow", "max_sumTUiv", "fineSed", 
                         "pastures", "arableLand", "FUS", "ls_hetero")

write.csv(env_jn, file.path(datadir, "environmental_data/env_predictors_june.csv"), row.names = FALSE)

