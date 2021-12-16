## ------------------------------------ Script for creating index of landscape heterogeneity ----------------- ##

## WORKFLOW
## 1. Download catchments and land use
## 2. Intersect land use with catchments
## 3. Reclassify land use classes to Forest, Pastures, Arable land and Other
## 4. make raster from each polygon and save as file, once for clc18 data and once for lu_class
## 5. load rasters and run analysis of landscape heterogeneity on rasters, 
##    use package landscapemetrics
## 6. check collinearity of calculated metrics and between datasets.

###### Section 1 Load packages and set paths ######
library(data.table)
library(dplyr)
library(sf)
library(RPostgreSQL)
library(stars)
library(raster)
library(tidyr)


wd <- "/home/moritz/Nextcloud/PHD/Projects/Romania/UKL_RO_Macroinv"
## check the wd
getwd()

## create path to data folders
datadir <- '/home/moritz/Nextcloud/PHD/Projects/Romania/UKL_RO_Macroinv/data'

## enable connection
source("/home/moritz/Nextcloud/PHD/Projects/Romania/Daten/Scripts/00_Connection_Data.R")
drv = dbDriver("PostgreSQL")
con = dbConnect(drv, dbname = NEW_DB, user = DBuser_m, host = DBhost, port = DBport, password = DBpassword_m)

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


###### Section 2 Load the CLC data for catchments ######

## //// Load the sites and the catchments
## sampling sites
smpl_sites <- st_read(con, query = 'SELECT * FROM spatial.sites_snapped')
smpl_sites <- smpl_sites[!smpl_sites$site %in% c("X", "Y", "Z"),]

## catchments
catchments <- st_read(con, query = 'SELECT * FROM spatial_derived.catchments')
catchments <- catchments[!catchments$site %in% c("X", "Y", "Z"),]

## download landuse file from server
land_use <- st_read(con, query = 'select * from spatial_derived.corine_land_cover_2018')

## check for catchments with invalid geometries (e.g. self-intersecting polygons)
inval_catch <- catchments[!st_is_valid(catchments$geometry), ]
catchments_valid = NULL
todo <- catchments$site
for(i in seq_along(todo)){
  sub_catch <- catchments[catchments$site == todo[i],]
  if(st_is_valid(sub_catch) == FALSE){
    sub_catch <- st_buffer(sub_catch, 0.1)
  }
  catchments_valid <- rbind(catchments_valid, sub_catch)
}
inval_catch2 <- catchments[!st_is_valid(catchments_valid$geometry), ]
## --> all catchments valid now!

## intersect with catchments
catch_clc18 <- st_intersection(catchments_valid, land_use)

## create new column with reclassified land use properties
catch_clc18$lu_class <- ifelse(catch_clc18$c18 %in% c(311,312,313,324),
                               1,
                               ifelse(catch_clc18$c18 %in% c(211,222),
                                      2,
                                      ifelse(catch_clc18$c18 %in% c(231, 243),
                                             3,
                                             ifelse(! catch_clc18$c18 %in% c(311,312,313,324,211,222,231,243,511,512),
                                                    4,
                                                    NA))))


###### Section 3 Turn Polygons into Raster ######

## run a loop to rasterize each catchments individually
## save two files per catchments, one with clc codes, one with land use codes
todo <- unique(catch_clc18$site)
for(i in seq_along(todo)){
  s <- todo[i]
  catch <- catch_clc18[catch_clc18$site == s,]
  catch_rast <- st_rasterize(catch,
                             dx = 5, dy = 5)
  name1 <- paste0("catchclc_", s, "_raster")
  name2 <- paste0("catchlu_", s, "_raster")
  write_stars(catch_rast, layer = 5, file.path(datadir, "/raster_data", paste0(name1,".tif")))
  write_stars(catch_rast, layer = 6, file.path(datadir, "/raster_data", paste0(name2,".tif")))
}


###### Section 4 Load Raster and run analysis ######

## load packages landscapemetrics and landscapetools
library(landscapemetrics)
library(landscapetools)

## load rasters of with clc codes
raster_files <- list.files("./data/raster_data", pattern = "^catchclc*")
# raster_files <- raster_files[1]
lsm_clc = NULL
pb = txtProgressBar(min = 0, max = length(raster_files), initial = 0) 

## takes a few minutes to run...
for(i in seq_along(raster_files))
{
  setTxtProgressBar(pb,i)
  ## ------------------------- ##
  ## -- read file by file  --- ##
  ## ------------------------- ##
  file <- raster_files[i]
  temp_file <- raster(file.path("./data/raster_data", file))
  
  lsm_measures <- calculate_lsm(temp_file, 
                                what = c("lsm_l_area_mn", ## mean patch area
                                         "lsm_l_frac_cv", ## 'shape metric' gives variation of fractal dimension in an area
                                         "lsm_l_mesh", ## Aggregation metric
                                         "lsm_l_pd", ## patch density in area
                                         "lsm_l_lsi", ## aggregation metric, measures ratio of total edge length to theoretical min edge length
                                         "lsm_l_prd"), ## Patch richness density
                                full_name = TRUE)
  lsm_measures$site <- c(rep(substr(file, start = 10, stop = 10), nrow(lsm_measures)))
  
  lsm_clc <- rbind(lsm_clc, lsm_measures)
}    


## load rasters with land use classes
raster_files2 <- list.files("./data/raster_data", pattern = "^catchlu*")
lsm_lu = NULL
pb = txtProgressBar(min = 0, max = length(raster_files2), initial = 0) 

## takes a few minutes to run
for(i in seq_along(raster_files2))
{
  setTxtProgressBar(pb,i)
  ## ------------------------- ##
  ## -- read file by file  --- ##
  ## ------------------------- ##
  file <- raster_files2[i]
  temp_file <- raster(file.path("./data/raster_data", file))
  
  lsm_measures <- calculate_lsm(temp_file, 
                                what = c("lsm_l_area_mn", ## mean patch area
                                         "lsm_l_frac_cv", ## 'shape metric' gives variation of fractal dimension in an area
                                         "lsm_l_mesh", ## Aggregation metric
                                         "lsm_l_pd", ## patch density in area
                                         "lsm_l_lsi", ## aggregation metric, measures ratio of total edge length to theoretical min edge length
                                         "lsm_l_prd"), ## Patch richness density
                                full_name = TRUE)
  lsm_measures$site <- c(rep(substr(file, start = 9, stop = 9), nrow(lsm_measures)))
  
  lsm_lu <- rbind(lsm_lu, lsm_measures)
}    

## cast lsm_clc and lsm_lu into wide format
lsm_clc <- pivot_wider(lsm_clc[,c(5, 6,10)], names_from = "metric",  values_from = "value")
lsm_lu <- pivot_wider(lsm_lu[,c(5, 6,10)], names_from = "metric",  values_from = "value")

## check collinearity of variables
pairs(lsm_clc[,-1], lower.panel = panel.smooth, upper.panel = panel.cor)
## check collinearity between datasets
par(mfrow = c(2,3))
plot(lsm_clc$area_mn ~ lsm_lu$area_mn)
plot(lsm_clc$frac_cv ~ lsm_lu$frac_cv)
plot(lsm_clc$mesh ~ lsm_lu$mesh)
plot(lsm_clc$pd ~ lsm_lu$pd)
plot(lsm_clc$lsi ~ lsm_lu$lsi)
plot(lsm_clc$prd ~ lsm_lu$prd)
par(mfrow = c(1,1))
## ====> correlation of metrics looks good for different aggregation levels, 
## ====> I will use the landuse measures on aggregated level

## write both dfs to csv
write.csv(lsm_clc, "./data/raster_data/lsm_metrics_clc.csv", row.names = FALSE)
write.csv(lsm_lu, "./data/raster_data/lsm_metrics_lu.csv", row.names = FALSE)
