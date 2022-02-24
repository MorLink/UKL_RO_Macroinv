## ------------------------- Plots and tables for the Paper ----------------------- ##

## Plots for the Paper
## 1. Map of sampling area
## 2. Abundance of macroinvertebrates as heatmap
## 3. Boxplot of catchment size and presence of refugia
## 4. Scatterplot of max sumTU with log_agriarea
## 5. Scatterplot of max_sumTU with (log?)catchment size
## 6. Boxplot of SPEAR and presence of refugium
## 7. RDA plots of April an June
## 

###### Section 1 ######
library(data.table)
library(dplyr)
library(ggplot2)
require(vegan)
library(reshape2)
library(viridis)

## load connections
wd <- setwd("/home/moritz/Nextcloud/PHD/Projects/Romania/UKL_RO_Macroinv")
## check the wd
getwd()

## create path to data folders
datadir <- '/home/moritz/Nextcloud/PHD/Projects/Romania/UKL_RO_Macroinv/data'


###### Section 2 Plotting Abundance as Heatmap ######

## Macroinvertebrate data in long format
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
                              paste0(mzb_dat$genus, " sp."),
                              ifelse(!is.na(mzb_dat$family),
                                     paste0(mzb_dat$family),
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
mzb_april <- mzb_dat[month == "April"]

## transform Abundance
mzb_dat$log_abun <- log10(mzb_dat$abundance +1)
mzb_june$log_abun <- log10(mzb_june$abundance +1)

## Plotting April-heat map, sites ordered with decreasing toxicity
# Color Brewer palette
library(hrbrthemes)
heatmap <- ggplot(mzb_june, aes(site, taxa, fill= log_abun)) + 
  geom_tile() +
  #scale_y_discrete(limits=rev) +
  scale_fill_distiller(palette = "Oranges", name = "Abundance",
                       direction = "reverse", breaks = c(0.3,0.7,1,1.3,1.7,2,2.3,2.7,3,3.3),
                       labels = c("2","5","10","20","50","100","200","500","1000","2000")) +
  xlab("Sampling sites") + ylab("") +
  #facet_grid(cols = vars(month)) +
  theme_bw() +
  theme(#legend.position="bottom",
        strip.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size=6),
        strip.background = element_rect(color="black", fill="white", linetype = "blank"),
        panel.grid = element_blank(),
        legend.text = element_text(angle = 0)
        #plot.margin=grid::unit(c(2,2,0,0), "mm")
        )

heatmap

ggsave(file= "./plots/Heatmap.tiff", plot=heatmap, width = 6, height = 8)
ggsave(file= "./plots/Heatmap.png", plot=heatmap, width = 6, height = 8)
ggsave(file= "./plots/Heatmap.svg", plot=heatmap, width = 6, height = 8)



## ------------- some simple stats ----------- ##
## how many taxa in April
unique(mzb_april$taxa)
sum(mzb_april$abundance)

## how many taxa in June
unique(mzb_june$taxa)
sum(mzb_june$abundance)

###### Section 3 Plotting Map of sampling sites ######
## -------------------- load packages and connection
library(sf)
library(RColorBrewer)
library(maps)
library(RPostgreSQL)
library(ggspatial)
library(rnaturalearth)
library(svglite)

## enable connection
source("/home/moritz/Nextcloud/PHD/Projects/Romania/Daten/Scripts/00_Connection_Data.R")
drv = dbDriver("PostgreSQL")
con = dbConnect(drv, dbname = NEW_DB, user = DBuser_m, host = DBhost, port = DBport, password = DBpassword_m)


## -------------------- load spatial points, streams and catchments
## sampling sites
smpl_sites <- st_read(con, query = 'SELECT * FROM spatial.sites_snapped')
smpl_sites <- smpl_sites[!smpl_sites$site %in% c("X", "Y", "Z"),]

## catchments
catchments <- st_read(con, query = 'SELECT * FROM spatial_derived.catchments')
catchments <- catchments[!catchments$site %in% c("X", "Y", "Z"),]

## /////////////////////////////////////////////////////////////////
## ///////////// data for the background coloring of the map ///////
## /////////////////////////////////////////////////////////////////
## -------------------- create rectangle for land use of whole map
## create rectangle for sampling area
rect = st_sfc(st_polygon(list(cbind(c(5297691,5364916,5364916,5297691,5297691),
                                   c(2801019,2801019,2699119,2699119,2801019)))))
rect <- st_set_crs(rect, 3035)
rect_b <- st_buffer(rect, 50000)
## upload to database
dbSendQuery(con, "drop table IF EXISTS spatial_derived.sampling_area")
st_write(rect_b, con, layer = c('spatial_derived', 'sampling_area'))
### Add comment on units to columns
dbSendQuery(con, sprintf("
            COMMENT ON TABLE spatial_derived.sampling_area
              IS 'Rectange around the sampling area,
              CRS is EPSG:3035';
            ALTER TABLE spatial_derived.sampling_area
              OWNER to %s;
            GRANT ALL ON TABLE spatial_derived.sampling_area TO %s;
            GRANT ALL ON TABLE spatial_derived.sampling_area TO bfg;
            GRANT SELECT ON TABLE spatial_derived.sampling_area TO bfg_read;",
                         DBuser_m, DBuser_m))

## intersect rectangle of sampling area with corine land cover 2018 on the server
rect_clc18 <- st_read(con, query = "SELECT c18,
                                 st_multi(st_intersection(clc.geometry, rect.geom))
                                 As geometry
                                 FROM spatial_derived.corine_land_cover_2018 AS clc
 	                               INNER JOIN spatial_derived.sampling_area as rect
 	                               ON st_intersects(clc.geometry, rect.geom)")
## group the land use categories
forest_all <- rect_clc18[rect_clc18$c18 %in% c(311,312,313,324),]
agri_all <- rect_clc18[rect_clc18$c18 %in% c(211,222),]
pastures_all <- rect_clc18[rect_clc18$c18 %in% c(231, 243),]
water_all <- rect_clc18[rect_clc18$c18 %in% c(511,512),]
other_all <- rect_clc18[! rect_clc18$c18 %in% c(311,312,313,324,211,222,231,243,511,512),]

## combine land use in rectanlge into one object
agri_all$type <- c(rep("Agriculture", nrow(agri_all)))
forest_all$type <- c(rep("Forest", nrow(forest_all)))
pastures_all$type <- c(rep("Pastures", nrow(pastures_all)))
water_all$type <- c(rep("Streams", nrow(water_all)))
other_all$type <- c(rep("Other", nrow(other_all)))

land_use_all <- rbind(agri_all, forest_all, pastures_all, water_all, other_all)
## Factorize the type column (Agriculture, Pastures, Forest, Other)
land_use_all$type <- factor(land_use_all$type, levels = c("Agriculture", "Pastures", "Forest", "Streams", "Other"))

## intersect rectangle of sampling area with stream network of Romania on the server
streams_ro <- st_read(con, query = "SELECT cat_, stream, prev_str01, prev_str02, length, cum_length,
                                           out_dist, next_str, flow_accu,
                                       st_multi(st_intersection(str.geometry, rect.geom))
                                      As geometry
                                FROM spatial_derived.stream_network AS str
	                              INNER JOIN spatial_derived.sampling_area as rect
	                              ON st_intersects(str.geometry, rect.geom)")
streams_ro$dummy_col <- c(rep(1,nrow(streams_ro)))



## /////////////////////////////////////////////////////////////////
## ///////////// data for the coloring of the catchments     ///////
## /////////////////////////////////////////////////////////////////
## -------------------- intersect the stream network with the calculated catchments
## check size of files
dbGetQuery(con, paste0("SELECT pg_size_pretty(pg_total_relation_size('spatial_derived.stream_network'))"))

## run intersectioin of catchments and streams on server
streams_catch <- st_read(con, query = "SELECT cat_, stream, prev_str01, prev_str02, length, cum_length,
                                           out_dist, next_str, flow_accu,
                                       st_multi(st_intersection(str.geometry, catch.geometry))
                                      As geometry
                                FROM spatial_derived.stream_network AS str
	                              INNER JOIN spatial_derived.catchments
	                              as catch ON st_intersects(str.geometry, catch.geometry)
                                WHERE site NOT IN ('X', 'Y', 'Z')")
streams_catch$dummy_col <- c(rep(1,nrow(streams_catch)))


## -------------------- get land use in catchment
## run intersection of CLC18 and catchments on server
catch_clc18 <- st_read(con, query = "SELECT c18,
                                 st_multi(st_intersection(clc.geometry, catch.geometry))
                                 As geometry
                                 FROM spatial_derived.corine_land_cover_2018 AS clc
 	                               INNER JOIN spatial_derived.catchments as catch
 	                               ON st_intersects(clc.geometry, catch.geometry)
                                 WHERE site NOT IN ('X', 'Y', 'Z')")


forest <- catch_clc18[catch_clc18$c18 %in% c(311,312,313,324),]
agri <- catch_clc18[catch_clc18$c18 %in% c(211,222),]
pastures <- catch_clc18[catch_clc18$c18 %in% c(231, 243),]
water <- catch_clc18[catch_clc18$c18 %in% c(511,512),]
other <- catch_clc18[!catch_clc18$c18 %in% c(311,312,313,324,211,222,231,243,511,512),]

## combine landuse in catchments into one object
agri$type <- c(rep("Agriculture", nrow(agri)))
forest$type <- c(rep("Forest", nrow(forest)))
pastures$type <- c(rep("Pastures", nrow(pastures)))
water$type <- c(rep("Streams", nrow(water)))
other$type <- c(rep("Other", nrow(other)))

land_use <- rbind(agri, forest, pastures, water, other)
## Factorize the type column (Agriculture, Pastures, Forest, Other)
land_use$type <- factor(land_use$type, levels = c("Agriculture", "Pastures", "Forest", "Streams", "Other"))


# Local file location to save to
romania_rds <- "./data/spatial/gadm36_ROU_1_sf.rds"

# download.file(
#   url = paste0(u_remote, p_remote, f_name),
#   destfile = romania_rds,
#   method = "auto")

romania_sf <- readRDS(romania_rds)
romania_sf <- st_transform(romania_sf, 3035)

pol = st_sfc(st_polygon(list(cbind(c(5297691,5364916,5364916,5297691,5297691),
                                   c(2801019,2801019,2699119,2699119,2801019)))))
pol <- st_set_crs(pol, 3035)
pol_b <- st_buffer(pol, 50000)

romania_base <- st_intersection(romania_sf, pol_b)

## get coordinates for cities Cluj Napoca, Turda, Zalau, Jibou
cities <- data.frame(name = c("Cluj-Napoca ", "Zalău", "Jibou"),
                     x = c(23.59228,23.05550,23.25638),
                     y = c(46.77103,47.18216,47.26123),
                     stringsAsFactors = FALSE)
cities_sf <- st_as_sf(cities, coords = c("x", "y"),
                      crs = 4326)

cities_sf <- st_transform(cities_sf, crs = 3035)

smpl_points <- data.frame(st_coordinates(smpl_sites))


## set colors
discrete_palettes <- c("tan3", "yellow2", "palegreen3", "deepskyblue", "gray80")

catchment_plot <-  
  ggplot() +
  geom_sf(data = romania_base, color = "gray98", fill = "gray98") +
  # geom_sf(data = land_use_map, aes(fill = type), color = NA) +
  # geom_sf(data = catch_buf, color = NA, fill = "gray80") + # fills the catch_buf with same color as "Other"
  geom_sf(data = catchments, fill = "gray80") + # fills the catchment with same color as "Other"
  geom_sf(data = land_use_all, aes(fill = type), color = NA, alpha = 0.4) +
  geom_sf(data = land_use, aes(fill = type), color = NA, alpha = 0.6) +
  scale_fill_discrete(type = discrete_palettes, labels = c('Arable land ', 'Pastures', 'Forest', 'Streams', 'Other')) +
  geom_sf(data = catchments, fill = NA) + # for the lines around the catchments
  geom_sf(data = streams_ro, color = "deepskyblue") + 
  geom_point(data = smpl_points, aes(x = X, y = Y), size = 4, color = 'red', alpha = 0.9) +
  geom_point(data = smpl_points, aes(x = X, y = Y), size  = 4, shape = 1,  color = 'red') +
  geom_sf_text(data = smpl_sites,
               aes(label = site),
               size = 3) +
  geom_sf(data = cities_sf, color = "black", size = 2) +
  geom_sf_text(data = cities_sf[2:3,],
               aes(label = name),
               size = 4.5,
               hjust = 1.2) +
  geom_sf_text(data = cities_sf[1,],
               aes(label = name),
               size = 4.5,
               hjust = 1.1) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.35, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(5297691, 5364916), ylim = c(2699119, 2801019), crs = 3035, expand = TRUE) +
  #guides(color=guide_legend("Species"), fill = FALSE) + 
  theme_bw() + theme(axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.text.y = element_blank(),
                     axis.ticks.y = element_blank(),
                     legend.position = c(0.8965, 0.695),
                     legend.background = element_rect(fill = "white", 
                                                      color = "black"),
                     #panel.spacing = unit(c(0,0,0,0), "cm"), 
                     plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")
  ) +
  labs(x = "", y = "", #title = "Study area", caption = "Data from GADM, https://gadm.org",
       fill = "Land use", alpha = 1)

#catchment_plot

## ----------------- Create main map of Romania ------------- ##
## load level 0 map of romania
romania_all <- "./data/spatial/gadm36_ROU_0_sf.rds"

romania_all <- readRDS(romania_all)
romania_all <- st_transform(romania_all, 3035)

## create rectangle polygon for the map
df <- data.frame(
  x = c(5297691, 5297691, 5364916, 5364916),
  y = c(2699119, 2801019, 2801019, 2699119)
)
library(tidyverse)
#> Linking to GEOS 3.6.1, GDAL 2.2.3, proj.4 4.9.3
polygon <- df %>%
  st_as_sf(coords = c("x", "y"), crs = 3035) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")
polygon

## creat sf object for the Romania label
rmn <- data.frame(name = c("Romania"),
                     x = c(24.924202),
                     y = c(45.700000),
                     stringsAsFactors = FALSE)
rmn_sf <- st_as_sf(rmn, coords = c("x", "y"),
                      crs = 4326)
rmn_sf <- st_transform(rmn_sf, crs = 3035)

## create box around the romania map
rmn_box <- data.frame(x = c(5870000, 5870000, 5100000, 5100000),
                      y = c(2370000, 2950000, 2950000, 2370000))
rmn_pol<- rmn_box %>%
  st_as_sf(coords = c("x", "y"), crs = 3035) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")
# rmn_pol <- st_transform(rmn_pol, crs = 3035)


romania_plot <- ggplot() +
  geom_sf(data = rmn_pol, color = "black", fill = "white") +
  geom_sf(data = romania_all,  color = "black", fill = "gray98") + ## map
  geom_sf(data = polygon, color = "black", fill = NA) + ## xlim = c(5297691, 5364916), ylim = c(2699119, 2801019), crs = 3035 use this for rectangle
  geom_sf(data = cities_sf[1,], color = "black", size = 0.8) +
  # geom_sf(data = rmn_pol, color = "black") +
  geom_sf_text(data = cities_sf[1,],
               aes(label = name),
               size = 4.0,
               hjust = -0.07) + ## city name
  geom_sf_text(data = rmn_sf[1,],
               aes(label = name),
               size = 5.5) + ## Country label
  theme_minimal() + theme(axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.text.y = element_blank(),
                     axis.ticks.y = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank()) +
  labs(x = "", y = "")

romania_plot

library(grid)
rect <- rectGrob(
  width = unit(2.35, "in"),
  height = unit(1.8, "in"),
  gp = gpar(fill = "white")
)

## add both plots ontop of each other
library(cowplot)
catchment_fin <- ggdraw(catchment_plot) +
    # draw_grob(
    #   {rect},
    #   x = 0.08,
    #   y = 0.3715) +
    draw_plot(
    {
      romania_plot},
    x = 0.5302,
    y = 0.7102,
    width = 0.44,
    height = 0.32
  )

catchment_fin

ggsave(file= "./plots/Map.svg", plot=catchment_fin, width = 6.5, height = 9)
ggsave(file= "./plots/Map.png", plot=catchment_fin, width = 6.5, height = 9)
ggsave(file= "./plots/Map.tiff", plot=catchment_fin, width = 6.5, height = 9)




###### Section 4 Plot Effect of Refugium on SPEAR ######
## load taxa/spear data
inv_met <- read.csv(file.path(datadir, "invertebrate_data/macroinv_metrics.csv"), stringsAsFactors = FALSE)
## remove April data
inv_met <- inv_met[inv_met$month == "June",]
spear <- inv_met[,c(1,12)]
#spear$spear_pest <- spear$spear_pest/100

## load environmental data
env_jn <- read.csv(file.path(datadir, "/environmental_data/env_predictors_june.csv"), stringsAsFactors = FALSE)


## merge spear values to env_j
mod_dat <- merge(env_jn[,c(1,10)], spear, by = "site")
mod_dat <- mod_dat[,-1]

library(ggpubr)

summary(lm(spear_pest ~ FUS, data = mod_dat))

gr1 <- mod_dat[mod_dat$FUS == 0, "spear_pest"]
gr2 <- mod_dat[mod_dat$FUS == 1, "spear_pest"]
t.test(gr1,gr2)


## boxplot of differences in mean spear values
spear_plot <- ggplot(data=mod_dat, aes(x=as.factor(FUS), y=spear_pest, 
                                       fill = as.factor(FUS))) +
  # geom_violin(draw_quantiles = c(0.5)) + 
  stat_boxplot(geom = "errorbar", width = 0.1) +
  geom_boxplot(width = 0.2, fill = "lightgrey") +
  geom_point(size = 2) + 
  ylab(expression(SPEAR[pesticides]*" (%)")) + 
  theme_minimal() +
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank()) +
  # scale_fill_manual(values=c("lightgrey", "gray98")) +
  theme(legend.position = "none",
        axis.text.x = element_text(face="bold"),
        panel.grid.major.x = element_blank()) +
  xlab("") +
  scale_x_discrete(labels= c("without\nFUS", "with\nFUS"))

spear_plot

ggsave(file= "./plots/SPEAR.png", plot= spear_plot, width = 3.5, height = 4.5)
ggsave(file= "./plots/SPEAR.svg", plot= spear_plot, width = 3.5, height = 4.5)
ggsave(file= "./plots/SPEAR.tiff", plot= spear_plot, width = 3.5, height = 4.5)


###### Section 5 RDA Plots for April and June data ######

## See RDA Scripts 

###### Section 6 Table with environmental data ######
library(tidyr)
jn <- read.csv(file.path(datadir, "environmental_data/env_predictors_june.csv"), stringsAsFactors = FALSE)

## multiply fine_sed, pastures, arable_land with 100 to get %
jn$fineSed <- jn$fineSed * 100
jn$pastures <- jn$pastures * 100
jn$arableLand <- jn$ arableLand * 100
jn <- jn[,c(1:9,11,10)]


##june
jn.sum <- jn %>%
  summarise(across(
    .cols = c(elCond:ls_hetero), 
    .fns = list(Min = min, Max =max, Mean = mean, SD = sd), na.rm = TRUE, 
    .names = "{col}&{fn}"
  ))

jn.stats.tidy <- jn.sum %>% gather(stat, val) %>%
  separate(stat, into = c("var", "stat"), sep = "&") %>%
  spread(stat, val) %>%
  select(var, Min, Max, Mean, SD) # reorder columns
names(jn.stats.tidy) <- c("Variable", "Minimum", "Maximum", "Mean", "SD")
jn.stats.tidy %>% mutate(across(where(is.numeric), round, 3))
jn.stats.tidy$Variable = c("% arable land", "elec. Conductivity", "% fine sediment", "flow velocity m/s",
                           "ls shape index", "max sumTU iv", "O2 mg/L", "ortho-Phosphate mg/L", "% pastures")
jn.stats.tidy <- jn.stats.tidy %>% mutate(across(where(is.numeric), round, 2))
## change order of rows
jn.stats.tidy <- jn.stats.tidy[c(1,9,5,3,6,2,7,8,4),]

jn.stats.tidy$unit <- c("%", "%", "NA", "%", "NA", "mS/cm", "mg/L", "mg/L","m/s")
jn.stats.tidy <- jn.stats.tidy[,c(1,6,2:5)]

ufs_df <- data.frame(Variable = "FUS",
                     unit = "NA",
                     Minimum = "no",
                     Maximum = "yes",
                     Mean = "NA",
                     SD = "NA",
                     stringsAsFactors = FALSE)

## rbind jn.stats.tidy and ufs_df
jn_fin <- rbind(jn.stats.tidy, ufs_df)
## write to csv
write.csv(jn_fin, file.path(datadir, "/environmental_data/env_dat_jn.csv"), row.names = FALSE)



###### Section 7 Table for invertebrate metrics and LLD ######
## invertebrate metrics
inv_metrics <- read.csv(file.path(datadir, "invertebrate_data/macroinv_metrics.csv"), stringsAsFactors = FALSE)
inv_metrics <- inv_metrics[inv_metrics$month %like% "June",]
inv_metrics <- inv_metrics[,c(1,5,8,9,12)]
inv_metrics$epttax_percent <- inv_metrics$epttax_percent
inv_metrics$ept_percent <- inv_metrics$ept_percent

## shredder metrics
shr_com <- read.csv(file.path(datadir, "invertebrate_data/shredder_community.csv"), stringsAsFactors = FALSE)
shr_com <- shr_com[,c(1,2,3,5)]

## leaf litter decomposition rate k
k_data <- read.csv(file.path(datadir, "leaf_breakdown/lld_rate.csv"), stringsAsFactors = FALSE)

## merge all tables by site and melt to long format
metrics_df <- merge(inv_metrics, shr_com, by = "site")
metrics_df <- merge(metrics_df, k_data, by = "site")

metrics_df <- melt(metrics_df, value.name = "value")
metrics_summary <- metrics_df %>% 
  group_by(variable) %>%
  summarize(min = round(min(value),2), mean = round(mean(value, na.rm = TRUE),2), 
            max = round(max(value),2), sd = round(sd(value, na.rm = TRUE),2))

write.csv(metrics_summary, file.path(datadir, "invertebrate_data/SI_metrics.csv"), row.names = FALSE)

###### Section 8 Plot for dominant shredder taxa ######

## load table with shredder abundances
shr_com <- read.csv(file.path(datadir, "invertebrate_data/shredder_community.csv"), stringsAsFactors = FALSE)
shr_com <- shr_com[,-c(6:8)]

plot(shr_com$gam_shr ~ shr_com$all_shr)
cor(shr_com$gam_shr, shr_com$all_shr)
mod <- lm(gam_shr~all_shr, data = shr_com)
summary(mod)

gamma_dominance <- ggplot(shr_com, aes(x = log10(all_shr + 1), y = log10(gam_shr + 1))) +
  geom_point(shape = 21, fill = NA) +
  stat_smooth(method = "lm",
              col = "#C42126",
              se = FALSE,
              size = 1) +
  scale_x_continuous(limits = c(0,3.5),
                     breaks = c(0,1,2,3),
                     labels= c("0", "10", "100", "1000")) +
  scale_y_continuous(limits = c(0,3.5),
                     breaks = c(0,1,2,3),
                     labels= c("0", "10", "100", "1000")) +
  # theme_bw() +
  xlab("Total shredder abundance") +
  ylab(expression(paste(italic("G. balcanicus"), " abundance"))) +
  theme_minimal() +
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank()) +
  scale_fill_manual(values=c("#A6CEE3", "#1F78B4")) 
  # theme(axis.text.x = element_text(face="bold"),
  #       axis.text.y = element_text(face="bold"))
  
gamma_dominance

ggsave(file= "./plots/Gamma_shred.svg", plot=gamma_dominance, width = 4, height = 3 )
ggsave(file= "./plots/Gamma_shred.png", plot=gamma_dominance, width = 4, height = 3 )
ggsave(file= "./plots/Gamma_shred.tiff", plot=gamma_dominance, width = 4, height = 3 )

###### Section 14 Table for the environmental variables ######
## SERVER ACCESS REQUIRED!
library(RPostgreSQL)

## enable connection
source("/home/moritz/Nextcloud/PHD/Projects/Romania/Daten/Scripts/00_Connection_Data.R")
drv = dbDriver("PostgreSQL")
con = dbConnect(drv, dbname = NEW_DB, user = DBuser_m, host = DBhost, port = DBport, password = DBpassword_m)

## load tables with environmental data
reg_phch <- dbGetQuery(con, "SELECT * FROM phch.reg_phch")
sdb_phch <- dbGetQuery(con, "SELECT * FROM phch.sdb_phch")
sr_phch <- dbGetQuery(con, "SELECT * FROM phch.sr_phch")

## extract the frequently measured parameters (O2, pH, water_temp, conductivity, flow velocity)
vars <- c("o2_perc", "water_temp", "ph", "el_cond", "flow_1", "flow_2", "flow_3", "flow1_sdb", "flow2_sdb", "flow_sr")

## extract from reg_phch
reg_sub <- reg_phch[reg_phch$variable %in% vars, c(1,2,4,6)]
## remove sites x, y and z
reg_sub <- reg_sub[!reg_sub$site %in% c("X", "Y", "Z"),]
## remove measurements from April
reg_sub <- reg_sub[reg_sub$date >= "2016-05-01",]
## combine flow measures
reg_flow <- reg_sub[reg_sub$variable %like% "flow",]
reg_flow$value <- as.numeric(reg_flow$value)
reg_avgflow <- reg_flow %>%
  group_by(site, date) %>%
  summarize(value = round(mean(value, na.rm = TRUE),2))
reg_avgflow$variable <- c(rep("flow", nrow(reg_avgflow)))
## rbind values to reg_phch, remove the replicate flow measures
reg_sub <- reg_sub[!reg_sub$variable %in% c("flow_1", "flow_2", "flow_3"),]
reg_sub <- rbind(reg_sub, reg_avgflow)
reg_sub$source <- c(rep("reg", nrow(reg_sub)))
reg_sub$value <- as.numeric(reg_sub$value)

## extract from sdb_phch
sdb_sub <- sdb_phch[sdb_phch$variable %in% vars, c(1,2,6,7)]
sdb_flow <- sdb_sub[sdb_sub$variable %like% "flow",]
sdb_flow$value <- as.numeric(sdb_flow$value)
sdb_avgflow <- sdb_flow %>%
  group_by(site, date) %>%
  summarize(value = round(mean(value, na.rm = TRUE),2))
sdb_avgflow$variable <- c(rep("flow", nrow(sdb_avgflow)))
## rbind values to sdb_phch, remove the replicate flow measures
sdb_sub <- sdb_sub[!sdb_sub$variable %in% c("flow1_sdb", "flow2_sdb"),]
sdb_sub <- rbind(sdb_sub, sdb_avgflow)
sdb_sub$source <- c(rep("sdb", nrow(sdb_sub)))
sdb_sub$value <- as.numeric(sdb_sub$value)
sdb_sub[is.na(sdb_sub$value),]

## extract from sr
sr_sub <- sr_phch[sr_phch$variable %in% vars, c(1,2,5,6)]
sr_flow <- sr_sub[sr_sub$variable %like% "flow",]
sr_flow$value <- as.numeric(sr_flow$value)
sr_avgflow <- sr_flow %>%
  group_by(site, date) %>%
  summarize(value = round(mean(value, na.rm = TRUE),2))
sr_avgflow$variable <- c(rep("flow", nrow(sr_avgflow)))
## rbind values to sr_sub, remove replicate flow measures
sr_sub <- sr_sub[!sr_sub$variable %like% "flow_sr",]
sr_sub <- rbind(sr_sub, sr_avgflow)
sr_sub$source <- c(rep("sr", nrow(sr_sub)))
sr_sub$value <- as.numeric(sr_sub$value)


## rbind all three tables
phys_comb <- rbind(reg_sub, sdb_sub, sr_sub)
## there are duplicates in the data frame, the second measurement of SR is always a duplicate
## of the third measurment of sdb
## find and remove
phys_clean <- phys_comb %>%
  distinct(site, date, variable, .keep_all = TRUE)
## check for plausible values
range(phys_clean[phys_clean$variable == "water_temp",4], na.rm = TRUE)
range(phys_clean[phys_clean$variable == "ph",4], na.rm = TRUE) ## 816 is likely 8.16, check
#phys_clean[phys_clean$variable == "ph" & phys_clean$value > 9,]
range(phys_clean[phys_clean$variable == "el_cond",4], na.rm = TRUE)
range(phys_clean[phys_clean$variable == "o2_perc",4], na.rm = TRUE)
range(phys_clean[phys_clean$variable == "flow",4], na.rm = TRUE)


## count the measurements per site
observations <- phys_clean %>%
  group_by(site) %>%
  count(variable)
## bring into wide format
obs_wide <- dcast(observations, variable ~ site, value.var = "n")
names(obs_wide)[2] <- "n"

## get min, max mean and sd values 
sum_df <- phys_clean %>%
  group_by(site, variable) %>%
  summarize(min = min(value, na.rm = TRUE),
            max = max(value, na.rm = TRUE),
            mean = round(mean(value, na.rm = TRUE), 2),
            sd = round(sd(value, na.rm = TRUE), 2))
sum_df
## change into long format
sum_lf <- melt(sum_df, id = c("site", "variable"))
names(sum_lf)[3] <- "measure"
## cast into wide format with sites as columns
sum_wide <- dcast(sum_lf, variable + measure ~ site, value.var = "value")
## add row that indicates number of observations
sum_wide <- merge(obs_wide[,1:2], sum_wide, by = "variable", all.y = TRUE)
sum_wide$unit <- c(rep("mS/cm", 4), rep("m/s", 4), rep("%",4), rep(NA, 4), rep("°C", 4))
write.csv(sum_wide, file.path(datadir, "/env_data/summary_stats_wide.csv"), row.names = FALSE)

## ///////////////////// Add nutrient variables to table /////////////
vars2 <- c("no3", "no2", "nh4", "po4", "so4", "cl")
reg_nut <- reg_phch[reg_phch$variable %in% vars2, ]
## drop sites X, Y, Z and samples from April
reg_nut <- reg_nut[!reg_nut$site %in% c("X", "Y", "Z"),]
reg_nut <- reg_nut[reg_nut$date > "2016-05-01",]


reg_nut$value <- ifelse(!is.na(reg_nut$qualifier),
                        paste0(reg_nut$qualifier, " ", reg_nut$value),
                        paste0(reg_nut$value))

## change into wide format
nut_wide <- dcast(reg_nut, variable ~ site, value.var = "value")
nut_wide$unit <- c(rep("mg/L", nrow(nut_wide)))
nut_wide$n <- c(rep(1,nrow(nut_wide)))
nut_wide$measure <- c(rep(NA, nrow(nut_wide)))
        
## combine both tables
env_comb <- rbind(sum_wide, nut_wide)


## ///////////////////// Add habitat variables to table /////////////

## load table with habitat data
habi_dat <- read.csv("/home/moritz/Nextcloud/PHD/Projects/Romania/Daten/Offline_Backups/habitat.csv", stringsAsFactors = FALSE)
head(habi_dat)
## check how many varibles
unique(habi_dat$variable)
## drop the reference sites
habi_dat <- habi_dat[!habi_dat$site %in% c("X", "Y", "Z"),]
## rename O_UP and drop O_DOWN
habi_dat <- habi_dat[habi_dat$site != "O_down",]
habi_dat[habi_dat$site == "O_up", 1] <- c(rep("O",17))
habi_dat$month <- ifelse(habi_dat$date %like% "-04-",
                         paste0("April"),
                         ifelse(habi_dat$date %like% "-06-",
                                paste0("June"),
                                NA))
## drop April data
habi_dat <- habi_dat[habi_dat$month == "June",]


## Following variables for a table
## --> buffer width: min of dist_field
## --> shading
## --> main vegetation: max of bank_covers
## --> height of vegetation: max of bank_height
## --> width
## --> depth
## --> velocity

## buffer width
buffer <- habi_dat[habi_dat$variable %like% "dist_field",]
buffer$value <- as.numeric(buffer$value)
## keep the minimum width
buffer_min <- buffer %>% 
  group_by(site) %>%
  summarize(min = min(value))
names(buffer_min) <- c("site", "buffer")


## shading
shading <- habi_dat[habi_dat$variable %like% "shading", c(1,7)]
shading$value <- as.numeric(shading$value)
names(shading)[2] <- "shading"
## merge with buffer
habitat <- merge(buffer_min, shading, by = "site", all = TRUE)

## main vegetation
veg <- habi_dat[habi_dat$variable %like% "_cover_",]
veg$value <- as.numeric(veg$value)
## get maximum per site, month and direction
main_veg <- veg %>% 
  group_by(site) %>%
  summarize(maxcover = max(value), maxcovname = variable[which.max(value)])
main_veg$maxcovname <- gsub("^(\\w*?)(_)(\\w*?)(_)(\\w*?)$", "\\5", main_veg$maxcovname)
## merge with habitat
habitat <- merge(habitat, main_veg, by = "site", all = TRUE)


## average height of bank vegetation
veg_height <- habi_dat[habi_dat$variable %like% "bank_height_",]
veg_height$value <- as.numeric(veg_height$value)
## drop all entries with 0.0
veg_height <- veg_height[veg_height$value > 0,]
avg_height <- veg_height %>%
  group_by(site) %>%
  summarize(maxheight = max(value), maxheightname = variable[which.max(value)])
avg_height$maxheightname <- gsub("^(\\w*?)(_)(\\w*?)(_)(\\w*?)$", "\\5", avg_height$maxheightname)
## merge with habitat
habitat <- merge(habitat, avg_height, by = "site", all = TRUE)


## stream width
width <- habi_dat[habi_dat$variable %in% c("min_width", "max_width"), c(1,4,7)]
width$value <- as.numeric(width$value)
min_width <- width[width$variable == "min_width",]
names(min_width)[3] <- "min_width"
max_width <- width[width$variable == "max_width",]
names(max_width)[3] <- "max_width"
width <- merge(min_width, max_width, by = "site", all = TRUE)
## merge with habitat
habitat <- merge(habitat, width[,c(1,3,5)], by = "site", all = TRUE)


## stream depth
depth <- habi_dat[habi_dat$variable %like% "depth", c(1,4,7)]
depth$value <- as.numeric(depth$value)
avg_depth <- depth %>%
  group_by(site) %>%
  summarize(depth = round(mean(value, na.rm = TRUE),2))
## merge with habitat
habitat <- merge(habitat, avg_depth, by = "site", all = TRUE)


## curvature
curv <- habi_dat[habi_dat$variable %like% "curvature", c(1,4,7)]
curv$value <- as.numeric(curv$value)
curv$curvature <- ifelse(curv$value == 2,
                           paste0("stretched"),
                           ifelse(curv$value == 3,
                                  paste0("slightly curved"),
                                  ifelse(curv$value == 4,
                                         paste0("moderatly curved"),
                                         ifelse(curv$value == 5,
                                                paste0("strongly curved"),
                                                ifelse(curv$value == 6,
                                                       paste0("winding"),
                                                       ifelse(curv$value == 7,
                                                              paste0("meandering"),
                                                              NA))))))

## merge with habitat
habitat <- merge(habitat, curv[,c(1, 4)], by = "site", all = TRUE)


## riffles
riff <- habi_dat[habi_dat$variable %like% "riffels", c(1,4,7)]
riff$value <- as.numeric(riff$value)
names(riff)[3] <- "riffles"
## merge with habitat
habitat <- merge(habitat, riff[,c(1,3)], by = "site", all = TRUE)


## pools
pool <- habi_dat[habi_dat$variable %like% "pools", c(1,4,7)]
pool$value <- as.numeric(pool$value)
names(pool)[3] <- "pools"
## merge with habitat
habitat <- merge(habitat, pool[,c(1,3)], by = "site", all = TRUE)

## substrate
substrate <- dbGetQuery(con, "SELECT * FROM habitat.sites_substrate")
## drop sites X, Y, Z
substrate <- substrate[!substrate$site %in% c( "X", "Y", "Z"),]
## rename site O_up to O and drop O_down
substrate[substrate$site == "O_up", ]$site <- c(rep("O", 19))
substrate <- substrate[substrate$site != "O_down", ]
## second sampling here
substrate <- substrate[substrate$date >= "2016-06-01",]
substrate <- substrate[substrate$variable %in% c("Psammal", "Argyllal"),]
substrate$value <- as.numeric(substrate$value)
fine_substrate <- substrate %>%
  group_by(site) %>%
  summarize(fine_sed = sum(value))

## merge with habitat
habitat <- merge(habitat, fine_substrate, by = "site", all = TRUE)
## melt table
habi_melt <- melt(habitat, id = c("site"), value.name = "value")
# cast into dataframe with sites as column names
habi_wide <- dcast(habi_melt, variable ~ site, value.var = "value")
## add columns with units, number of measurements etc
habi_wide$unit <- c("m", "%", "%", NA, "m", NA, "m", "m", "cm", "NA", "%", "%", "%")
habi_wide$n <- c(rep("1", nrow(habi_wide)))
habi_wide$measure <- c(rep(NA, nrow(habi_wide)))
## change order of columns

## combine with env_comb
table_fin <- rbind(env_comb, habi_wide)
table_fin <- table_fin[,c(1,23,3,2,4:22)]
## save table
write.csv(table_fin, file.path(datadir, "/env_data/SI_maintable.csv"), row.names = FALSE)


###### Section 15 Plot for the SPEAR values in April and June ######
## load macroinvertebrate metric
ma_metric <- read.csv("./data/invertebrate_data/macroinv_metrics.csv", stringsAsFactors = FALSE)
ma_metric$spear_pest <- ma_metric$spear_pest

boxplot(ma_metric[ma_metric$month == "April",]$spear_pest, ma_metric[ma_metric$month == "June",]$spear_pest)


## get the envrionmental data, especially FUS
env <- read.csv("./data/environmental_data/env_predictors_june.csv", stringsAsFactors = FALSE)

## merge both dfs
df <- merge(ma_metric[,c(1,2,12)], env[,c(1,10)], by = "site")


## boxplots
boxplot(df[df$month == "April" & df$FUS == 0,]$spear_pest,df[df$month == "April" & df$FUS == 1,]$spear_pest,
        df[df$month == "June" & df$FUS == 0,]$spear_pest,df[df$month == "June" & df$FUS == 1,]$spear_pest,
        names = c("no FUS April", "FUS April", "no FUS June", "FUS June"), main = "SPEAR values")

mean(df[df$month == "April" & df$FUS == 0,]$spear_pest)
mean(df[df$month == "April" & df$FUS == 1,]$spear_pest)

## same with ggplot!

## boxplot of differences in mean spear values
df$dummy <- paste0(df$month, df$FUS)
df <- df[order(df$dummy),]
apriljune_plot <- ggplot(data=df, aes(x=as.factor(dummy), y=spear_pest, 
                                       fill = as.factor(dummy))) +
  #geom_violin(draw_quantiles = c(0.5)) + 
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(width = 0.3, fill = "lightgrey") +
  geom_point(size = 2) + 
  ylab(expression(SPEAR[pesticides])) + 
  theme_minimal() +
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank()) +
  # scale_fill_manual(values=c("lightgrey", "gray98", "lightgrey", "gray98")) +
  theme(legend.position = "none",
        axis.text.x = element_text(face="bold"),
        panel.grid.major.x = element_blank()) +
  xlab("") +
  scale_x_discrete(labels= c("April\nwithout FUS", "April\nwith FUS", "June\nwithout FUS", "June\nwith FUS"))

apriljune_plot


ggsave(file= "./plots/ApJnSPEAR.png", plot= apriljune_plot, width = 3.5, height = 4.5)
ggsave(file= "./plots/ApJnSPEAR.svg", plot= apriljune_plot, width = 3.5, height = 4.5)


###### Section 16 Barplot for Parasite Prevalence ######
## load csv file
prev <- read.csv("./data/parasite_data/para_prevalence.csv", stringsAsFactors = FALSE)
## prepare data for ggplot
prev_clean <- prev[-17,]
## melt into long format
prev_long <- melt(prev_clean, id.vars = "Site", value.name = "prevalence")
## change value column to numeric
prev_long$prevalence <- as.numeric(prev_long$prevalence)
## split into two df
prev_gam <- prev_long[prev_long$variable %like% "_Gam",]
prev_bae <- prev_long[prev_long$variable %like% "_Bae",]


## set colors 
colors_gam <- c("Microsporidia_Gam" = "lightgrey", "Nematoda_Gam" = "grey60", "Trematoda_Gam" = "slategray4")

## plot for Gammarids
ggplot(data = prev_gam, aes(x = Site, y = prevalence, fill = variable)) +
  geom_bar(stat = "identity", width = 1.7, position = position_dodge(width = -0.5)) +
  scale_fill_manual(values = colors_gam, name = "",
                    labels = c("Microsporidia", "Nematoda", "Trematoda")) +
  theme_minimal() +
  theme(#legend.position = "none",
        #axis.text.x = element_text(face="bold"),
        panel.grid.major.x = element_blank()) +
  xlab("") + ylab("Parasite prevalence in %")



## set colors
colors_bae <- c("Microsporidia_Bae" = "lightgrey", "Nematoda_Bae" = "grey60", "Trematoda_Bae" = "slategray4")

## plot for Baetis
ggplot(data = prev_bae, aes(x = Site, y = prevalence, fill = variable)) +
  geom_bar(stat = "identity", width = 1.7, position = position_dodge(width = -0.5)) +
  scale_fill_manual(values = colors_bae, name = "",
                    labels = c("Microsporidia", "Nematoda", "Trematoda")) +
  theme_minimal() +
  theme(#legend.position = "none",
    #axis.text.x = element_text(face="bold"),
    panel.grid.major.x = element_blank()) +
  xlab("") + ylab("Parasite prevalence in %")



## plot with prev_long and facet grid
## add dummy id to prev_long
prev_long$dummy <- factor(gsub(".*_", "",prev_long$variable), levels = c("Gam", "Bae"))
prev_long$variable <- gsub("_.*", "",prev_long$variable)

  

## set colors
colors <- c("Microsporidia" = "lightgrey", "Nematoda" = "grey60", "Trematoda" = "slategray4")
## set labels for grids
prev_long$d2 <- factor(prev_long$dummy, labels = c(expression(paste(italic("G. balcanicus"))), 
                                                              expression(paste(italic("Baetis")," sp."))))
## set text labels
text_df <- data.frame(Site = c("B", "H", "I", "S", "G", "P"), prevalence = 35, lab=c(rep("No host individuals", 6)),
                      dummy = factor(c(rep("Gam",4), rep("Bae", 2)), levels = c("Gam", "Bae")), 
                      variable = c(rep("Nematoda",6)), 
                      d2 = factor(c(rep("Gam", 4), rep("Bae", 2)), 
                                  labels = c(expression(paste(italic("G. balcanicus"))), 
                                             expression(paste(italic("Baetis")," sp."))),
                                  levels = c("Gam", "Bae")),
                      stringsAsFactors = FALSE)


## plot for both
paras <- ggplot(data = prev_long, aes(x = Site, y = prevalence, fill = variable)) +
  geom_bar(stat = "identity", width = 1.7, position = position_dodge(width = -0.5)) +
  scale_fill_manual(values = colors, name = "",
                    labels = c("Microsporidia", "Nematoda", "Trematoda")) +
  theme_minimal() +
  xlab("") + ylab("Parasite prevalence in %") + 
  facet_wrap(~d2, labeller = label_parsed, ncol = 1) +
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        strip.text = element_text(face = "bold", size = 14)) +
  geom_text(data = text_df, mapping = aes(x = Site, y = prevalence, label = lab),
            angle = 90)
    
paras

# ## create plot for the mean values per site
# mean_prev <- data.frame(prev_long %>%
#   group_by(dummy, Site) %>%
#   summarize(mean_prev = round(mean(prevalence, na.rm = TRUE),2)),
#   stringsAsFactors = FALSE)
# 
# mean_prev$d2 <- factor(mean_prev$dummy, labels = c(expression(paste(italic("G. balcanicus"))), 
#                                                    expression(paste(italic("Baetis")," sp."))))
# 
# 
# mean_plot <- ggplot(data = mean_prev, aes(x = Site, y = mean_prev)) +
#   geom_point(shape = 5) + theme_minimal() +
#   xlab("") + ylab("") + ylim(c(0,100)) + 
#   facet_wrap(~d2, ncol = 1, labeller = label_parsed) +
#   theme(panel.grid.major.x = element_blank(),
#         panel.grid.major.y = element_blank(),
#         panel.grid.minor.y = element_blank(),
#         axis.text = element_text(size = 12),
#         axis.title.y = element_text(size = 12),
#         legend.text = element_text(size = 12),
#         strip.text = element_text(face = "bold", size = 14))
# 
# mean_plot
# 
# library(cowplot)
# 
# aligend_plots <- align_plots(paras, mean_plot, align = "hv", axis = "tblr")
# paras_fin <- ggdraw(aligend_plots[[1]]) + draw_plot(aligend_plots[[2]])
# 
# ggsave(file= "./model_plots/parasite_prevalence.svg", plot=paras_fin, width = 8.3, height = 5.85)
# ggsave(file= "./model_plots/parasite_prevalence.png", plot=paras_fin, width = 8.3, height = 5.85)
# ggsave(file= "./model_plots/parasite_prevalence.tiff", plot=paras_fin, width = 8.3, height = 5.85)

## plot with mean values removed
ggsave(file= "./plots/parasite_prevalence.svg", plot=paras, width = 8.3, height = 5.85)
ggsave(file= "./plots/parasite_prevalence.png", plot=paras, width = 8.3, height = 5.85)
ggsave(file= "./plots/parasite_prevalence.tiff", plot=paras, width = 8.3, height = 5.85)


###### Section 17 checking effects of revised SPEAR classification ######
## Macroinvertebrate data in long format
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
                              paste0(mzb_dat$genus, " sp."),
                              ifelse(!is.na(mzb_dat$family),
                                     paste0(mzb_dat$family),
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

## Taxa reclassified in knillmann et al: Anabolia nervosa, Ephemerellidae, Ephemeridae, Leptophlebiidae
## Leptoceridae, Nemouridae, Sericostomatidae, Tabanidae

anabolia <- mzb_june[mzb_june$genus %like% "Anabolia",]
ephemerellidae <- mzb_june[mzb_june$family %like% "Ephemerellidae",]
ephemeridae <- mzb_june[mzb_june$family %like% "Ephemeridae",]
ironoquia <- mzb_june[mzb_june$genus %like% "Ironoquia",]
leptophlebiidae <- mzb_june[mzb_june$family %like% "Leptophlebiidae",]
leptoceridae <- mzb_june[mzb_june$family %like% "Leptoceridae",]
nemouridae <- mzb_june[mzb_june$family %like% "Nemouridae",]
sericostomatidae <- mzb_june[mzb_june$family %like% "Sericostomatidae",]
tabanidae <- mzb_june[mzb_june$family %like% "Tabanidae",]

## rbind tables 
not_spear <- rbind(anabolia, ephemerellidae, leptophlebiidae, leptoceridae, nemouridae, tabanidae)
## count individuals per site
ns_site <- not_spear %>%
  group_by(site) %>%
  summarize(recla = sum(abundance))

## load taxa/spear data
spear_dat <- read.csv(file.path(datadir, "invertebrate_data/spear_taxa.csv"), stringsAsFactors = FALSE)
spear_pos <- spear_dat[spear_dat$spear_class == 1,]
## get numbers per site
spear_numbers <- spear_pos %>%
  group_by(site) %>%
  summarize(spear = sum(abundance))

## merge with not_spear
numbers <- merge(spear_numbers, ns_site, by = "site")
sp <- sum(numbers$spear)
not_sp <- sum(numbers$recla)
606/1182


