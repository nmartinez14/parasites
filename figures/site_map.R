## **********************************************************
## Load libraries and source files
## **********************************************************
rm(list=ls())
setwd("~/")
source("lab_paths.R")
setwd(file.path(local.path))
setwd("parasites/figures")

load("../parasitesData.Rdata")

## packages
library(sf)
library(ggplot2)
library(terra)
library(ggspatial)
library(basemaps)
library(tidyverse)

## **********************************************************
## Make map
## **********************************************************

ggplot() +
  geom_sf(data = sites_sf) +
  ggtitle("Map of Plot Locations")


## set defaults for the basemap
set_defaults(map_service = "esri", map_type = "usa_topo_maps")
## Change the CRS to match the basemaps
site_points <- st_transform(sites_sf, crs = st_crs(3857))
## Create a new bounding box to avoid points in the corners
bbox_new <- st_bbox(site_points) # current bounding box

xrange <- bbox_new$xmax - bbox_new$xmin # range of x values
yrange <- bbox_new$ymax - bbox_new$ymin # range of y values

bbox_new[1] <- bbox_new[1] - (0.1 * xrange) # xmin - left
bbox_new[3] <- bbox_new[3] + (0.1 * xrange) # xmax - right
bbox_new[2] <- bbox_new[2] - (0.1 * yrange) # ymin - bottom
bbox_new[4] <- bbox_new[4] + (0.1 * yrange) # ymax - top
bbox_new <- bbox_new %>%  # take the bounding box make it a spatial object
  st_as_sfc()

map <- ggplot() +
  basemap_gglayer(bbox_new) + # Use new bbox to download the basemap
  geom_sf(data = site_points, color = "black") +
  coord_sf(xlim = st_coordinates(bbox_new)[c(1,2),1], # min & max of x values
           ylim = st_coordinates(bbox_new)[c(2,3),2], expand = FALSE) +
  geom_sf_text(data = site_points, aes(label = MtRange), size = 2.5,
               color = "black", nudge_y = 19500) +
  scale_fill_identity()+ 
  xlab("Longitude") + ylab("Latitude") +
  annotation_north_arrow(location = "tl", 
                         style = north_arrow_fancy_orienteering)+
  annotation_scale(location = "br")+
  theme_minimal()

## Toggle on to save out figure, updating filepath to your desired save location
save.fig = FALSE
if (save.fig == TRUE){
  setwd("parasites/figures")
  ggsave(map, file="map.pdf", height=6, width=4)
}

