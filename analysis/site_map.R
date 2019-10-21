# Making map of study sites - Peru
# John Godlee (johngodlee@gmail.com)
# 2017_03_09

# Packages ----
library(maps)  # Draw Geographical Maps. Display of maps.
library(mapdata)  # Projects map data
library(maptools)  # Tools for reading and handling spatial objects
library(gpclib)  # to clip polygons from shapefiles
library(raster)  # Reading, writing, manipulating, analyzing and modeling of gridded spatial data
library(rgdal)  # Projection/transformation operations for GDAL raster and OGR vector map data
library(ggmap)
library(ggsn)
library(plyr)  # join()
library(dplyr)
library(viridis)
library(ggrepel)

#EPSG <- make_EPSG()

# set working directory ----
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load data ----

# Seedlings
seedlings <- read.csv("data/seedlings.csv")

# Manu shapefile
manu_wdpa <- readOGR("data/manu_wdpa_outline", "manu_wdpa")
proj4string(manu_wdpa)

# Peru shapefile
peru <- readOGR("data/peru_outline", "peru")
proj4string(peru)

# Elevation data for Peru using getData() from {raster}
elevation_peru <- getData("alt", country = "Peru")

# Create objects for map of Peru with manu outline ----
# convert peru raster to points for plotting
elev_peru_p <- rasterToPoints(elevation_peru)
elev_peru_df <- data.frame(elev_peru_p)
colnames(elev_peru_df) <- c("lon", "lat", "elev")

# Convert manu polygon to df for plotting
manu_wdpa@data$id <- rownames(manu_wdpa@data)
manu_p = fortify(manu_wdpa, region="id")
manu_df = join(manu_p, manu_wdpa@data, by="id")

# Create objects for map of study site ----

site_loc <- seedlings %>%
  group_by(Site) %>%
  dplyr::summarise(lat_mean = mean(Lat.DD),
    lon_mean = mean(Long.DD),
    elev = mean(Elevation))

loc_extent <- c(
  min(site_loc$lon_mean) - 0.05, min(site_loc$lat_mean) - 0.05, 
  max(site_loc$lon_mean) + 0.05, max(site_loc$lat_mean) + 0.05)

plot_poly <- as(raster::extent(loc_extent[1], loc_extent[3], loc_extent[2], loc_extent[4]), "SpatialPolygons")

plot(plot_poly)

plot_poly_f <- fortify(plot_poly)

# Plot of peru with manu outline on top ----
pdf("../manuscript/img/peru_map.pdf", width = 10, height = 12)
(peru_map_ggplot <- ggplot(data = elev_peru_df, aes(x = lon, y = lat)) +
  geom_tile(aes(fill = elev)) + 
  geom_polygon(data = manu_df, aes(x = long, y = lat, group = group), 
    fill = "grey", colour = "white", alpha = 0.4) + 
  geom_polygon(data = plot_poly_f,  aes(x = long, y = lat), fill = "#FC7C7C", alpha = 0.5, colour = "red", size = 2) + 
  coord_quickmap() + 
  scale_fill_viridis() + 
  theme_classic() + 
  guides(fill = guide_legend(title="Elevation (m)")) + 
  labs(y = "Latitude", x = "Longitude")) 
dev.off()

# Get map
manu_gg <- get_map(location = loc_extent, source = "google", maptype = "terrain", color = "bw", zoom = 12, scale = 2)

pdf("../manuscript/img/plot_map.pdf", width = 10, height = 12)
ggmap(manu_gg) +
  geom_point(data = site_loc,
    aes(x = lon_mean, y = lat_mean), 
    shape = 21, colour = "black", fill = "#E06614", size = 4) + 
  geom_label_repel(data = site_loc,
    aes(x = lon_mean, y = lat_mean, label = Site),
    box.padding = 1, nudge_x = 0.03) + 
  theme_classic() + 
  labs(x = "Longitude", y = "Latitude") + 
  coord_map()
dev.off()
