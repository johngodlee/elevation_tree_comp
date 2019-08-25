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
# plot(peru, axes = T)
proj4string(peru)

# Elevation data for Manu in raster format from USGS
elevation_1 <- raster("data/usgs_dem/ASTGTM2_S13W071_dem.tif")
elevation_2 <- raster("data/usgs_dem/ASTGTM2_S13W072_dem.tif")
elevation_3 <- raster("data/usgs_dem/ASTGTM2_S13W073_dem.tif")
elevation_4 <- raster("data/usgs_dem/ASTGTM2_S14W071_dem.tif")
elevation_5 <- raster("data/usgs_dem/ASTGTM2_S14W072_dem.tif")
elevation_6 <- raster("data/usgs_dem/ASTGTM2_S14W073_dem.tif")
elevation_all <- merge(elevation_1, elevation_2, elevation_3, elevation_4, elevation_5, elevation_6)

# Elevation data for Peru using getData() from {raster}
elevation_peru <- getData("alt", country = "Peru")

# Peru rivers lines
water <- readOGR("data/peru_water", "sa_riv_15s")

# Create objects for map of Peru with manu outline ----
# convert peru raster to points for plotting
elev_peru_p <- rasterToPoints(elevation_peru)
elev_peru_df <- data.frame(elev_peru_p)
colnames(elev_peru_df) <- c("lon", "lat", "elev")

# Convert manu polygon to df for plotting
manu_wdpa@data$id <- rownames(manu_wdpa@data)
manu_p = fortify(manu_wdpa, region="id")
manu_df = join(manu_p, manu_wdpa@data, by="id")

# Plot of peru with manu outline on top ----

peru_map_ggplot <- ggplot(data = elev_peru_df, aes(x = lon, y = lat)) +
  geom_tile(aes(fill = elev)) + 
  geom_polygon(data = manu_df, aes(x = long, y = lat, group = group), 
    fill = "grey", colour = "white", alpha = 0.4) + 
  coord_quickmap() + 
  scale_fill_viridis() + 
  theme_classic() + 
  guides(fill = guide_legend(title="Elevation (m)")) + 
  labs(y = "Latitude", x = "Longitude")

ggsave(filename = "../manuscript/img/peru_map.pdf", plot = peru_map_ggplot, width = 5, height = 9)

# Create objects for map of study site ----
# Find out the position of each plot from mean of seedling locations
site_loc <- seedlings %>%
  group_by(Site) %>%
  dplyr::summarise(lat_mean = mean(Lat.DD),
    lon_mean = mean(Long.DD),
    elev = mean(Elevation))

# Create bounding box for map based on extent of plot locations
site_bbox <- as(raster::extent(-71.75, -71.25, -13.25, -12.5), "SpatialPolygons")

# Crop elevation dem to site locations
elevation_clipped <- crop(elevation_all, extent(site_bbox))
elevation_clipped <- mask(elevation_clipped, site_bbox)

# Convert cropped elevation dem to data frame for ggplot
elev_spdf <- as(elevation_clipped, "SpatialPixelsDataFrame")
elev_df <- as.data.frame(elev_spdf)
colnames(elev_df) <- c("elev", "lon", "lat")

# Crop water to site locations
# water_clipped <- crop(water, raster::extent(site_bbox))
# Write to object
# writeSpatialShape(x = water_clipped, fn = "water_clipped")

# Convert cropped water lines to data frame for ggplot
water_sldf <- as(water, "SpatialLinesDataFrame")
water_clipped <- raster::intersect(water_sldf, site_bbox)
water_df <- fortify(water_clipped)


# Crop Manu outline to site_loc bbox
manu_wdpa_clipped <- crop(manu_wdpa, extent(site_bbox))
# Convert manu polygon to df for plotting
manu_clipped_p = fortify(manu_wdpa_clipped, region="id")
manu_clipped_df = join(manu_clipped_p, manu_wdpa@data, by="id")


# ggmap of study sites, zoomed in on Manu ----

site_map_ggplot <- ggplot() +
	geom_raster(data = elev_df, aes(x = lon, y = lat, fill=elev)) +  # DEM
	geom_polygon(data = manu_clipped_df, aes(x = long, y = lat, group = group), 
							 fill = "grey", colour = "white", alpha = 0.3) + 
	geom_point(data = site_loc, aes(x = lon_mean, y = lat_mean), 
						 size = 5, alpha = 0.8, fill = "#D64E0F", colour = "black", pch=21) +  # Site locations
	#geom_path(data = water_df, aes(x = long, y = lat, group=group), colour = "blue") +  # River locations
	scale_fill_viridis() +  # Colour for dem
	theme_classic() +
	xlab(expression(paste("Longitude (",degree,")"))) + 
	ylab(expression(paste("Latitude (",degree,")")))
ggsave(file = "map_test.pdf", plot = site_map_ggplot)

# Would be good to have river on the map
# Map is all off kilter now, fuck it

# Save all plots ----

# Create list of plots using grep
plot_pattern<-grep("_ggplot", names(.GlobalEnv), value=TRUE)
ggplot_list<-do.call("list", mget(plot_pattern))

# Compile
mapply(ggsave, 
			 file = paste0("img/", names(ggplot_list), ".pdf"), 
			 plot = ggplot_list, 
			 width = 30, height = 25, units = "cm")

############################################
######### OLD MAPS #########################
############################################

# Plot maps
myLocation <- c(min(camp_loc$Site_long_mean.Long.DD_mean), min(camp_loc$Lat.DD_mean), max(camp_loc$Site_long_mean.Long.DD_mean), max(camp_loc$Lat.DD_mean))
Manu_gg <- get_map(location = myLocation, zoom = 10, maptype = "terrain", )
ggmap(Manu_gg) + scale_shape_identity() + 
  geom_point(aes(x = lon, y = lat), data = camp_loc, color = "red2", fill = "red2", shape = 20, size = 5) +
  theme(axis.title.x = element_text(size = 20)) +
  theme(axis.title.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.text.y = element_text(size = 20)) +
  scale_x_continuous(limits = c(-71.65, -71.2)) + 
  scale_y_continuous(limits = c(-13.2, -12.6)) +
  xlab(expression(paste("Longitude (",degree,")"))) + 
  ylab(expression(paste("Latitude (",degree,")"))) +
  scaleBar(lon = -71.3, lat = -13.2, distanceLon = 10, distanceLat = 1.5, distanceLegend = 10000, dist.unit = "km")  #10km scale bar

Trocha_camp <- filter(camp_loc, Site %in% c("TRU02", "TRU04", "TRU06", "TRU07", "TRU08", "SP1500", "SP1750"))
Trocha_location <- c(min(Trocha_camp$lon), min(Trocha_camp$lat), max(Trocha_camp$lon), max(Trocha_camp$lat))
Trocha_gg <- get_map(location = Trocha_location, zoom = 13, maptype = "terrain")


  ggmap(Trocha_gg) + scale_shape_identity() + geom_point(aes(x = lon, y = lat), data = Trocha_camp, color = "red", shape = 43, size = 15) +
    theme(axis.title.x = element_text(size = 25)) +
    theme(axis.title.y = element_text(size = 25)) +
    theme(axis.text.x = element_text(size = 25)) +
    theme(axis.text.y = element_text(size = 25)) +
    xlab(expression(paste("Longitude (",degree,")"))) + 
    ylab(expression(paste("Latitude (",degree,")"))) +  
    scaleBar(lon = -71.540, lat = -13.125, distanceLon = 1, distanceLat = 0.1, distanceLegend = 0, dist.unit = "km")


dev.off()
#For Presentation

peru_google <- get_map(location = "Peru", maptype = "terrain", source = "google", zoom = 4)

data(world.cities)
Peru_outline <- data.frame(map("world", "Peru", plot=FALSE)[c("x","y")])
    ggplot(Peru, aes(x=x, y=y)) +
      geom_path(colour = 'green') +
      coord_map() + theme_bw()

peru_google_ggmap <- ggmap(peru_google)
peru_google_ggmap + geom_polygon(data = Peru_outline, aes(x=x, y=y), fill = "black", alpha = .75, size = .2)

PER_rds <- readRDS("/Users/johngodlee/Downloads/PER_adm0.rds")
PER_rds
plot(PER_rds, col = 'lightgrey', border = 'darkgrey')

peru_google_ggmap + 
  geom_point(data = PER_rds, aes(x=long, y=lat), colour = "red2", alpha = 1, size = 0.8) + 
  xlab(expression(paste("Longitude ("*degree*")"))) + ylab(expression(paste("Latitude ("*degree*")"))) +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20))


