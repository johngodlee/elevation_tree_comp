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
manu <- readOGR("data/manu_outline", "manu")
# plot(manu, axes = T)
proj4string(manu)

# Peru shapefile
peru <- readOGR("data/peru_outline", "peru")
# plot(peru, axes = T)
proj4string(peru)

# ggmap of manu outline
bb_manu <- bbox(manu)
manu_ggmap <- ggmap(get_map(location = bb_manu, zoom = 9))
manu_ggmap

# Map of Peru with Elevation data using getData() from {raster}
elevation <- getData("alt", country = "Peru")
# convert raster to points for plotting
elev_p <- rasterToPoints(elevation)

# Make the points a dataframe for ggplot
elev_df <- data.frame(elev_p)

# Make appropriate column headings
colnames(elev_df) <- c("lon", "lat", "elev")

manu@data$id <- rownames(manu@data)
manu_p = fortify(manu, region="id")
manu_df = join(manu_p, manu@data, by="id")


# Plot with manu outline on top
ggplot(data = elev_df, aes(x = lon, y = lat)) +
	geom_raster(aes(fill=elev)) + 
	geom_polygon(data = manu_df, aes(x = long, y = lat, group = group), 
							 fill = "grey", colour = "white", alpha = 0.5) + 
	coord_cartesian() + 
	scale_fill_viridis() + 
	theme_classic() + 
	guides(fill = guide_legend(title="Elevation (m)"))

# Find out the position of each plot from mean of seedling locations ----
site_loc <- seedlings %>%
	group_by(Site) %>%
	summarise(lat_mean = mean(Lat.DD),
						lon_mean = mean(Long.DD),
						elev = mean(Elevation))


# Maps of study sites, zoomed in on Manu



#Plot maps
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

