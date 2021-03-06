# Making a map of plot locations
# John Godlee (johngodlee@gmail.com)
# 2019_10_21

# Packages
library(ggmap)
library(dplyr)

# setwd ----
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Remade maps for new pres
seedlings <- read.csv("R_Seedlings.csv")
plot_locations <- seedlings %>%
  group_by(., Site) %>%
  summarise(lat = mean(Lat.DD), lon = mean(Long.DD), elev = mean(Elevation)) %>%
  mutate(forest_type_intuitive = c("lowland", "lowland", "transition", "transition", "montane", "montane", "montane", "montane", "montane", "lowland"),
         forest_type_malhi = c("lowland", "lowland", "lowland", "montane", "montane", "montane", "montane", "montane", "montane", "lowland"))


bbox = c(min(plot_locations$lon) - .5,
         min(plot_locations$lat) - .5,
         max(plot_locations$lon) + .5,
         max(plot_locations$lat) + .5)

plot_map <- get_map(location = bbox, maptype = "terrain", source = "google", color = "bw")

ggmap(plot_map) + 
  geom_point(data = plot_locations, 
             aes(x = lon,
                 y = lat, 
                 colour = forest_type_malhi), 
             size = 8, shape = 18, alpha = 0.9) + 
  scale_colour_brewer(palette = "Dark2") + 
  theme(legend.position = "none") + 
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 15)) + 
  xlab(expression(paste("Longitude (",degree,")"))) + 
  ylab(expression(paste("Latitude (",degree,")"))) + 
  scale_x_continuous(limits = c(-71.65, -71.15)) + 
  scale_y_continuous(limits = c(-13.2, -12.6))

## Country map
country_map <- get_map(location = "Peru", maptype = "terrain", source = "google", color = "bw", zoom = 6)
ggmap(country_map) +
  geom_point(aes(x = -71.54309, y = -13.04688), size = 10, shape = 19, color = "#CF8700") + 
  geom_point(aes(x = -71.54309, y = -13.04688), size = 13, shape = 1, color = "#CF8700") + 
  geom_point(aes(x = -71.54309, y = -13.04688), size = 15, shape = 1, color = "#CF8700") + 
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 15)) + 
  xlab(expression(paste("Longitude (",degree,")"))) + 
  ylab(expression(paste("Latitude (",degree,")"))) 
