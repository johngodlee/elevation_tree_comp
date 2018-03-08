# Dissertation manuscript
# John Godlee (johngodlee@gmail.com)
# 2017_10_09

# Packages ----
library(ggplot2)
library(dplyr)
library(gridExtra)
library(tidyr)
library(stargazer)
library(reshape2)
library(lme4)
library(glmmADMB) 
library(usdm)

# Set working directory to the location of the source file ----
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load and fix data ----

## Seedlings
seedlings <- read.csv("data/seedlings.csv")
seedlings_comp <- filter(seedlings, Comp.Y.N. == "Y")

### Checking variables for collinearity

VIF <- vif(data.frame(seedlings_comp$LAI.4.ring, 
					 seedlings_comp$Comp.seed.total, 
					 seedlings_comp$Comp.adult.log.metric, 
					 seedlings_comp$Elevation))
VIF  # VIF is less than 4 for all, nearly 1, so no multicollinearity

### Removing all NA rows for analysis

seedlings_rem_na <- seedlings_comp %>%
	filter(!is.na(LAI.4.ring)) %>%
	filter(!is.na(Comp.seed.total)) %>%
	filter(!is.na(Comp.adult.log.metric)) %>%
	filter(!is.na(Elevation)) %>%
	filter(!is.na(Species))

## Genus Level migration rates (Feeley et al. 2011)
genus_mig_rates <- read.csv("data/genus_mig_rates.csv")

## Rank abundance data
aberg_census <- read.csv("data/aberg.csv")

## Soil data 
soil <- read.csv("data/soil.csv")

# Plotting interaction plots of plant traits over elevation ----

## Add LMA (Leaf mass / Area)
seedlings$LMA <- (seedlings$Leaf.mass.dry.g / seedlings$Leaf.area)  

## Add Height:Leaf Ratio (Height / no. leaves)
seedlings$Height.leaf.ratio <- (seedlings$Height.cm / seedlings$No.leaves)

## Variables vs Elevation.code - Interaction plots - grouped by species
seedlings_melt_traits <- melt(seedlings, id.vars=c("Individual.code", 
																			 "Elevation.code", "Species", "Elevation"), 
						 measure.vars = c("Leaf.area", "Leaf.mass.dry.g", "Stem.volume.mm3", 
						 								 "LMA", "Height.leaf.ratio", 
						 								"Leaf.thick.mean.mm", "SPAD.mean", "D.FvFm"))

## Make the plot
ggplot(seedlings_melt_traits, aes(x=Elevation.code, y=value, colour=Species, group=Species)) + stat_summary(fun.y=mean, geom="point") +
	stat_summary(fun.y=mean, geom="line") + 
	scale_x_discrete(limits=c("Top","Middle","Bottom")) + 
	facet_wrap(~variable, scales = "free")

# Plotting variation in plant traits across species ----

## Make the plot
ggplot(seedlings_melt_traits, aes(x=Species, y=value, colour=Species)) + 
	geom_point() + 
	facet_wrap(~variable, scales = "free")

# Plotting all seedlings by elevation ----
seedlings_comp$Elevation_rank <- rank(seedlings_comp$Elevation, na.last=TRUE, ties.method="random")
ggplot(seedlings_comp, aes(x=Elevation_rank, y=Elevation, colour=Species)) + geom_point()

# Methods tables and plots ----

## What is the max and min range of each species
species <- c("HG","AV","CT","SP","CR","MS","TG","DL","ID")
max <- c(3250, 2250, 2375, 2250, 2875, 2250, 1800, 1800, 1500)
min <- c(1800, 1750, 1625, 1800, 1625, 595, 1250, 1500, 425)
ranges <- data.frame(species, max, min) %>%
	arrange(., species) %>%
	gather(key = "max_min", value = "range", max:min)

ranges_spread <- ranges %>%
	spread(max_min, range)

## What are the camp locations
camp_loc <- seedlings %>%
	group_by(Site) %>%
	summarise(lat_mean = mean(Lat.DD),
						lon_mean = mean(Long.DD)) %>%
	arrange(., Site) %>%
	mutate(elev_mean = c(406, 822, 1497, 1770, 3213, 2733, 2281, 2135, 1839, 861)) %>%
	arrange(., elev_mean)

## Where will each species be sampled
unique_elevation.code_species_Site <- unique(seedlings[c("Elevation.code", "Species", "Site")])
Top_range_site <- c("TRU02", "TRU07", "TRU06","TRU06", "TRU04", "TRU07", "TRU08", "TRU08", "VC")
Middle_range_site <- c("TRU04", "TRU08", "TRU07", "TRU07", "TRU06", "SP1500", NA, NA, NA)
Bottom_range_site <- c("TRU08", "SP1750", "SP1750", "TRU08", "TRU08/SP1750", "PA800", "SP1500", "SP1500", "PA400")
Species_sample_loc <- data.frame(species, Top_range_site, Middle_range_site, Bottom_range_site)
Species_sample_loc <- arrange(Species_sample_loc, species)
stargazer(Species_sample_loc, type = "latex", summary = F)

## Table of elevations of measurement
Species_all_3 <- c("AV", "AV", "AV", "CR", "CR", "CR", "CR", "CT", "CT", "CT", "DL","DL","DL", "HG", "HG", "HG", "ID", "ID", "ID", "MS", "MS", "MS", "SP", "SP", "SP", "TG", "TG", "TG")
Position_all_3 <- c("Top", "Middle", "Bottom","Top", "Middle", "Bottom", "Bottom", "Top", "Middle", "Bottom","Top", "Middle", "Bottom","Top", "Middle", "Bottom","Top", "Middle", "Bottom","Top", "Middle", "Bottom","Top", "Middle", "Bottom","Top", "Middle", "Bottom")
Range_all_3 <- c(2135, 1839, 1770, 2733, 2281, 1839, 1770, 2281, 2135, 1770, 1839, NA, 1497, 3212, 2733, 1839, 861, NA, 406, 2135, 1497, 822, 2281, 2135, 1839, 1839, NA, 1497)
Species_all_loc <- data.frame("species" = Species_all_3, "position" = Position_all_3, "range" = Range_all_3)

## How many species are at each elevation 
species_elevcode_tally <- data.frame(table(seedlings_comp$Species, seedlings_comp$Elevation.code) [,])
stargazer(species_elevcode_tally, type = "latex", summary = F)

# genus level migration rates plot
ggplot(Genus_mig_rates, aes(x = Genus, y = Migration.rate.m.yr.abund, fill = Genus)) + 
geom_bar(stat = "identity")

## Site level environmental variables (Whitaker?)
Site_Code <- camp_loc$Site
Elevation <- camp_loc$elev_mean
Annual_Precip <- c(NA, NA, 3087, 2631, 2631, 2472, 1827, NA, 2318, NA)
Annual_Air_Temp <- c(NA, NA, 20.7, 17.4, 15.8, 16, NA, 14.9, 11.1, 8.9)
Air_Temp_SD <- c(NA, NA, 0.02, 1.5, 1.3, 1.3, NA, 1.0, 1.0, 1.0)
Slope_deg <- c(NA, NA, NA, 22.7, 40.1, 41.8, 18, NA, 21.4, 11.8)
Total_C <- c(NA, NA, 16, 10.5, 26, 31, 37, NA, 28.5, 44.5)
Total_N <- c(NA, NA, 1.45, 1, 1.85, 2, 2.1, NA, 1.75, 2.6)
Soil_pH <- c(NA, NA, 3.85, 4.05, 4.25, 4.3, 4, NA, 3.9, 3.75)
Site_Char <- data.frame(Site_Code, Elevation, Annual_Precip, Annual_Air_Temp, Slope_deg, Total_C, Total_N, Soil_pH)
stargazer(Site_Char, summary = F, type = "latex", column.labels = c("Site Code", "Elevation (m.a.s.l.)", "Annual Precipitation (mm year-1)", "Annual Air Temperature (C)", "Slope (\textdegrees)", "Total C (%)", "Total N (%)", "Soil pH"), rownames = F)

# Plotting ranges and sample locations
ggplot() + 
	geom_abline(aes(intercept = elev_mean, slope = 0), linetype = 2, size = 0.5, data = camp_loc) +
	geom_point(aes(x = species, y = range, colour = position), size = 3, data = Species_all_loc) + 
	geom_point(aes(x = species, y = range), shape = 15, size = 4, data = ranges) + 
	geom_segment(aes(x = species, y = min, xend = species, yend = max), data = ranges_spread) + 
	xlab("Species") + 
	ylab("Elevation (m)")


# Variables vs Elevation Boxplots ----

seedlings_melt_env <- melt(seedlings_comp, id.vars=c("Individual.code", "Elevation.code", "Species", "Elevation"), 
													 measure.vars = c("Soil.temp.mean", "Soil.mois.mean", "LAI.4.ring", "Comp.seed.total", "Comp.seed.same.sp", "Comp.adult.total", "Comp.adult.metric", "Comp.adult.log.metric"))

ggplot(seedlings_melt_env, aes(x=Elevation.code, y=value)) + 
	geom_boxplot() + 
	facet_wrap(~variable, scales = "free")

# Variables vs Species separated by elevation code boxplots ----

ggplot(seedlings_melt_env, aes(x=Species, y=value, fill = interaction(Elevation.code, Species))) + 
	geom_boxplot(position = position_dodge(), outlier.shape = ) + 
	facet_wrap(~variable, scales = "free")

ggplot(seedlings_melt_traits, aes(x=Species, y=value, fill = interaction(Elevation.code, Species))) + 
	geom_boxplot(position = position_dodge(), outlier.shape = ) + 
	facet_wrap(~variable, scales = "free")

# Variables vs elevation scatterplots ----

ggplot(seedlings_melt_env, aes(x=Elevation, y=value, colour = Species)) + 
	geom_point() + 
	facet_wrap(~variable, scales = "free")

ggplot(seedlings_melt_env, aes(x=Elevation, y=value, colour = Species)) + 
	geom_point() + 
	geom_smooth(aes(fill = Species, colour = Species), method = lm, se = T) + 
	facet_wrap(~variable, scales = "free")

ggplot(seedlings_melt_traits, aes(x=Elevation, y=value, colour = Species)) + 
	geom_point() + 
	geom_smooth(aes(fill = Species, colour = Species), method = lm, se = T) + 
	facet_wrap(~variable, scales = "free")

# Competition radius calculation ----

## Do the calculation
k <- rep(2, 10)
trees_ha <- c(475, 690, 645, 860, 887, 954, 1060, 1101, 1287, 1417)
Competition_Radius <- data.frame("site" = camp_loc$Site, "elev_mean" = camp_loc$elev_mean, k, trees_ha)
Competition_Radius$C_R <- Competition_Radius$k * (sqrt(10000 / trees_ha))

## Linear regression of trees / ha vs elevation to get VC
Competition_Radius_VC <- filter(Competition_Radius,  site %in% c("PA400", "PA800", "SP1500", "SP1750", "TRU08", "TRU07", "TRU06", "TRU04", "TRU02"))

trees_ha_elev <- lm(trees_ha~elev_mean, data = Competition_Radius_VC)
summary(trees_ha_elev)

## Make a plot of with a linear regression
ggplot(Competition_Radius, aes(x = elev_mean, y = trees_ha)) + 
	geom_point() +
	geom_smooth(method = lm) + 
	geom_text(aes(label = site), 
						colour = c(rep("black", 2),"red",rep("black", 7)),
						hjust =-0.2, size = 4) + 
	scale_x_continuous(limits = c(0, 3600)) + 
	ylab(expression(Trees~ha^"-1")) + xlab("Elevation (m.a.s.l.)")

## Pretty table
Competition_Radius$C_R <- signif(Competition_Radius$C_R, digits = 1)
Competition_Radius_K <- dplyr::select(Competition_Radius, -k, -elev_mean)
stargazer(Competition_Radius_K, type = "latex", summary = F)

# Plotting rank abundance curve with our species highlighted ----

## Create summary data frame 
aberg_census_summ <- aberg_census %>%
	mutate(genus_species = paste(genus, specie)) %>%
	group_by(genus_species) %>%
	summarise(n = n()) %>%
	arrange(n) %>%
	mutate(id = seq(from =  length(.$n), to = 1, by = -1),
				 sampled = if_else(genus_species == "Alzatea verticillata" | 
				 										genus_species == "Tapirira guianensis" | 
				 										genus_species == "Clethra revoluta" | 
				 										genus_species == "Clusia thurifera" | 
				 										genus_species == "Hedyosmum goudotianum" | 
				 										genus_species == "Schefflera patula" |
				 										genus_species == "Iriartea deltoidea" |
				 										genus_species == "Dictyocaryum lamarckianum",
				 									T, F))
# NEED TO FIGURE OUT MYRCIA COLOURING
ggplot() + 
	geom_point(data = filter(aberg_census_summ, sampled == T), aes(x = id, y = n), colour = "#E03A50FE", alpha = 1, size = 5) +
	geom_point(data = filter(aberg_census_summ, sampled == F), aes(x = id, y = n), colour = "black", alpha = 0.5, size = 1) #+
	#geom_point(data = filter(aberg_census_summ, sampled == F), aes(x = id, y = n), colour = "black", alpha = 0.5, size = 1) +
	
# Soil Carbon and nitrogen content ----

soil$N <- soil$N....
soil$C <- soil$C....
soil <- select(soil, Site, ID.number, Year, Replicate, N, C)

ggplot(soil, aes(x = Site, y = N)) + geom_boxplot(aes(colour = Site))
ggplot(soil, aes(x = Site, y = C)) + geom_boxplot(aes(colour = Site))

# Linear mixed models, competition versus elevation ----

## LAI.4.ring
LAI4ring_vs_elev_lmer <- lmer(LAI.4.ring ~ Elevation + (1|Site), data = seedlings_rem_na, REML = F)
LAI4ring_vs_elev_lmer_slope <- lmer(LAI.4.ring ~ Elevation + (Elevation|Site), data = seedlings_rem_na, REML = F)
LAI4ring_vs_elev_lmer_rand <- lmer(LAI.4.ring ~ (1|Site), data = seedlings_rem_na, REML = F)
LAI4ring_vs_elev_lmer_null <- lm(LAI.4.ring ~ 1, data = seedlings_rem_na)

## Total seedling abundance
comp_seed_vs_elev_lmer <- lmer(Comp.seed.total ~ Elevation + (Elevation|Site), data = seedlings_rem_na, REML = F)
comp_seed_vs_elev_lmer_slope <- lmer(Comp.seed.total ~ Elevation + (Elevation|Site), data = seedlings_rem_na, REML = F)
comp_seed_vs_elev_glmer_pois <- glmer(Comp.seed.total ~ Elevation + (1|Site), data = seedlings_rem_na, family = poisson, REML = F)
comp_seed_vs_elev_glmer_negbi <- glmer(Comp.seed.total ~ Elevation + (1|Site), data = seedlings_D.FvFm, REML = F, family = negative.binomial)
comp_seed_vs_elev_glmer_rand <- glmer(Comp.seed.total ~ (1|Site), data = seedlings_rem_na, family = poisson, REML = F)
comp_seed_vs_elev_lmer_null <- lm(Comp.seed.total ~ 1, data = seedlings_rem_na)

## Conspecific seedling abundance
comp_seed_same_vs_elev_glmmADMB_pois <- glmmadmb(Comp.seed.same.sp ~ Elevation + (Species|Site), data = seedlings_rem_na, family = "poisson", zeroInflation = T)
comp_seed_same_vs_elev_glmmADMB_nbinom <- glmmadmb(Comp.seed.same.sp ~ Elevation + (Species|Site), data = seedlings_rem_na, family = "nbinom", zeroInflation = T)
comp_seed_same_vs_elev_glmmADMB_rand <- glmmadmb(Comp.seed.same.sp ~ (1|Site), data = seedlings_rem_na, family = "nbinom", zeroInflation = T)

## Adult competition metric
comp_adult_metric_vs_elev_lmer <- lmer(Comp.adult.metric ~ Elevation + (1|Site), data = seedlings_rem_na, REML = F)
comp_adult_log_metric_vs_elev_lmer <- lmer(Comp.adult.log.metric ~ Elevation + (1|Site), data = seedlings_rem_na, REML = F)
comp_adult_log_metric_vs_elev_lmer_slope <- lmer(Comp.adult.log.metric ~ Elevation + (Elevation|Site), data = seedlings_rem_na, REML = F)
comp_adult_log_metric_vs_elev_lmer_null <- lm(Comp.adult.log.metric ~ 1, data = seedlings_rem_na)
comp_adult_log_metric_vs_elev_lmer_rand <- lmer(Comp.adult.log.metric ~ (1|Site), data = seedlings_rem_na)

## Soil temperature
soil_temp_vs_elev_lmer <- lmer(Soil.temp.mean ~ Elevation + (1|Site), data = seedlings_rem_na, REML = F)
soil_temp_vs_elev_lmer_day <- lmer(Soil.temp.mean ~ Elevation + (1|Site) + (1|Collection.day), data = seedlings_rem_na, REML = F)

## Soil moisture
soil_mois_vs_elev_lmer <-  lmer(Soil.mois.mean ~ Elevation + (1|Site), data=seedlings_rem_na, REML = F)
soil_mois_vs_elev_lmer_null <- lmer(Soil.temp.mean ~ (1|Site), data = seedlings_rem_na, REML = F)
soil_mois_vs_elev_date_lmer <- lmer(Soil.mois.mean ~ Elevation + (1|Site) + (1|Collection.date), data = seedlings_rem_na, REML = F)

stargazer(LAI4ring_vs_elev_lmer, comp_seed_vs_elev_glmer_negbi, comp_adult_log_metric_vs_elev_lmer, summary = F)   

# Linear model fits for env. vs elev.

## LAI.4.ring
LAI4ring_pred <- expand.grid(Elevation = c(seq(from=378, to=3224, by=1)), Site = c("PA400", "PA800", "VC", "SP1500", "SP1750", "TRU08", "TRU07", "TRU06", "TRU04", "TRU02"), LAI.4.ring = 0)
matrix_LAI4ring = model.matrix(terms(LAI4ring_vs_elev_lmer), data = LAI4ring_pred)
LAI4ring_pred$LAI.4.ring = matrix_LAI4ring %*% fixef(LAI4ring_vs_elev_lmer)
LAI4ring_var <- diag(matrix_LAI4ring %*% tcrossprod(vcov(LAI4ring_vs_elev_lmer), matrix_LAI4ring))
LAI4ring_pred <- data.frame(LAI4ring_pred, plo = LAI4ring_pred$LAI.4.ring-2*sqrt(LAI4ring_var), phi = LAI4ring_pred$LAI.4.ring+2*sqrt(LAI4ring_var))

LAI4ring_plot <- ggplot(seedlings_rem_na, aes(x=Elevation, y=LAI.4.ring)) + 
	geom_point(aes(colour=Site), alpha = 0.8) +
	layer(data = LAI4ring_pred, geom = "line", mapping=aes(x=Elevation, y=LAI.4.ring), stat = "identity", position = "identity", params=list(na.rm = F)) + 
	layer(data = LAI4ring_pred, geom = "ribbon", mapping = aes(x=Elevation, ymin=plo, ymax=phi), stat="identity", position="identity", params=list(na.rm = F, alpha = 0.6))

## Total seedling abundance
Total_seed_pred <- expand.grid(Elevation = c(seq(from=378, to=3224, by=1)), Site = c("PA400", "PA800", "VC", "SP1500", "SP1750", "TRU08", "TRU07", "TRU06", "TRU04", "TRU02"), Comp.seed.total = 0)
matrix_Total_seed = model.matrix(terms(comp_seed_vs_elev_lmer), data = Total_seed_pred)
Total_seed_pred$Comp.seed.total = matrix_Total_seed %*% fixef(comp_seed_vs_elev_lmer)
Total_seed_var <- diag(matrix_Total_seed %*% tcrossprod(vcov(comp_seed_vs_elev_lmer), matrix_Total_seed))
Total_seed_pred <- data.frame(Total_seed_pred, plo = Total_seed_pred$Comp.seed.total-2*sqrt(Total_seed_var), phi = Total_seed_pred$Comp.seed.total+2*sqrt(Total_seed_var))

Total_seed_plot <- ggplot(seedlings_rem_na, aes(x = Elevation, y = Comp.seed.total)) + 
	geom_point(aes(colour=Site), alpha = 0.8) +
	layer(data = Total_seed_pred, geom='line', mapping = aes(x=Elevation, y = Comp.seed.total), stat = "identity", position = "identity", params = list(na.rm = F)) + 
	layer(data = Total_seed_pred, geom = "ribbon", mapping = aes(x=Elevation, ymin = plo, ymax = phi), stat = "identity", position = "identity", params = list(na.rm = F, alpha = 0.6))

## Conspecific seedling abundance 

Same_seed_pred <- expand.grid(Elevation = c(seq(from=378, to=3224, by=1)), Site = c("PA400", "PA800", "VC", "SP1500", "SP1750", "TRU08", "TRU07", "TRU06", "TRU04", "TRU02"), Comp.seed.same.sp = 0)
matrix_Same_seed = model.matrix(terms(comp_seed_same_vs_elev_glmmADMB_nbinom), data = Same_seed_pred)
Same_seed_pred$Comp.seed.same.sp = matrix_Same_seed %*% fixef(comp_seed_same_vs_elev_glmmADMB_nbinom)
Same_seed_var <- diag(matrix_Same_seed %*% tcrossprod(vcov(comp_seed_same_vs_elev_glmmADMB_nbinom), matrix_Same_seed))
Same_seed_pred <- data.frame(Same_seed_pred, plo = Same_seed_pred$Comp.seed.same.sp-2*sqrt(Same_seed_var), phi = Same_seed_pred$Comp.seed.same.sp+2*sqrt(Same_seed_var))

Same_seed_plot <- ggplot(seedlings_rem_na, aes(x = Elevation, y = Comp.seed.same.sp)) + 
	geom_point(aes(colour=Site), alpha = 0.8) +
	layer(data = Same_seed_pred, geom='line', mapping = aes(x=Elevation, y = Comp.seed.same.sp), stat = "identity", position = "identity", params = list(na.rm = F)) + 
	layer(data = Same_seed_pred, geom = "ribbon", mapping = aes(x=Elevation, ymin = plo, ymax = phi), stat = "identity", position = "identity", params = list(na.rm = F, alpha = 0.6))

## Adult competition metric

comp_adult_log_metric_pred <- expand.grid(Elevation = c(seq(from=378, to=3224, by=1)), Site = c("PA400", "PA800", "VC", "SP1500", "SP1750", "TRU08", "TRU07", "TRU06", "TRU04", "TRU02"), Comp.adult.log.metric = 0)
matrix_comp_adult_log_metric = model.matrix(terms(comp_adult_log_metric_vs_elev_lmer), data = comp_adult_log_metric_pred)
comp_adult_log_metric_pred$comp_adult_log_metric = matrix_comp_adult_log_metric %*% fixef(comp_adult_log_metric_vs_elev_lmer)
comp_adult_log_var <- diag(matrix_comp_adult_log_metric %*% tcrossprod(vcov(comp_adult_log_metric_vs_elev_lmer), matrix_comp_adult_log_metric))
comp_adult_log_metric_pred <- data.frame(comp_adult_log_metric_pred, plo = comp_adult_log_metric_pred$comp_adult_log_metric-2*sqrt(comp_adult_log_var), phi = comp_adult_log_metric_pred$comp_adult_log_metric+2*sqrt(comp_adult_log_var))

isi_plot <- ggplot(seedlings_rem_na, aes(x = Elevation, y = Comp.adult.log.metric)) + 
	geom_point(aes(colour=Site), alpha = 0.8) +
	layer(data = comp_adult_log_metric_pred, geom='line', mapping = aes(x=Elevation, y = (plo + phi) / 2 ), stat = "identity", position = "identity", params = list(na.rm = F)) + 
	layer(data = comp_adult_log_metric_pred, geom = "ribbon", mapping = aes(x=Elevation, ymin = plo, ymax = phi), stat = "identity", position = "identity", params = list(na.rm = F, alpha = 0.6))

## Soil temperature

soil_temp_pred <- expand.grid(Elevation = c(seq(from=378, to=3224, by=1)), Site = c("PA400", "PA800", "VC", "SP1500", "SP1750", "TRU08", "TRU07", "TRU06", "TRU04", "TRU02"), Soil.temp.mean = 0)
matrix_soil_temp = model.matrix(terms(soil_temp_vs_elev_lmer), data = soil_temp_pred)
soil_temp_pred$soil_temp = matrix_soil_temp %*% fixef(soil_temp_vs_elev_lmer)
soil_temp_var <- diag(matrix_soil_temp %*% tcrossprod(vcov(soil_temp_vs_elev_lmer), matrix_soil_temp))
soil_temp_pred <- data.frame(soil_temp_pred, plo = soil_temp_pred$soil_temp-2*sqrt(soil_temp_var), phi = soil_temp_pred$soil_temp+2*sqrt(soil_temp_var))

temp_plot <- ggplot(seedlings_rem_na, aes(x = Elevation, y = Soil.temp.mean)) + 
	geom_point(aes(colour=Site), alpha = 0.8) +
	layer(data = soil_temp_pred, geom='line', mapping = aes(x=Elevation, y = (plo + phi) / 2 ), stat = "identity", position = "identity", params = list(na.rm = F)) + 
	layer(data = soil_temp_pred, geom = "ribbon", mapping = aes(x=Elevation, ymin = plo, ymax = phi), stat = "identity", position = "identity", params = list(na.rm = F, alpha = 0.6))

## Soil moisture

soil_mois_pred <- expand.grid(Elevation = c(seq(from=378, to=3224, by=1)), Site = c("PA400", "PA800", "VC", "SP1500", "SP1750", "TRU08", "TRU07", "TRU06", "TRU04", "TRU02"), Soil.mois.mean = 0)
matrix_soil_mois <- model.matrix(terms(soil_mois_vs_elev_lmer), data = soil_mois_pred)
soil_mois_pred$soil_mois = matrix_soil_mois %*% fixef(soil_mois_vs_elev_lmer)
soil_mois_var <- diag(matrix_soil_mois %*% tcrossprod(vcov(soil_mois_vs_elev_lmer), matrix_soil_mois))
soil_mois_pred <- data.frame(soil_mois_pred, plo = soil_mois_pred$soil_mois-2*sqrt(soil_mois_var), phi = soil_mois_pred$soil_mois+2*sqrt(soil_mois_var))

mois_plot <- ggplot(seedlings_rem_na, aes(x = Elevation, y = Soil.mois.mean)) + 
	geom_point(aes(colour=Site), alpha = 0.8) +
	layer(data = soil_mois_pred, geom='line', mapping = aes(x=Elevation, y = (plo + phi) / 2 ), stat = "identity", position = "identity", params = list(na.rm = F)) + 
	layer(data = soil_mois_pred, geom = "ribbon", mapping = aes(x=Elevation, ymin = plo, ymax = phi), stat = "identity", position = "identity", params = list(na.rm = F, alpha = 0.6))
mois_plot

# Interval plots of linear model slopes for each species

fvfm_slope <- seedlings_rem_na %>%
	group_by(Species) %>%
	do(mod = lm(D.FvFm ~ Elevation, data = .)) %>%
	mutate(Slope = summary(mod)$coeff[2],
				 SE = summary(mod)$coeff[, 2][2],
				 trait = rep("D. Fv/Fm")) %>%
	select(-mod)

spad_slope <- seedlings_rem_na %>%
	group_by(Species) %>%
	do(mod = lm(SPAD.mean ~ Elevation, data = .)) %>%
	mutate(Slope = summary(mod)$coeff[2],
				 SE = summary(mod)$coeff[, 2][2],
				 trait = rep("SPAD")) %>%
	select(-mod)

leaf_area_slope <- seedlings_rem_na %>%
	group_by(Species) %>%
	do(mod = lm(log(Leaf.area) ~ Elevation, data = .)) %>%
	mutate(Slope = summary(mod)$coeff[2],
				 SE = summary(mod)$coeff[, 2][2],
				 trait = rep("Leaf Area")) %>%
	select(-mod)

thick_slope <- seedlings_rem_na %>%
	group_by(Species) %>%
	do(mod = lm(Leaf.thick.mean.mm ~ Elevation, data = .)) %>%
	mutate(Slope = summary(mod)$coeff[2],
				 SE = summary(mod)$coeff[, 2][2],
				 trait = rep("Leaf thickness")) %>%
	select(-mod)

hl_ratio_slope <- seedlings_rem_na %>%
	group_by(Species) %>%
	do(mod = lm(Height.leaf.ratio ~ Elevation, data = .)) %>%
	mutate(Slope = summary(mod)$coeff[2],
				 SE = summary(mod)$coeff[, 2][2],
				 trait = rep("Height:leaf ratio")) %>%
	select(-mod)

stemvol_slope <- seedlings_rem_na %>%
	group_by(Species) %>%
	do(mod = lm(Stem.volume.cm3 ~ Elevation, data = .)) %>%
	mutate(Slope = summary(mod)$coeff[2],
				 SE = summary(mod)$coeff[, 2][2],
				 trait = rep("Stem volume")) %>%
	select(-mod)

lm_sp_slope <- rbind(fvfm_slope, 
										 spad_slope, 
										 leaf_area_slope, 
										 thick_slope,
										 hl_ratio_slope, 
										 stemvol_slope)

ggplot(lm_sp_slope, aes(x = Species)) + 
	geom_point(aes(y = Slope, colour = Species), size = 2) +
	geom_errorbar(aes(ymin = Slope - SE, ymax = Slope + SE, colour = Species), width = 1) +
	facet_wrap(~trait, scales = "free")

# Comparing random intercept and random slope models 

stdz.mod_fvfm_elev_ri <- standardize(lmer(D.FvFm ~ Elevation + (1|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)
stdz.mod_fvfm_elev_rs <- standardize(lmer(D.FvFm ~ Elevation + (Elevation|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)

stdz.mod_fvfm_lai_ri <- standardize(lmer(D.FvFm ~ LAI.4.ring + (1|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)
stdz.mod_fvfm_lai_rs <- standardize(lmer(D.FvFm ~ LAI.4.ring + (LAI.4.ring|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)

stdz.mod_fvfm_seed_ri <- standardize(lmer(D.FvFm ~ Comp.seed.total + (1|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)
stdz.mod_fvfm_seed_rs <- standardize(lmer(D.FvFm ~ Comp.seed.total + (Comp.seed.total|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)

stdz.mod_fvfm_isi_ri <- standardize(lmer(D.FvFm ~ Comp.adult.log.metric + (1|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)
stdz.mod_fvfm_isi_rs <- standardize(lmer(D.FvFm ~ Comp.adult.log.metric + (Comp.adult.log.metric|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)

stdz.mod_spad_elev_ri <- standardize(lmer(SPAD.mean ~ Elevation + (1|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)
stdz.mod_spad_elev_rs <- standardize(lmer(SPAD.mean ~ Elevation + (Elevation|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)

stdz.mod_spad_lai_ri <- standardize(lmer(SPAD.mean ~ LAI.4.ring + (1|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)
stdz.mod_spad_lai_rs <- standardize(lmer(SPAD.mean ~ LAI.4.ring + (LAI.4.ring|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)

stdz.mod_spad_seed_ri <- standardize(lmer(SPAD.mean ~ Comp.seed.total + (1|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)
stdz.mod_spad_seed_rs <- standardize(lmer(SPAD.mean ~ Comp.seed.total + (Comp.seed.total|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)

stdz.mod_spad_isi_ri <- standardize(lmer(SPAD.mean ~ Comp.adult.log.metric + (1|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)
stdz.mod_spad_isi_rs <- standardize(lmer(SPAD.mean ~ Comp.adult.log.metric + (Comp.adult.log.metric|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)

stdz.mod_thick_elev_ri <- standardize(lmer(Leaf.thick.mean.mm ~ Elevation + (1|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)
stdz.mod_thick_elev_rs <- standardize(lmer(Leaf.thick.mean.mm ~ Elevation + (Elevation|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)

stdz.mod_thick_lai_ri <- standardize(lmer(Leaf.thick.mean.mm ~ LAI.4.ring + (1|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)
stdz.mod_thick_lai_rs <- standardize(lmer(Leaf.thick.mean.mm ~ LAI.4.ring + (LAI.4.ring|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)

stdz.mod_thick_seed_ri <- standardize(lmer(Leaf.thick.mean.mm ~ Comp.seed.total + (1|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)
stdz.mod_thick_seed_rs <- standardize(lmer(Leaf.thick.mean.mm ~ Comp.seed.total + (Comp.seed.total|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)

stdz.mod_thick_isi_ri <- standardize(lmer(Leaf.thick.mean.mm ~ Comp.adult.log.metric + (1|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)
stdz.mod_thick_isi_rs <- standardize(lmer(Leaf.thick.mean.mm ~ Comp.adult.log.metric + (Comp.adult.log.metric|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)

stdz.mod_hlratio_elev_ri <- standardize(lmer(Height.leaf.ratio ~ Elevation + (1|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)
stdz.mod_hlratio_elev_rs <- standardize(lmer(Height.leaf.ratio ~ Elevation + (Elevation|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)

stdz.mod_hlratio_lai_ri <- standardize(lmer(Height.leaf.ratio ~ LAI.4.ring + (1|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)
stdz.mod_hlratio_lai_rs <- standardize(lmer(Height.leaf.ratio ~ LAI.4.ring + (LAI.4.ring|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)

stdz.mod_hlratio_seed_ri <- standardize(lmer(Height.leaf.ratio ~ Comp.seed.total + (1|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)
stdz.mod_hlratio_seed_rs <- standardize(lmer(Height.leaf.ratio ~ Comp.seed.total + (Comp.seed.total|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)

stdz.mod_hlratio_isi_ri <- standardize(lmer(Height.leaf.ratio ~ Comp.adult.log.metric + (1|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)
stdz.mod_hlratio_isi_rs <- standardize(lmer(Height.leaf.ratio ~ Comp.adult.log.metric + (Comp.adult.log.metric|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)

stdz.mod_area_elev_ri <- standardize(lmer(Leaf.area ~ Elevation + (1|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)
stdz.mod_area_elev_rs <- standardize(lmer(Leaf.area ~ Elevation + (Elevation|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)

stdz.mod_area_lai_ri <- standardize(lmer(Leaf.area ~ LAI.4.ring + (1|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)
stdz.mod_area_lai_rs <- standardize(lmer(Leaf.area ~ LAI.4.ring + (LAI.4.ring|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)

stdz.mod_area_seed_ri <- standardize(lmer(Leaf.area ~ Comp.seed.total + (1|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)
stdz.mod_area_seed_rs <- standardize(lmer(Leaf.area ~ Comp.seed.total + (Comp.seed.total|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)

stdz.mod_area_isi_ri <- standardize(lmer(Leaf.area ~ Comp.adult.log.metric + (1|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)
stdz.mod_area_isi_rs <- standardize(lmer(Leaf.area ~ Comp.adult.log.metric + (Comp.adult.log.metric|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)

stdz.mod_stemvol_elev_ri <- standardize(lmer(Stem.volume.cm3 ~ Elevation + (1|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)
stdz.mod_stemvol_elev_rs <- standardize(lmer(Stem.volume.cm3 ~ Elevation + (Elevation|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)

stdz.mod_stemvol_lai_ri <- standardize(lmer(Stem.volume.cm3 ~ LAI.4.ring + (1|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)
stdz.mod_stemvol_lai_rs <- standardize(lmer(Stem.volume.cm3 ~ LAI.4.ring + (LAI.4.ring|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)

stdz.mod_stemvol_seed_ri <- standardize(lmer(Stem.volume.cm3 ~ Comp.seed.total + (1|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)
stdz.mod_stemvol_seed_rs <- standardize(lmer(Stem.volume.cm3 ~ Comp.seed.total + (Comp.seed.total|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)

stdz.mod_stemvol_isi_ri <- standardize(lmer(Stem.volume.cm3 ~ Comp.adult.log.metric + (1|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)
stdz.mod_stemvol_isi_rs <- standardize(lmer(Stem.volume.cm3 ~ Comp.adult.log.metric + (Comp.adult.log.metric|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)

stdz.mod_fvfm_null <- standardize(lmer(D.FvFm ~ (1|Species) + (1|Site), data=seedlings_rem_na), standardize.y = T)
stdz.mod_spad_null <- standardize(lmer(SPAD.mean ~ (1|Species) + (1|Site), data=seedlings_rem_na), standardize.y = T)
stdz.mod_thick_null <- standardize(lmer(Leaf.thick.mean.mm ~ (1|Species) + (1|Site), data=seedlings_rem_na), standardize.y = T)
stdz.mod_hlratio_null <- standardize(lmer(Height.leaf.ratio ~ (1|Species) + (1|Site), data=seedlings_rem_na), standardize.y = T)
stdz.mod_area_null <- standardize(lmer(Leaf.area ~ (1|Species) + (1|Site), data=seedlings_rem_na), standardize.y = T)
stdz.mod_stemvol_null <- standardize(lmer(Stem.volume.cm3 ~ (1|Species) + (1|Site), data=seedlings_rem_na), standardize.y = T)
AIC_ri <- c(AIC(stdz.mod_fvfm_elev_ri),
						AIC(stdz.mod_fvfm_lai_ri),
						AIC(stdz.mod_fvfm_seed_ri),
						AIC(stdz.mod_fvfm_isi_ri),
						AIC(stdz.mod_spad_elev_ri),
						AIC(stdz.mod_spad_lai_ri),
						AIC(stdz.mod_spad_seed_ri),
						AIC(stdz.mod_spad_isi_ri),
						AIC(stdz.mod_thick_elev_ri),
						AIC(stdz.mod_thick_lai_ri),
						AIC(stdz.mod_thick_seed_ri),
						AIC(stdz.mod_thick_isi_ri),
						AIC(stdz.mod_hlratio_elev_ri),
						AIC(stdz.mod_hlratio_lai_ri),
						AIC(stdz.mod_hlratio_seed_ri),
						AIC(stdz.mod_hlratio_isi_ri),
						AIC(stdz.mod_area_elev_ri),
						AIC(stdz.mod_area_lai_ri),
						AIC(stdz.mod_area_seed_ri),
						AIC(stdz.mod_area_isi_ri),
						AIC(stdz.mod_stemvol_elev_ri),
						AIC(stdz.mod_stemvol_lai_ri),
						AIC(stdz.mod_stemvol_seed_ri),
						AIC(stdz.mod_stemvol_isi_ri))

AIC_rs <- c(AIC(stdz.mod_fvfm_elev_rs),
						AIC(stdz.mod_fvfm_lai_rs),
						AIC(stdz.mod_fvfm_seed_rs),
						AIC(stdz.mod_fvfm_isi_rs),
						AIC(stdz.mod_spad_elev_rs),
						AIC(stdz.mod_spad_lai_rs),
						AIC(stdz.mod_spad_seed_rs),
						AIC(stdz.mod_spad_isi_rs),
						AIC(stdz.mod_thick_elev_rs),
						AIC(stdz.mod_thick_lai_rs),
						AIC(stdz.mod_thick_seed_rs),
						AIC(stdz.mod_thick_isi_rs),
						AIC(stdz.mod_hlratio_elev_rs),
						AIC(stdz.mod_hlratio_lai_rs),
						AIC(stdz.mod_hlratio_seed_rs),
						AIC(stdz.mod_hlratio_isi_rs),
						AIC(stdz.mod_area_elev_rs),
						AIC(stdz.mod_area_lai_rs),
						AIC(stdz.mod_area_seed_rs),
						AIC(stdz.mod_area_isi_rs),
						AIC(stdz.mod_stemvol_elev_rs),
						AIC(stdz.mod_stemvol_lai_rs),
						AIC(stdz.mod_stemvol_seed_rs),
						AIC(stdz.mod_stemvol_isi_rs))

AIC_ri - AIC_rs

daic1 <-  AIC(stdz.mod_fvfm_null) - AIC(stdz.mod_fvfm_elev_ri) 
daic2 <-  AIC(stdz.mod_fvfm_null) - AIC(stdz.mod_fvfm_lai_ri) 
daic3 <-  AIC(stdz.mod_fvfm_null) - AIC(stdz.mod_fvfm_seed_ri) 
daic4 <-  AIC(stdz.mod_fvfm_null) - AIC(stdz.mod_fvfm_isi_ri) 
daic5 <-  AIC(stdz.mod_spad_null) - AIC(stdz.mod_spad_elev_ri) 
daic6 <-  AIC(stdz.mod_spad_null) - AIC(stdz.mod_spad_lai_ri) 
daic7 <-  AIC(stdz.mod_spad_null) - AIC(stdz.mod_spad_seed_ri) 
daic8 <-  AIC(stdz.mod_spad_null) - AIC(stdz.mod_spad_isi_ri) 
daic9 <-  AIC(stdz.mod_thick_null) - AIC(stdz.mod_thick_elev_ri) 
daic10 <-  AIC(stdz.mod_thick_null) - AIC(stdz.mod_thick_lai_rs) 
daic11<-  AIC(stdz.mod_thick_null) - AIC(stdz.mod_thick_seed_ri) 
daic12<-  AIC(stdz.mod_thick_null) -  AIC(stdz.mod_thick_isi_ri) 
daic13<-  AIC(stdz.mod_hlratio_null) - AIC(stdz.mod_hlratio_elev_rs) 
daic14<-  AIC(stdz.mod_hlratio_null) - AIC(stdz.mod_hlratio_lai_ri) 
daic15<-  AIC(stdz.mod_hlratio_null) - AIC(stdz.mod_hlratio_seed_rs) 
daic16<-  AIC(stdz.mod_hlratio_null) - AIC(stdz.mod_hlratio_isi_ri) 
daic17<-  AIC(stdz.mod_area_null) - AIC(stdz.mod_area_elev_rs)
daic18<-  AIC(stdz.mod_area_null) - AIC(stdz.mod_area_lai_rs) 
daic19<-  AIC(stdz.mod_area_null) - AIC(stdz.mod_area_seed_ri) 
daic20<-  AIC(stdz.mod_area_null) - AIC(stdz.mod_area_isi_rs) 
daic21<-  AIC(stdz.mod_stemvol_null) - AIC(stdz.mod_stemvol_elev_ri) 
daic22<-  AIC(stdz.mod_stemvol_null) - AIC(stdz.mod_stemvol_lai_ri) 
daic23<-  AIC(stdz.mod_stemvol_null) - AIC(stdz.mod_stemvol_seed_ri) 
daic24<-  AIC(stdz.mod_stemvol_null) - AIC(stdz.mod_stemvol_isi_rs) 


r1 <-   unname(r.squaredGLMM(stdz.mod_fvfm_elev_ri))
r2 <-   unname(r.squaredGLMM(stdz.mod_fvfm_lai_ri))
r3 <-   unname(r.squaredGLMM(stdz.mod_fvfm_seed_ri))
r4 <-   unname(r.squaredGLMM(stdz.mod_fvfm_isi_ri))
r5 <-   unname(r.squaredGLMM(stdz.mod_spad_elev_rs))
r6 <-   unname(r.squaredGLMM(stdz.mod_spad_lai_ri))
r7 <-   unname(r.squaredGLMM(stdz.mod_spad_seed_ri))
r8 <-   unname(r.squaredGLMM(stdz.mod_spad_isi_ri))
r9 <-   unname(r.squaredGLMM(stdz.mod_thick_elev_ri))
r10 <-   unname(r.squaredGLMM(stdz.mod_thick_lai_rs))
r11<-   unname(r.squaredGLMM(stdz.mod_thick_seed_ri))
r12<-   unname(r.squaredGLMM(stdz.mod_thick_isi_ri))
r13<-   unname(r.squaredGLMM(stdz.mod_hlratio_elev_rs))
r14<-   unname(r.squaredGLMM(stdz.mod_hlratio_lai_ri))
r15<-   unname(r.squaredGLMM(stdz.mod_hlratio_seed_rs))
r16<-   unname(r.squaredGLMM(stdz.mod_hlratio_isi_ri))
r17<-   unname(r.squaredGLMM(stdz.mod_area_elev_rs))
r18<-   unname(r.squaredGLMM(stdz.mod_area_lai_rs))
r19<-   unname(r.squaredGLMM(stdz.mod_area_seed_ri))
r20<-   unname(r.squaredGLMM(stdz.mod_area_isi_rs))
r21<-   unname(r.squaredGLMM(stdz.mod_stemvol_elev_ri))
r22<-   unname(r.squaredGLMM(stdz.mod_stemvol_lai_ri))
r23<-   unname(r.squaredGLMM(stdz.mod_stemvol_seed_ri))
r24<-   unname(r.squaredGLMM(stdz.mod_stemvol_isi_rs))

slope1 <-   unname(coef(summary(stdz.mod_fvfm_elev_ri)) [, "Estimate"])
slope2 <-   unname(coef(summary(stdz.mod_fvfm_lai_ri)) [, "Estimate"])
slope3 <-   unname(coef(summary(stdz.mod_fvfm_seed_ri)) [, "Estimate"])
slope4 <-   unname(coef(summary(stdz.mod_fvfm_isi_ri)) [, "Estimate"])
slope5 <-   unname(coef(summary(stdz.mod_spad_elev_ri)) [, "Estimate"])
slope6 <-   unname(coef(summary(stdz.mod_spad_lai_ri)) [, "Estimate"])
slope7 <-   unname(coef(summary(stdz.mod_spad_seed_ri)) [, "Estimate"])
slope8 <-   unname(coef(summary(stdz.mod_spad_isi_ri)) [, "Estimate"])
slope9 <-   unname(coef(summary(stdz.mod_thick_elev_ri)) [, "Estimate"])
slope10 <-   unname(coef(summary(stdz.mod_thick_lai_rs)) [, "Estimate"])
slope11<-   unname(coef(summary(stdz.mod_thick_seed_ri)) [, "Estimate"])
slope12<-   unname(coef(summary(stdz.mod_thick_isi_ri)) [, "Estimate"])
slope13<-   unname(coef(summary(stdz.mod_hlratio_elev_rs)) [, "Estimate"])
slope14<-   unname(coef(summary(stdz.mod_hlratio_lai_ri)) [, "Estimate"])
slope15<-   unname(coef(summary(stdz.mod_hlratio_seed_rs)) [, "Estimate"])
slope16<-   unname(coef(summary(stdz.mod_hlratio_isi_ri)) [, "Estimate"])
slope17<-   unname(coef(summary(stdz.mod_area_elev_rs)) [, "Estimate"])
slope18<-   unname(coef(summary(stdz.mod_area_lai_rs)) [, "Estimate"])
slope19<-   unname(coef(summary(stdz.mod_area_seed_ri)) [, "Estimate"])
slope20<-   unname(coef(summary(stdz.mod_area_isi_rs)) [, "Estimate"])
slope21<-   unname(coef(summary(stdz.mod_stemvol_elev_ri)) [, "Estimate"])
slope22<-   unname(coef(summary(stdz.mod_stemvol_lai_ri)) [, "Estimate"])
slope23<-   unname(coef(summary(stdz.mod_stemvol_seed_ri)) [, "Estimate"])
slope24<-   unname(coef(summary(stdz.mod_stemvol_isi_rs)) [, "Estimate"])

ss1 <-   sqrt(diag(vcov(stdz.mod_fvfm_elev_ri)))
ss2 <-   sqrt(diag(vcov(stdz.mod_fvfm_lai_ri)))
ss3 <-   sqrt(diag(vcov(stdz.mod_fvfm_seed_ri)))
ss4 <-   sqrt(diag(vcov(stdz.mod_fvfm_isi_ri)))
ss5 <-   sqrt(diag(vcov(stdz.mod_spad_elev_ri)))
ss6 <-   sqrt(diag(vcov(stdz.mod_spad_lai_ri)))
ss7 <-   sqrt(diag(vcov(stdz.mod_spad_seed_ri)))
ss8 <-   sqrt(diag(vcov(stdz.mod_spad_isi_ri)))
ss9 <-   sqrt(diag(vcov(stdz.mod_thick_elev_ri)))
ss10 <-   sqrt(diag(vcov(stdz.mod_thick_lai_rs)))
ss11<-   sqrt(diag(vcov(stdz.mod_thick_seed_ri)))
ss12<-   sqrt(diag(vcov(stdz.mod_thick_isi_ri)))
ss13<-   sqrt(diag(vcov(stdz.mod_hlratio_elev_rs)))
ss14<-   sqrt(diag(vcov(stdz.mod_hlratio_lai_ri)))
ss15<-   sqrt(diag(vcov(stdz.mod_hlratio_seed_rs)))
ss16<-   sqrt(diag(vcov(stdz.mod_hlratio_isi_ri)))
ss17<-   sqrt(diag(vcov(stdz.mod_area_elev_rs)))
ss18<-   sqrt(diag(vcov(stdz.mod_area_lai_rs)))
ss19<-   sqrt(diag(vcov(stdz.mod_area_seed_ri)))
ss20<-   sqrt(diag(vcov(stdz.mod_area_isi_rs)))
ss21<-   sqrt(diag(vcov(stdz.mod_stemvol_elev_ri)))
ss22<-   sqrt(diag(vcov(stdz.mod_stemvol_lai_ri)))
ss23<-   sqrt(diag(vcov(stdz.mod_stemvol_seed_ri)))
ss24<-   sqrt(diag(vcov(stdz.mod_stemvol_isi_rs)))

model_name_rirs <- c("elev", "lai", "seed", "isi","elev", "lai", "seed", "isi","elev", "lai", "seed", "isi","elev", "lai", "seed", "isi","elev", "lai", "seed", "isi","elev", "lai", "seed", "isi")
delta_aic_rirs <- AIC_ri - AIC_rs
r2c_best <- c(r1[2], r2[2], r3[2], r4[2],r5[2],r6[2],r7[2],r8[2],r9[2],r10[2],r11[2],r12[2],r13[2],r14[2],r15[2],r16[2],r17[2],r18[2],r19[2],r20[2],r21[2],r22[2],r23[2],r24[2])
r2m_best <- c(r1[1], r2[1], r3[1], r4[1],r5[1],r6[1],r7[1],r8[1],r9[1],r10[1],r11[1],r12[1],r13[1],r14[1],r15[1],r16[1],r17[1],r18[1],r19[1],r20[1],r21[1],r22[1],r23[1],r24[1])
slope_best <- c(slope1[2], slope2[2], slope3[2], slope4[2],slope5[2],slope6[2],slope7[2],slope8[2],slope9[2],slope10[2],slope11[2],slope12[2],slope13[2],slope14[2],slope15[2],slope16[2],slope17[2],slope18[2],slope19[2],slope20[2],slope21[2],slope22[2],slope23[2],slope24[2])
ss_best <- c(ss1[2], ss2[2], ss3[2], ss4[2],ss5[2],ss6[2],ss7[2],ss8[2],ss9[2],ss10[2],ss11[2],ss12[2],ss13[2],ss14[2],ss15[2],ss16[2],ss17[2],ss18[2],ss19[2],ss20[2],ss21[2],ss22[2],ss23[2],ss24[2])
dAIC_best <- c(daic1 , daic2 , daic3 , daic4 ,daic5 ,daic6 ,daic7 ,daic8 ,daic9 ,daic10 ,daic11 ,daic12 ,daic13 ,daic14 ,daic15 ,daic16 ,daic17 ,daic18 ,daic19 ,daic20 ,daic21 ,daic22 ,daic23 ,daic24 )
rsri <- c("RI","RI","RI","RI","RI","RS","RI","RS","RI","RS","RI","RI","RI","RI","RI","RS","RI","RS","RS","RS","RI","RI","RI","RS")
stdz.mod_fvfm_cov <- data.frame(model_name_rirs, delta_aic_rirs, rsri, dAIC_best, slope_best, ss_best, r2c_best)
stargazer(stdz.mod_fvfm_cov, summary = F)

# Plot of dAIC, Slope and R2 for all best fit lmers 

trait_name_best <- c("D.FvFm","D.FvFm","D.FvFm","D.FvFm",
										 "SPAD","SPAD","SPAD","SPAD",
										 "Leaf Thickness","Leaf Thickness","Leaf Thickness","Leaf Thickness",
										 "Height:Leaf Ratio", "Height:Leaf Ratio", "Height:Leaf Ratio", "Height:Leaf Ratio",
										 "Leaf Area","Leaf Area","Leaf Area","Leaf Area",
										 "Stem Volume", "Stem Volume", "Stem Volume", "Stem Volume")
model_name_best <- c("elev", "lai", "seed", "isi","elev", "lai", "seed", "isi","elev", "lai", "seed", "isi","elev", "lai", "seed", "isi","elev", "lai", "seed", "isi","elev", "lai", "seed", "isi")
rsri_best <- c("RI","RI","RI","RI","RI","RS","RI","RS","RI","RS","RI","RI","RI","RI","RI","RS","RI","RS","RS","RS","RI","RI","RI","RS")
slope_best
ss_best
dAIC_best
r2c_best
explan_best <- data.frame(trait_name_best, model_name_best, slope_best, ss_best, dAIC_best, r2c_best, r2m_best)
explan_best
explan_best$trait_name_best <- factor(explan_best$trait_name_best, levels = explan_best$trait_name_best)

explan_best_phys <- filter(explan_best, trait_name_best %in% c("D.FvFm", "SPAD"))

# dAIC
daicplot<- ggplot(explan_best_phys, aes(x = trait_name_best, y = dAIC_best, group = factor(model_name_best))) +geom_bar(stat = "identity", position = position_dodge(), aes(fill = model_name_best)) +
	theme(legend.position = "right") + 
	theme(legend.title = element_blank()) +
	scale_fill_discrete(breaks = c("elev","isi","lai","seed"), labels = c("Elevation", "ISI", "LAI", "Herbaceous Plant Abundance")) +
	scale_x_discrete(breaks = c("D.FvFm", "SPAD", "Leaf Thickness", "Height:Leaf Ratio", "Leaf Area", "Stem Volume"), labels = c("Photosynthetic\nEfficiency", "Chlorophyll\nContent", "Leaf Thickness", "Leaf:Height Ratio", "Leaf Area", "Stem Volume")) +
	ylab(expression(paste(Delta,"AIC"[r]))) +
	xlab("") +
	geom_hline(aes(yintercept = 2), linetype = 5, colour = "red")  + theme_bw() +theme(plot.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +theme(panel.border= element_blank())+theme(axis.line.x = element_line(color="black", size = 1),axis.line.y = element_line(color="black", size = 1)) +
	theme(axis.text = element_text(size = 25), axis.title = element_text(size = 25))

# R^2
r2plot <- ggplot(explan_best_phys, aes(x = trait_name_best, group = factor(model_name_best)))  +
	geom_bar(stat = "identity", position = position_dodge(), aes(y = r2m_best, fill = model_name_best)) +
	theme(legend.position = "right") + 
	theme(legend.title = element_blank()) +
	scale_fill_discrete(breaks = c("elev","isi","lai","seed"), labels = c("Elevation", "ISI", "LAI", "Herbaceous Plant Abundance")) +
	ylab(expression(paste("Variance Explained (",R[M]^2,")"))) +
	scale_x_discrete(breaks = c("D.FvFm", "SPAD", "Leaf Thickness", "Height:Leaf Ratio", "Leaf Area", "Stem Volume"), labels = c("Photosynthetic\nEfficiency", "Chlorophyll\nContent", "Leaf Thickness", "Leaf:Height Ratio", "Leaf Area", "Stem Volume")) +
	xlab("") + theme_bw() +theme(plot.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +theme(panel.border= element_blank())+theme(axis.line.x = element_line(color="black", size = 1),axis.line.y = element_line(color="black", size = 1)) +
	theme(axis.text = element_text(size = 25), axis.title = element_text(size = 25))

daicplot
r2plot

# Effect size graphs for each trait single predictor model 

ggplot(data = explan_best) + 
	geom_errorbar(aes(x = model_name_best, ymin = slope_best - ss_best, ymax = slope_best + ss_best, colour = factor(model_name_best))) +
	geom_point(aes(x = model_name_best, y = slope_best, colour = factor(model_name_best)), size = 5) + 
	geom_hline(yintercept = 0, linetype = 5) + 
	scale_x_discrete("model_name_best", labels = c("elev" = "Elevation", "isi" = "ISI", "lai" = "LAI", "seed" = "Herbaceous Plants")) +
	xlab("Fixed Effects") + 
	ylab("Fixed Effect Slope") +
	theme(legend.position = "null") +
	facet_wrap(~trait_name_best, scales = "fixed")
	
# Combined predictor linear mixed models ----

# D.FvFm
stdz.mod_fvfm_full <- standardize(lmer(D.FvFm ~ LAI.4.ring + Comp.seed.total + Comp.adult.log.metric + Elevation + (1|Species) + (1|Site) , data=seedlings_rem_na,  REML = F), standardize.y = T)

stdz.mod_fvfm_lai_seed_isi <- standardize(lmer(D.FvFm ~ LAI.4.ring + Comp.seed.total + Comp.adult.log.metric + (1|Species) + (1|Site) , data=seedlings_rem_na,  REML = F), standardize.y = T)

stdz.mod_fvfm_lai_seed_elev <- standardize(lmer(D.FvFm ~ LAI.4.ring + Comp.seed.total + Elevation + (1|Species) + (1|Site) , data=seedlings_rem_na,  REML = F), standardize.y = T)

stdz.mod_fvfm_lai_isi_elev <- standardize(lmer(D.FvFm ~ LAI.4.ring + Comp.adult.log.metric + Elevation + (1|Species) + (1|Site), data=seedlings_rem_na,  REML = F), standardize.y = T)

stdz.mod_fvfm_seed_isi_elev <- standardize(lmer(D.FvFm ~ Comp.seed.total + Comp.adult.log.metric + Elevation + (1|Species) + (1|Site), data=seedlings_rem_na,  REML = F), standardize.y = T)

stdz.mod_fvfm_rand <- standardize(lmer(D.FvFm ~ (1|Species) + (1|Site), data=seedlings_rem_na, REML = F), standardize.y = T)

stdz.mod_fvfm_null <- lm(D.FvFm ~ 1, data=seedlings_rem_na)

stdz.mod_fvfm_aic_val <- c(AIC(stdz.mod_fvfm_full),
													 AIC(stdz.mod_fvfm_lai_seed_isi),
													 AIC(stdz.mod_fvfm_lai_seed_elev),
													 AIC(stdz.mod_fvfm_lai_isi_elev),
													 AIC(stdz.mod_fvfm_seed_isi_elev),
													 AIC(stdz.mod_fvfm_rand),
													 AIC(stdz.mod_fvfm_elev_ri),
													 AIC(stdz.mod_fvfm_lai_ri),
													 AIC(stdz.mod_fvfm_isi_ri),
													 AIC(stdz.mod_fvfm_seed_ri))

model_name_stdz.mod_fvfm_aic <- c("stdz.mod_fvfm_full",
																	"stdz.mod_fvfm_lai_seed_isi",
																	"stdz.mod_fvfm_lai_seed_elev",
																	"stdz.mod_fvfm_lai_isi_elev",
																	"stdz.mod_fvfm_seed_isi_elev",
																	"stdz.mod_fvfm_rand",
																	"stdz.mod_fvfm_elev_ri",
																	"stdz.mod_fvfm_lai_ri",
																	"stdz.mod_fvfm_isi_ri",
																	"stdz.mod_fvfm_seed_ri")

stdz.mod_fvfm_aic <- data.frame(model_name_stdz.mod_fvfm_aic, stdz.mod_fvfm_aic_val)
stdz.mod_fvfm_aic_arr <- arrange(stdz.mod_fvfm_aic, stdz.mod_fvfm_aic_val)
stdz.mod_fvfm_aic_arr
stdz.mod_fvfm_aic_arr$dAIC <- c(AIC(stdz.mod_fvfm_lai_isi_elev) - AIC(stdz.mod_fvfm_rand), 
																AIC(stdz.mod_fvfm_full) - AIC(stdz.mod_fvfm_rand), 
																AIC(stdz.mod_fvfm_seed_isi_elev) - AIC(stdz.mod_fvfm_rand),
																AIC(stdz.mod_fvfm_elev_ri) - AIC(stdz.mod_fvfm_rand),
																AIC(stdz.mod_fvfm_lai_seed_elev) - AIC(stdz.mod_fvfm_rand),
																AIC(stdz.mod_fvfm_isi_ri) - AIC(stdz.mod_fvfm_rand),
																AIC(stdz.mod_fvfm_lai_seed_isi) - AIC(stdz.mod_fvfm_rand),
																AIC(stdz.mod_fvfm_lai_ri) - AIC(stdz.mod_fvfm_rand),
																AIC(stdz.mod_fvfm_rand) - AIC(stdz.mod_fvfm_rand),
																AIC(stdz.mod_fvfm_seed_ri) - AIC(stdz.mod_fvfm_rand))

stdz.mod_fvfm_wi <- akaike.weights(stdz.mod_fvfm_aic_arr$stdz.mod_fvfm_aic_val)
stdz.mod_fvfm_wi
stdz.mod_fvfm_aic_arr$wi <- stdz.mod_fvfm_wi$weights
stdz.mod_fvfm_aic_arr <- stdz.mod_fvfm_aic_arr %>% dplyr::select(model_name_stdz.mod_fvfm_aic, stdz.mod_fvfm_aic_val, dAIC, wi)
r1 <- unname(r.squaredGLMM(stdz.mod_fvfm_lai_isi_elev))
r2 <-   unname(r.squaredGLMM(stdz.mod_fvfm_full))
r3 <-   unname(r.squaredGLMM(stdz.mod_fvfm_seed_isi_elev))
r4 <-   unname(r.squaredGLMM(stdz.mod_fvfm_elev_ri))
r5 <-   unname(r.squaredGLMM(stdz.mod_fvfm_lai_seed_elev))
r6 <-   unname(r.squaredGLMM(stdz.mod_fvfm_isi_ri))
r7 <-   unname(r.squaredGLMM(stdz.mod_fvfm_lai_seed_isi))
r8 <-   unname(r.squaredGLMM(stdz.mod_fvfm_lai_ri))
r9 <-   unname(r.squaredGLMM(stdz.mod_fvfm_rand))
r10 <-   unname(r.squaredGLMM(stdz.mod_fvfm_seed_ri))

stdz.mod_fvfm_aic_arr$R2c <- c(r1[2],r2[2], r3[2], r4[2], r5[2], r6[2], r7[2], r8[2], r9[2], r10[2])
stdz.mod_fvfm_aic_arr$r2m <- c(r1[1], r2[1], r3[1], r4[1], r5[1], r6[1], r7[1], r8[1], r9[1], r10[1])

stdz.mod_fvfm_aic_arr
stargazer(stdz.mod_fvfm_aic_arr, summary = F)

summary(stdz.mod_fvfm_lai_isi_elev)
pamer.fnc(stdz.mod_fvfm_lai_isi_elev)
# Creating fixed effect slope graph with standard error bars
stdz.mod_fvfm_fix_name <- c("LAI", "ISI", "Elevation", "Comp.seed")                      
stdz.mod_fvfm_fix_slope <- c(0.11454, -0.12582, 0.36059, NA)
stdz.mod_fvfm_fix_se <- c(0.06910, 0.06756, 0.10562, NA)
stdz.mod_fvfm_fix_trait <- rep("fvfm", 4)
stdz.mod_fvfm_fix <- data.frame("factor" = stdz.mod_fvfm_fix_name, "slope" = stdz.mod_fvfm_fix_slope, "se" = stdz.mod_fvfm_fix_se, "trait" = stdz.mod_fvfm_fix_trait)
stdz.mod_fvfm_fix
stdz.mod_fvfm_sort <- c("Elevation", "ISI", "LAI", "Comp.seed")
stdz.mod_fvfm_fix <- stdz.mod_fvfm_fix %>% slice(match(stdz.mod_fvfm_sort, stdz.mod_fvfm_fix_name))
stdz.mod_fvfm_fix$stdz.mod_fvfm_fix_name <- factor(stdz.mod_fvfm_fix$stdz.mod_fvfm_fix_name, levels = stdz.mod_fvfm_fix$stdz.mod_fvfm_fix_name)
stdz.mod_fvfm_fix

# SPAD, working out the error structure
#SPAD.mean, Comparison of models, 'ri' was mostly more parsimonious (lower AIC) so using that throughout the models in this section
stdz.mod_spad_full <- standardize(lmer(SPAD.mean ~ LAI.4.ring + Comp.seed.total + Comp.adult.log.metric + Elevation + (1|Species) + (1|Site) , data=seedlings_rem_na,  REML = F), standardize.y = T)

stdz.mod_spad_lai_seed_isi <- standardize(lmer(SPAD.mean ~ LAI.4.ring + Comp.seed.total + Comp.adult.log.metric + (1|Species) + (1|Site) , data=seedlings_rem_na,  REML = F), standardize.y = T)

stdz.mod_spad_lai_seed_elev <- standardize(lmer(SPAD.mean ~ LAI.4.ring + Comp.seed.total + Elevation + (1|Species) + (1|Site) , data=seedlings_rem_na,  REML = F), standardize.y = T)

stdz.mod_spad_lai_isi_elev <- standardize(lmer(SPAD.mean ~ LAI.4.ring + Comp.adult.log.metric + Elevation + (1|Species) + (1|Site), data=seedlings_rem_na,  REML = F), standardize.y = T)

stdz.mod_spad_seed_isi_elev <- standardize(lmer(SPAD.mean ~ Comp.seed.total + Comp.adult.log.metric + Elevation + (1|Species) + (1|Site), data=seedlings_rem_na,  REML = F), standardize.y = T)

stdz.mod_spad_rand <- standardize(lmer(SPAD.mean ~ (1|Species) + (1|Site), data=seedlings_rem_na, REML = F), standardize.y = T)

stdz.mod_spad_null <- lm(SPAD.mean ~ 1, data=seedlings_rem_na)

stdz.mod_spad_aic_val <- c(AIC(stdz.mod_spad_full),
													 AIC(stdz.mod_spad_lai_seed_isi),
													 AIC(stdz.mod_spad_lai_seed_elev),
													 AIC(stdz.mod_spad_lai_isi_elev),
													 AIC(stdz.mod_spad_seed_isi_elev),
													 AIC(stdz.mod_spad_rand),
													 AIC(stdz.mod_spad_elev_ri),
													 AIC(stdz.mod_spad_isi_ri),
													 AIC(stdz.mod_spad_lai_ri),
													 AIC(stdz.mod_spad_seed_ri))

model_name_stdz.mod_spad_aic <- c("stdz.mod_spad_full",
																	"stdz.mod_spad_lai_seed_isi",
																	"stdz.mod_spad_lai_seed_elev",
																	"stdz.mod_spad_lai_isi_elev",
																	"stdz.mod_spad_seed_isi_elev",
																	"stdz.mod_spad_rand",
																	"stdz.mod_spad_elev_ri",
																	"stdz.mod_spad_isi_ri",
																	"stdz.mod_spad_lai_ri",
																	"stdz.mod_spad_seed_ri")

stdz.mod_spad_aic <- data.frame(model_name_stdz.mod_spad_aic, stdz.mod_spad_aic_val)
stdz.mod_spad_aic_arr <- arrange(stdz.mod_spad_aic, stdz.mod_spad_aic_val)
stdz.mod_spad_aic_arr
stdz.mod_spad_aic_arr$dAIC <- c(AIC(stdz.mod_spad_rand) - AIC(stdz.mod_spad_rand), 
																AIC(stdz.mod_spad_seed_ri) - AIC(stdz.mod_spad_rand), 
																AIC(stdz.mod_spad_lai_ri) - AIC(stdz.mod_spad_rand),
																AIC(stdz.mod_spad_elev_ri) - AIC(stdz.mod_spad_rand),
																AIC(stdz.mod_spad_isi_ri) - AIC(stdz.mod_spad_rand),
																AIC(stdz.mod_spad_lai_seed_elev) - AIC(stdz.mod_spad_rand),
																AIC(stdz.mod_spad_lai_seed_isi) - AIC(stdz.mod_spad_rand),
																AIC(stdz.mod_spad_seed_isi_elev) - AIC(stdz.mod_spad_rand),
																AIC(stdz.mod_spad_full) - AIC(stdz.mod_spad_rand),
																AIC(stdz.mod_spad_lai_isi_elev) - AIC(stdz.mod_spad_rand))


stdz.mod_spad_aic_arr
stdz.mod_spad_wi <- akaike.weights(stdz.mod_spad_aic_arr$stdz.mod_spad_aic_val)
stdz.mod_spad_wi
stdz.mod_spad_aic_arr$wi <- stdz.mod_spad_wi$weights
stdz.mod_spad_aic_arr <- stdz.mod_spad_aic_arr %>% dplyr::select(model_name_stdz.mod_spad_aic, stdz.mod_spad_aic_val, dAIC, wi)
r1 <- unname(r.squaredGLMM(stdz.mod_spad_rand))
r2 <-   unname(r.squaredGLMM(stdz.mod_spad_seed_ri))
r3 <-   unname(r.squaredGLMM(stdz.mod_spad_lai_ri))
r4 <-   unname(r.squaredGLMM(stdz.mod_spad_elev_ri))
r5 <-   unname(r.squaredGLMM(stdz.mod_spad_isi_ri))
r6 <-   unname(r.squaredGLMM(stdz.mod_spad_lai_seed_elev))
r7 <-   unname(r.squaredGLMM(stdz.mod_spad_lai_seed_isi))
r8 <-   unname(r.squaredGLMM(stdz.mod_spad_seed_isi_elev))
r9 <-   unname(r.squaredGLMM(stdz.mod_spad_full))
r10 <-   unname(r.squaredGLMM(stdz.mod_spad_lai_isi_elev))

stdz.mod_spad_aic_arr$R2c <- c(r1[2],r2[2], r3[2], r4[2], r5[2], r6[2], r7[2], r8[2], r9[2], r10[2])
stdz.mod_spad_aic_arr$r2m <- c(r1[1], r2[1], r3[1], r4[1], r5[1], r6[1], r7[1], r8[1], r9[1], r10[1])

stdz.mod_spad_aic_arr
stargazer(stdz.mod_spad_aic_arr, summary = F)

summary(stdz.mod_spad_lai_ri)
pamer.fnc(stdz.mod_spad_lai_seed_isi)
#Creating fixed effect slope graph with standard error bars
stdz.mod_spad_fix_name <- c("LAI", "ISI", "Elevation", "Comp.seed")                      
stdz.mod_spad_fix_slope <- c(-0.06607, NA,  NA, NA)
stdz.mod_spad_fix_se <- c(0.06725,NA , NA, NA)
stdz.mod_spad_fix_trait <- rep("spad", 4)
stdz.mod_spad_fix <- data.frame("factor" = stdz.mod_spad_fix_name, "slope" = stdz.mod_spad_fix_slope, "se" = stdz.mod_spad_fix_se, "trait" = stdz.mod_spad_fix_trait)
stdz.mod_spad_fix
stdz.mod_spad_sort <- c("Elevation", "ISI", "LAI", "Comp.seed")
stdz.mod_spad_fix <- stdz.mod_spad_fix %>% slice(match(stdz.mod_spad_sort, stdz.mod_spad_fix_name))
stdz.mod_spad_fix$stdz.mod_spad_fix_name <- factor(stdz.mod_spad_fix$stdz.mod_spad_fix_name, levels = stdz.mod_spad_fix$stdz.mod_spad_fix_name)
stdz.mod_spad_fix

#Height.leaf.ratio
stdz.mod_hlratio_full <- standardize(lmer(Height.leaf.ratio ~ LAI.4.ring + Comp.seed.total + Comp.adult.log.metric + Elevation + (1|Species) + (1|Site) , data=seedlings_rem_na,  REML = F), standardize.y = T)

stdz.mod_hlratio_int <- standardize(lmer(Height.leaf.ratio ~ LAI.4.ring + Comp.seed.total + Comp.adult.log.metric + Elevation + (Elevation + Comp.seed.total|Species) + (1|Site) , data=seedlings_rem_na,  REML = F), standardize.y = T)

stdz.mod_hlratio_lai_seed_isi <- standardize(lmer(Height.leaf.ratio ~ LAI.4.ring + Comp.seed.total + Comp.adult.log.metric + (1|Species) + (1|Site) , data=seedlings_rem_na,  REML = F), standardize.y = T)

stdz.mod_hlratio_lai_seed_elev <- standardize(lmer(Height.leaf.ratio ~ LAI.4.ring + Comp.seed.total + Elevation + (1|Species) + (1|Site) , data=seedlings_rem_na,  REML = F), standardize.y = T)

stdz.mod_hlratio_lai_isi_elev <- standardize(lmer(Height.leaf.ratio ~ LAI.4.ring + Comp.adult.log.metric + Elevation + (1|Species) + (1|Site), data=seedlings_rem_na,  REML = F), standardize.y = T)

stdz.mod_hlratio_seed_isi_elev <- standardize(lmer(Height.leaf.ratio ~ Comp.seed.total + Comp.adult.log.metric + Elevation + (1|Species) + (1|Site), data=seedlings_rem_na,  REML = F), standardize.y = T)

stdz.mod_hlratio_rand <- standardize(lmer(Height.leaf.ratio ~ (1|Species) + (1|Site), data=seedlings_rem_na, REML = F), standardize.y = T)

stdz.mod_hlratio_null <- lm(Height.leaf.ratio ~ 1, data=seedlings_rem_na)

stdz.mod_hlratio_aic_val <- c(AIC(stdz.mod_hlratio_full),
															AIC(stdz.mod_hlratio_int),
															AIC(stdz.mod_hlratio_lai_seed_isi),
															AIC(stdz.mod_hlratio_lai_seed_elev),
															AIC(stdz.mod_hlratio_lai_isi_elev),
															AIC(stdz.mod_hlratio_seed_isi_elev),
															AIC(stdz.mod_hlratio_rand),
															AIC(stdz.mod_hlratio_elev_rs),
															AIC(stdz.mod_hlratio_isi_ri),
															AIC(stdz.mod_hlratio_lai_ri),
															AIC(stdz.mod_hlratio_seed_ri))

model_name_stdz.mod_hlratio_aic <- c("stdz.mod_hlratio_full",
																		 "stdz.mod_hlratio_int",
																		 "stdz.mod_hlratio_lai_seed_isi",
																		 "stdz.mod_hlratio_lai_seed_elev",
																		 "stdz.mod_hlratio_lai_isi_elev",
																		 "stdz.mod_hlratio_seed_isi_elev",
																		 "stdz.mod_hlratio_rand",
																		 "stdz.mod_hlratio_elev_rs",
																		 "stdz.mod_hlratio_isi_ri",
																		 "stdz.mod_hlratio_lai_ri",
																		 "stdz.mod_hlratio_seed_ri")

stdz.mod_hlratio_aic <- data.frame(model_name_stdz.mod_hlratio_aic, stdz.mod_hlratio_aic_val)
stdz.mod_hlratio_aic_arr <- arrange(stdz.mod_hlratio_aic, stdz.mod_hlratio_aic_val)
stdz.mod_hlratio_aic_arr
stdz.mod_hlratio_aic_arr$dAIC <- c(AIC(stdz.mod_hlratio_int) - AIC(stdz.mod_hlratio_rand),
																	 AIC(stdz.mod_hlratio_elev_rs) - AIC(stdz.mod_hlratio_rand),
																	 AIC(stdz.mod_hlratio_lai_seed_elev) - AIC(stdz.mod_hlratio_rand),
																	 AIC(stdz.mod_hlratio_lai_isi_elev) - AIC(stdz.mod_hlratio_rand),
																	 AIC(stdz.mod_hlratio_lai_ri) - AIC(stdz.mod_hlratio_rand),
																	 AIC(stdz.mod_hlratio_full) - AIC(stdz.mod_hlratio_rand),
																	 AIC(stdz.mod_hlratio_rand) - AIC(stdz.mod_hlratio_rand),
																	 AIC(stdz.mod_hlratio_seed_isi_elev) - AIC(stdz.mod_hlratio_rand),
																	 AIC(stdz.mod_hlratio_seed_ri) - AIC(stdz.mod_hlratio_rand),
																	 AIC(stdz.mod_hlratio_isi_ri) - AIC(stdz.mod_hlratio_rand),
																	 AIC(stdz.mod_hlratio_lai_seed_isi) - AIC(stdz.mod_hlratio_rand))

stdz.mod_hlratio_aic_arr
stdz.mod_hlratio_wi <- Weights(stdz.mod_hlratio_aic_arr$stdz.mod_hlratio_aic_val)
stdz.mod_hlratio_wi
stdz.mod_hlratio_aic_arr$wi <- stdz.mod_hlratio_wi
stdz.mod_hlratio_aic_arr <- stdz.mod_hlratio_aic_arr %>% dplyr::select(model_name_stdz.mod_hlratio_aic, stdz.mod_hlratio_aic_val, dAIC, wi)
r1 <- unname(r.squaredGLMM(stdz.mod_hlratio_int))
r2 <-   unname(r.squaredGLMM(stdz.mod_hlratio_elev_rs))
r3 <-   unname(r.squaredGLMM(stdz.mod_hlratio_lai_seed_elev))
r4 <-   unname(r.squaredGLMM(stdz.mod_hlratio_lai_isi_elev))
r5 <-   unname(r.squaredGLMM(stdz.mod_hlratio_lai_ri))
r6 <-   unname(r.squaredGLMM(stdz.mod_hlratio_full))
r7 <-   unname(r.squaredGLMM(stdz.mod_hlratio_rand))
r8 <- unname(r.squaredGLMM(stdz.mod_hlratio_seed_isi_elev))
r9 <- unname(r.squaredGLMM(stdz.mod_hlratio_seed_ri))
r10 <- unname(r.squaredGLMM(stdz.mod_hlratio_isi_ri))
r11 <- unname(r.squaredGLMM(stdz.mod_hlratio_lai_seed_isi))

stdz.mod_hlratio_aic_arr$R2c <- c(r1[2],r2[2], r3[2], r4[2], r5[2], r6[2], r7[2], r7[2], r9[2], r10[2], r11[2])
stdz.mod_hlratio_aic_arr$r2m <- c(r1[1], r2[1], r3[1], r4[1], r5[1], r6[1], r7[1], r8[1], r9[1], r10[1], r11[1])

stdz.mod_hlratio_aic_arr
stargazer(stdz.mod_hlratio_aic_arr, summary = F)
summary(stdz.mod_hlratio_seed_isi_elev)


summary(stdz.mod_hlratio_int)
#Creating fixed effect slope graph with standard error bars
stdz.mod_hlratio_fix_name <- c("LAI", "ISI", "Elevation", "Comp.seed")                      
stdz.mod_hlratio_fix_slope <- c(-0.10258, 0.01219, -0.14145 , 0.02634)
stdz.mod_hlratio_fix_se <- c(0.05014, 0.04830, 0.12770, 0.08709)
stdz.mod_hlratio_fix_trait <- rep("hlratio", 4)
stdz.mod_hlratio_fix <- data.frame("factor" = stdz.mod_hlratio_fix_name, "slope" = stdz.mod_hlratio_fix_slope, "se" = stdz.mod_hlratio_fix_se, "trait" = stdz.mod_hlratio_fix_trait)
stdz.mod_hlratio_fix
stdz.mod_hlratio_sort <- c("Elevation", "ISI", "LAI", "Comp.seed")
stdz.mod_hlratio_fix <- stdz.mod_hlratio_fix %>% slice(match(stdz.mod_hlratio_sort, stdz.mod_hlratio_fix_name))
stdz.mod_hlratio_fix$stdz.mod_hlratio_fix_name <- factor(stdz.mod_hlratio_fix$stdz.mod_hlratio_fix_name, levels = stdz.mod_hlratio_fix$stdz.mod_hlratio_fix_name)
stdz.mod_hlratio_fix

#Leaf.area
stdz.mod_area_full <- standardize(lmer(Leaf.area ~ LAI.4.ring + Comp.seed.total + Comp.adult.log.metric + Elevation + (1|Species) + (1|Site) , data=seedlings_rem_na,  REML = F), standardize.y = T)

stdz.mod_area_int <- standardize(lmer(Leaf.area ~ LAI.4.ring + Comp.seed.total + Comp.adult.log.metric + Elevation + (1 + LAI.4.ring + Comp.seed.total + Elevation|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)

stdz.mod_area_lai_seed_isi <- standardize(lmer(Leaf.area ~ LAI.4.ring + Comp.seed.total + Comp.adult.log.metric + (1|Species) + (1|Site) , data=seedlings_rem_na,  REML = F), standardize.y = T)

stdz.mod_area_lai_seed_elev <- standardize(lmer(Leaf.area ~ LAI.4.ring + Comp.seed.total + Elevation + (1|Species) + (1|Site) , data=seedlings_rem_na,  REML = F), standardize.y = T)

stdz.mod_area_lai_isi_elev <- standardize(lmer(Leaf.area ~ LAI.4.ring + Comp.adult.log.metric + Elevation + (1|Species) + (1|Site), data=seedlings_rem_na,  REML = F), standardize.y = T)

stdz.mod_area_seed_isi_elev <- standardize(lmer(Leaf.area ~ Comp.seed.total + Comp.adult.log.metric + Elevation + (1|Species) + (1|Site), data=seedlings_rem_na,  REML = F), standardize.y = T)

stdz.mod_area_rand <- standardize(lmer(Leaf.area ~ (1|Species) + (1|Site), data=seedlings_rem_na, REML = F),standardize.y = T)

stdz.mod_area_null <- lm(Leaf.area ~ 1, data=seedlings_rem_na)

stdz.mod_area_aic_val <- c(AIC(stdz.mod_area_full),
													 AIC(stdz.mod_area_int),
													 AIC(stdz.mod_area_lai_seed_isi),
													 AIC(stdz.mod_area_lai_seed_elev),
													 AIC(stdz.mod_area_lai_isi_elev),
													 AIC(stdz.mod_area_seed_isi_elev),
													 AIC(stdz.mod_area_rand),
													 AIC(stdz.mod_area_elev_rs),
													 AIC(stdz.mod_area_isi_ri),
													 AIC(stdz.mod_area_lai_rs),
													 AIC(stdz.mod_area_seed_rs))


model_name_stdz.mod_area_aic <- c("stdz.mod_area_full",
	"stdz.mod_area_int",
	"stdz.mod_area_lai_seed_isi",
	"stdz.mod_area_lai_seed_elev",
	"stdz.mod_area_lai_isi_elev",
	"stdz.mod_area_seed_isi_elev",
	"stdz.mod_area_rand",
	"stdz.mod_area_elev_rs",
	"stdz.mod_area_isi_ri",
	"stdz.mod_area_lai_rs",
	"stdz.mod_area_seed_rs")

stdz.mod_area_aic <- data.frame(model_name_stdz.mod_area_aic, stdz.mod_area_aic_val)
stdz.mod_area_aic_arr <- arrange(stdz.mod_area_aic, stdz.mod_area_aic_val)
stdz.mod_area_aic_arr
stdz.mod_area_aic_arr$dAIC <- c(AIC(stdz.mod_area_lai_rs) - AIC(stdz.mod_area_rand), 
AIC(stdz.mod_area_int) - AIC(stdz.mod_area_rand), 
AIC(stdz.mod_area_elev_rs) - AIC(stdz.mod_area_rand),
AIC(stdz.mod_area_rand) - AIC(stdz.mod_area_rand),
AIC(stdz.mod_area_seed_rs) - AIC(stdz.mod_area_rand),
AIC(stdz.mod_area_isi_ri) - AIC(stdz.mod_area_rand),
AIC(stdz.mod_area_lai_isi_elev) - AIC(stdz.mod_area_rand),
AIC(stdz.mod_area_lai_seed_isi) - AIC(stdz.mod_area_rand),
AIC(stdz.mod_area_lai_seed_elev) - AIC(stdz.mod_area_rand),
AIC(stdz.mod_area_seed_isi_elev) - AIC(stdz.mod_area_rand),
AIC(stdz.mod_area_full) - AIC(stdz.mod_area_rand))

stdz.mod_area_aic_arr
stdz.mod_area_wi <- akaike.weights(stdz.mod_area_aic_arr$stdz.mod_area_aic_val)
stdz.mod_area_wi
stdz.mod_area_aic_arr$wi <- stdz.mod_area_wi$weights
stdz.mod_area_aic_arr <- stdz.mod_area_aic_arr %>% dplyr::select(model_name_stdz.mod_area_aic, stdz.mod_area_aic_val, dAIC, wi)
r1 <- unname(r.squaredGLMM(stdz.mod_area_lai_rs))
r2 <-   unname(r.squaredGLMM(stdz.mod_area_int))
r3 <-   unname(r.squaredGLMM(stdz.mod_area_elev_rs))
r4 <-   unname(r.squaredGLMM(stdz.mod_area_rand))
r5 <-   unname(r.squaredGLMM(stdz.mod_area_seed_rs))
r6 <-   unname(r.squaredGLMM(stdz.mod_area_isi_ri))
r7 <-   unname(r.squaredGLMM(stdz.mod_area_lai_isi_elev))
r8 <-   unname(r.squaredGLMM(stdz.mod_area_lai_seed_isi))
r9 <-   unname(r.squaredGLMM(stdz.mod_area_lai_seed_elev))
r10 <-   unname(r.squaredGLMM(stdz.mod_area_seed_isi_elev))
r11 <-   unname(r.squaredGLMM(stdz.mod_area_full))

stdz.mod_area_aic_arr$R2c <- c(r1[2],r2[2], r3[2], r4[2], r5[2], r6[2], r7[2], r8[2], r9[2], r10[2], r11[2])
stdz.mod_area_aic_arr$r2m <- c(r1[1], r2[1], r3[1], r4[1], r5[1], r6[1], r7[1], r8[1], r9[1], r10[1], r11[1])

stdz.mod_area_aic_arr
stargazer(stdz.mod_area_aic_arr, summary = F)

summary(stdz.mod_area_int)
pamer.fnc(stdz.mod_area_int)
#Creating fixed effect slope graph with standard error bars
stdz.mod_area_fix_name <- c("LAI", "ISI", "Elevation", "Comp.seed")                      
stdz.mod_area_fix_slope <- c(-0.06569, -0.01211, -0.37745, 0.10623)
stdz.mod_area_fix_se <- c(0.08959, 0.05213, 0.35761, 0.10254)
stdz.mod_area_fix_trait <- rep("area", 4)
stdz.mod_area_fix <- data.frame("factor" = stdz.mod_area_fix_name, "slope" = stdz.mod_area_fix_slope, "se" = stdz.mod_area_fix_se, "trait" = stdz.mod_area_fix_trait)
stdz.mod_area_fix
stdz.mod_area_sort <- c("Elevation", "ISI", "LAI", "Comp.seed")
stdz.mod_area_fix <- stdz.mod_area_fix %>% slice(match(stdz.mod_area_sort, stdz.mod_area_fix_name))
stdz.mod_area_fix$stdz.mod_area_fix_name <- factor(stdz.mod_area_fix$stdz.mod_area_fix_name, levels = stdz.mod_area_fix$stdz.mod_area_fix_name)

#Stem Volume
stdz.mod_stemvol_full <- standardize(lmer(Stem.volume.cm3 ~ LAI.4.ring + Comp.seed.total + Comp.adult.log.metric + Elevation + (1|Species) + (1|Site) , data=seedlings_rem_na,  REML = F), standardize.y = T)

stdz.mod_stemvol_int <- standardize(lmer(Stem.volume.cm3 ~ LAI.4.ring + Comp.seed.total + Comp.adult.log.metric + Elevation + (1 + Comp.adult.log.metric|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)

stdz.mod_stemvol_lai_seed_isi <- standardize(lmer(Stem.volume.cm3 ~ LAI.4.ring + Comp.seed.total + Comp.adult.log.metric + (1|Species) + (1|Site) , data=seedlings_rem_na,  REML = F), standardize.y = T)

stdz.mod_stemvol_lai_seed_elev <- standardize(lmer(Stem.volume.cm3 ~ LAI.4.ring + Comp.seed.total + Elevation + (1|Species) + (1|Site) , data=seedlings_rem_na,  REML = F), standardize.y = T)

stdz.mod_stemvol_lai_isi_elev <- standardize(lmer(Stem.volume.cm3 ~ LAI.4.ring + Comp.adult.log.metric + Elevation + (1|Species) + (1|Site), data=seedlings_rem_na,  REML = F), standardize.y = T)

stdz.mod_stemvol_seed_isi_elev <- standardize(lmer(Stem.volume.cm3 ~ Comp.seed.total + Comp.adult.log.metric + Elevation + (1|Species) + (1|Site), data=seedlings_rem_na,  REML = F), standardize.y = T)

stdz.mod_stemvol_rand <- standardize(lmer(Stem.volume.cm3 ~ (1|Species) + (1|Site), data=seedlings_rem_na, REML = F),standardize.y = T)

stdz.mod_stemvol_null <- lm(Stem.volume.cm3 ~ 1, data=seedlings_rem_na)

stdz.mod_stemvol_aic_val <- c(AIC(stdz.mod_stemvol_full),
															AIC(stdz.mod_stemvol_int),
															AIC(stdz.mod_stemvol_lai_seed_isi),
															AIC(stdz.mod_stemvol_lai_seed_elev),
															AIC(stdz.mod_stemvol_lai_isi_elev),
															AIC(stdz.mod_stemvol_seed_isi_elev),
															AIC(stdz.mod_stemvol_rand),
															AIC(stdz.mod_stemvol_elev_ri),
															AIC(stdz.mod_stemvol_isi_rs),
															AIC(stdz.mod_stemvol_lai_ri),
															AIC(stdz.mod_stemvol_seed_ri))

model_name_stdz.mod_stemvol_aic <- c("stdz.mod_stemvol_full",
		 "stdz.mod_stemvol_int",
		 "stdz.mod_stemvol_lai_seed_isi",
		 "stdz.mod_stemvol_lai_seed_elev",
		 "stdz.mod_stemvol_lai_isi_elev",
		 "stdz.mod_stemvol_seed_isi_elev",
		 "stdz.mod_stemvol_rand",
		 "stdz.mod_stemvol_elev_ri",
		 "stdz.mod_stemvol_isi_rs",
		 "stdz.mod_stemvol_lai_ri",
		 "stdz.mod_stemvol_seed_ri")

stdz.mod_stemvol_aic <- data.frame(model_name_stdz.mod_stemvol_aic, stdz.mod_stemvol_aic_val)
stdz.mod_stemvol_aic_arr <- arrange(stdz.mod_stemvol_aic, stdz.mod_stemvol_aic_val)
stdz.mod_stemvol_aic_arr
stdz.mod_stemvol_aic_arr$dAIC <- c(AIC(stdz.mod_stemvol_int) - AIC(stdz.mod_stemvol_rand), 
	 AIC(stdz.mod_stemvol_isi_rs) - AIC(stdz.mod_stemvol_rand), 
	 AIC(stdz.mod_stemvol_seed_isi_elev) - AIC(stdz.mod_stemvol_rand),
	 AIC(stdz.mod_stemvol_full) - AIC(stdz.mod_stemvol_rand),
	 AIC(stdz.mod_stemvol_lai_isi_elev) - AIC(stdz.mod_stemvol_rand),
	 AIC(stdz.mod_stemvol_lai_seed_isi) - AIC(stdz.mod_stemvol_rand),
	 AIC(stdz.mod_stemvol_elev_ri) - AIC(stdz.mod_stemvol_rand),
	 AIC(stdz.mod_stemvol_lai_seed_elev) - AIC(stdz.mod_stemvol_rand),
	 AIC(stdz.mod_stemvol_seed_ri) - AIC(stdz.mod_stemvol_rand),
	 AIC(stdz.mod_stemvol_rand) - AIC(stdz.mod_stemvol_rand),
	 AIC(stdz.mod_stemvol_lai_ri) - AIC(stdz.mod_stemvol_rand))


stdz.mod_stemvol_aic_arr
stdz.mod_stemvol_wi <- akaike.weights(stdz.mod_stemvol_aic_arr$stdz.mod_stemvol_aic_val)
stdz.mod_stemvol_wi
stdz.mod_stemvol_aic_arr$wi <- stdz.mod_stemvol_wi$weights
stdz.mod_stemvol_aic_arr <- stdz.mod_stemvol_aic_arr %>% dplyr::select(model_name_stdz.mod_stemvol_aic, stdz.mod_stemvol_aic_val, dAIC, wi)
r1 <- unname(r.squaredGLMM(stdz.mod_stemvol_int))
r2 <-   unname(r.squaredGLMM(stdz.mod_stemvol_isi_rs))
r3 <-   unname(r.squaredGLMM(stdz.mod_stemvol_seed_isi_elev))
r4 <-   unname(r.squaredGLMM(stdz.mod_stemvol_full))
r5 <-   unname(r.squaredGLMM(stdz.mod_stemvol_lai_isi_elev))
r6 <-   unname(r.squaredGLMM(stdz.mod_stemvol_lai_seed_isi))
r7 <-   unname(r.squaredGLMM(stdz.mod_stemvol_elev_ri))
r8 <-   unname(r.squaredGLMM(stdz.mod_stemvol_lai_seed_elev))
r9 <-   unname(r.squaredGLMM(stdz.mod_stemvol_seed_ri))
r10 <-   unname(r.squaredGLMM(stdz.mod_stemvol_rand))
r11 <-   unname(r.squaredGLMM(stdz.mod_stemvol_lai_ri))

stdz.mod_stemvol_aic_arr$R2c <- c(r1[2],r2[2], r3[2], r4[2], r5[2], r6[2], r7[2], r8[2], r9[2], r10[2], r11[2])
stdz.mod_stemvol_aic_arr$r2m <- c(r1[1], r2[1], r3[1], r4[1], r5[1], r6[1], r7[1], r8[1], r9[1], r10[1], r11[1])

stdz.mod_stemvol_aic_arr
stargazer(stdz.mod_stemvol_aic_arr, summary = F)

summary(stdz.mod_stemvol_int)
pamer.fnc(stdz.mod_stemvol_int)
#Creating fixed effect slope graph with standard error bars
stdz.mod_stemvol_fix_name <- c("LAI", "ISI", "Elevation", "Comp.seed")                      
stdz.mod_stemvol_fix_slope <- c(-0.05366, -0.19821, 0.31046, -0.06353)
stdz.mod_stemvol_fix_se <- c(0.05709, 0.11863, 0.13184, 0.05322)
stdz.mod_stemvol_fix_trait <- rep("stemvol", 4)
stdz.mod_stemvol_fix <- data.frame("factor" = stdz.mod_stemvol_fix_name, "slope" = stdz.mod_stemvol_fix_slope, "se" = stdz.mod_stemvol_fix_se, "trait" = stdz.mod_stemvol_fix_trait)
stdz.mod_stemvol_fix
stdz.mod_stemvol_sort <- c("Elevation", "ISI", "LAI", "Comp.seed")
stdz.mod_stemvol_fix <- stdz.mod_stemvol_fix %>% slice(match(stdz.mod_stemvol_sort, stdz.mod_stemvol_fix_name))
stdz.mod_stemvol_fix$stdz.mod_stemvol_fix_name <- factor(stdz.mod_stemvol_fix$stdz.mod_stemvol_fix_name, levels = stdz.mod_stemvol_fix$stdz.mod_stemvol_fix_name)
stdz.mod_stemvol_fix

#Leaf Thickness
stdz.mod_thick_full <- standardize(lmer(Leaf.thick.mean.mm ~ LAI.4.ring + Comp.seed.total + Comp.adult.log.metric + Elevation + (1|Species) + (1|Site) , data=seedlings_rem_na,  REML = F), standardize.y = T)

stdz.mod_thick_int <- standardize(lmer(Leaf.thick.mean.mm ~ LAI.4.ring + Comp.seed.total + Comp.adult.log.metric + Elevation + (1 + Comp.seed.total|Species) + (1|Site), data = seedlings_rem_na, REML = F), standardize.y = T)

stdz.mod_thick_lai_seed_isi <- standardize(lmer(Leaf.thick.mean.mm ~ LAI.4.ring + Comp.seed.total + Comp.adult.log.metric + (1|Species) + (1|Site) , data=seedlings_rem_na,  REML = F), standardize.y = T)

stdz.mod_thick_lai_seed_elev <- standardize(lmer(Leaf.thick.mean.mm ~ LAI.4.ring + Comp.seed.total + Elevation + (1|Species) + (1|Site) , data=seedlings_rem_na,  REML = F), standardize.y = T)

stdz.mod_thick_lai_isi_elev <- standardize(lmer(Leaf.thick.mean.mm ~ LAI.4.ring + Comp.adult.log.metric + Elevation + (1|Species) + (1|Site), data=seedlings_rem_na,  REML = F), standardize.y = T)

stdz.mod_thick_seed_isi_elev <- standardize(lmer(Leaf.thick.mean.mm ~ Comp.seed.total + Comp.adult.log.metric + Elevation + (1|Species) + (1|Site), data=seedlings_rem_na,  REML = F), standardize.y = T)

stdz.mod_thick_rand <- standardize(lmer(Leaf.thick.mean.mm ~ (1|Species) + (1|Site), data=seedlings_rem_na, REML = F),standardize.y = T)

stdz.mod_thick_null <- lm(Leaf.thick.mean.mm ~ 1, data=seedlings_rem_na)

stdz.mod_thick_aic_val <- c(AIC(stdz.mod_thick_full),
														AIC(stdz.mod_thick_int),
														AIC(stdz.mod_thick_lai_seed_isi),
														AIC(stdz.mod_thick_lai_seed_elev),
														AIC(stdz.mod_thick_lai_isi_elev),
														AIC(stdz.mod_thick_seed_isi_elev),
														AIC(stdz.mod_thick_rand),
														AIC(stdz.mod_thick_elev_ri),
														AIC(stdz.mod_thick_isi_ri),
														AIC(stdz.mod_thick_lai_ri),
														AIC(stdz.mod_thick_seed_rs))

model_name_stdz.mod_thick_aic <- c("stdz.mod_thick_full",
	 "stdz.mod_thick_int",
	 "stdz.mod_thick_lai_seed_isi",
	 "stdz.mod_thick_lai_seed_elev",
	 "stdz.mod_thick_lai_isi_elev",
	 "stdz.mod_thick_seed_isi_elev",
	 "stdz.mod_thick_rand",
	 "stdz.mod_thick_elev_ri",
	 "stdz.mod_thick_isi_ri",
	 "stdz.mod_thick_lai_ri",
	 "stdz.mod_thick_seed_rs")

stdz.mod_thick_aic <- data.frame(model_name_stdz.mod_thick_aic, stdz.mod_thick_aic_val)
stdz.mod_thick_aic_arr <- arrange(stdz.mod_thick_aic, stdz.mod_thick_aic_val)
stdz.mod_thick_aic_arr
stdz.mod_thick_aic_arr$dAIC <- c(AIC(stdz.mod_thick_lai_seed_elev) - AIC(stdz.mod_thick_rand), 
 AIC(stdz.mod_thick_lai_isi_elev) - AIC(stdz.mod_thick_rand), 
 AIC(stdz.mod_thick_full) - AIC(stdz.mod_thick_rand),
 AIC(stdz.mod_thick_lai_ri) - AIC(stdz.mod_thick_rand),
 AIC(stdz.mod_thick_int) - AIC(stdz.mod_thick_rand),
 AIC(stdz.mod_thick_elev_ri) - AIC(stdz.mod_thick_rand),
 AIC(stdz.mod_thick_lai_seed_isi) - AIC(stdz.mod_thick_rand),
 AIC(stdz.mod_thick_seed_isi_elev) - AIC(stdz.mod_thick_rand),
 AIC(stdz.mod_thick_rand) - AIC(stdz.mod_thick_rand),
 AIC(stdz.mod_thick_isi_ri) - AIC(stdz.mod_thick_rand),
 AIC(stdz.mod_thick_seed_rs) - AIC(stdz.mod_thick_rand))

stdz.mod_thick_aic_arr
stdz.mod_thick_wi <- akaike.weights(stdz.mod_thick_aic_arr$stdz.mod_thick_aic_val)
stdz.mod_thick_wi
stdz.mod_thick_aic_arr$wi <- stdz.mod_thick_wi$weights
stdz.mod_thick_aic_arr <- stdz.mod_thick_aic_arr %>% dplyr::select(model_name_stdz.mod_thick_aic, stdz.mod_thick_aic_val, dAIC, wi)
r1 <- unname(r.squaredGLMM(stdz.mod_thick_lai_seed_elev))
r2 <-   unname(r.squaredGLMM(stdz.mod_thick_lai_isi_elev))
r3 <-   unname(r.squaredGLMM(stdz.mod_thick_full))
r4 <-   unname(r.squaredGLMM(stdz.mod_thick_lai_ri))
r5 <-   unname(r.squaredGLMM(stdz.mod_thick_int))
r6 <-   unname(r.squaredGLMM(stdz.mod_thick_elev_ri))
r7 <-   unname(r.squaredGLMM(stdz.mod_thick_lai_seed_isi))
r8 <-   unname(r.squaredGLMM(stdz.mod_thick_seed_isi_elev))
r9 <-   unname(r.squaredGLMM(stdz.mod_thick_rand))
r10 <-   unname(r.squaredGLMM(stdz.mod_thick_isi_ri))
r11 <-   unname(r.squaredGLMM(stdz.mod_thick_seed_rs))

stdz.mod_thick_aic_arr$R2c <- c(r1[2],r2[2], r3[2], r4[2], r5[2], r6[2], r7[2], r8[2], r9[2], r10[2], r11[2])
stdz.mod_thick_aic_arr$r2m <- c(r1[1], r2[1], r3[1], r4[1], r5[1], r6[1], r7[1], r8[1], r9[1], r10[1], r11[1])

stdz.mod_thick_aic_arr
stargazer(stdz.mod_thick_aic_arr, summary = F)

summary(stdz.mod_thick_lai_isi_elev)
pamer.fnc(stdz.mod_thick_lai_isi_elev)
#Creating fixed effect slope graph with standard error bars
stdz.mod_thick_fix_name <- c("LAI", "ISI", "Elevation", "Comp.seed")                      
stdz.mod_thick_fix_slope <- c(-0.13643, -0.02387, 0.29271, NA)
stdz.mod_thick_fix_se <- c(0.04091, 0.03878, 0.08471, NA)
stdz.mod_thick_fix_trait <- rep("thick", 4)
stdz.mod_thick_fix <- data.frame("factor" = stdz.mod_thick_fix_name, "slope" = stdz.mod_thick_fix_slope, "se" = stdz.mod_thick_fix_se, "trait" = stdz.mod_thick_fix_trait)
stdz.mod_thick_sort <- c("Elevation", "ISI", "LAI", "Comp.seed")
stdz.mod_thick_fix <- stdz.mod_thick_fix %>% slice(match(stdz.mod_thick_sort, stdz.mod_thick_fix_name))
stdz.mod_thick_fix$stdz.mod_thick_fix_name <- factor(stdz.mod_thick_fix$stdz.mod_thick_fix_name, levels = stdz.mod_thick_fix$stdz.mod_thick_fix_name)
stdz.mod_thick_fix
stdz.mod_thick_fix_plot1 <- ggplot() + geom_errorbar(data = stdz.mod_thick_fix, aes(x = stdz.mod_thick_fix_name, ymin = stdz.mod_thick_fix_slope - stdz.mod_thick_fix_se, ymax = stdz.mod_thick_fix_slope + stdz.mod_thick_fix_se, colour = stdz.mod_thick_fix_name)) 
stdz.mod_thick_fix_plot2 <- stdz.mod_thick_fix_plot1 + geom_point(data=stdz.mod_thick_fix, aes(x = stdz.mod_thick_fix_name, y = stdz.mod_thick_fix_slope, colour = stdz.mod_thick_fix_name), size = 5)
stdz.mod_thick_fix_plot3 <- stdz.mod_thick_fix_plot2 + theme_bw() +theme(plot.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank() )+theme(panel.border= element_blank())+theme(axis.line.x = element_line(color="black", size = 1),axis.line.y = element_line(color="black", size = 1)) + theme(legend.position = "right")+ xlab("Fixed Effects") + ylab(expression(paste("Fixed Effect Slope")))
stdz.mod_thick_fix_plot4 <- stdz.mod_thick_fix_plot3 + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14)) + theme(legend.position = "none") + ggtitle("Leaf Thickness") + theme(axis.title.x = element_blank())+ geom_hline(yintercept = 0, linetype = 5) + annotate("text", x = 4, y = 0, label = "NA", size = 10) + ylim(-0.75, 0.75) +  scale_x_discrete("stdz.mod_fvfm_fix_name", labels = c("Elevation" = "Elevation", "ISI" = "ISI", "LAI" = "LAI", "Comp.seed" = "Herbaceous Plants"))
stdz.mod_thick_fix_plot4

best_model_name <- c("D.FvFm", "SPAD", "Leaf Thickness", "Height:Leaf Ratio", "Leaf Area", "Stem Volume")
best_model_lai <- c("","","","","","")
best_model_seed <- c("","","","","","")
best_model_isi <-  c("","","","","","")
best_model_elev <- c("","","","","","")
best_model_dAIC <- c(-15.849, -45.263, -214.505, -73.863, NA, NA)
best_model_wi <- c(0.495, 0.524, 0.413, 0.439, NA, NA)
best_model_r2c <- c(0.320, 0.321, 0.758, 0.404, NA, NA)
best_model_r2m <- c(0.140, 0, 0.123, 0, NA, NA)

stargazer(stdz.mod_fvfm_aic_arr, summary = F) 
stargazer(stdz.mod_spad_aic_arr, summary = F) 
stargazer(stdz.mod_thick_aic_arr, summary = F)
stargazer(stdz.mod_hlratio_aic_arr, summary = F)
stargazer(stdz.mod_area_aic_arr, summary = F)
stargazer(stdz.mod_stemvol_aic_arr, summary = F)

stargazer(stdz.mod_fvfm_lai_isi_elev, stdz.mod_spad_lai_ri, stdz.mod_thick_lai_isi_elev, 
					stdz.mod_hlratio_int, stdz.mod_area_int, stdz.mod_stemvol_int)

AIC(stdz.mod_fvfm_rand) -  AIC(stdz.mod_fvfm_lai_isi_elev)
AIC(stdz.mod_spad_rand) - AIC(stdz.mod_spad_lai_ri)
AIC(stdz.mod_thick_rand) -AIC(stdz.mod_thick_lai_isi_elev)
AIC(stdz.mod_hlratio_rand) -AIC(stdz.mod_hlratio_int)
AIC(stdz.mod_area_rand) -AIC(stdz.mod_area_int)
AIC(stdz.mod_stemvol_rand) -AIC(stdz.mod_stemvol_int)

r.squaredGLMM(stdz.mod_fvfm_lai_isi_elev)
r.squaredGLMM(stdz.mod_spad_lai_ri)
r.squaredGLMM(stdz.mod_thick_lai_isi_elev)
r.squaredGLMM(stdz.mod_hlratio_int)
r.squaredGLMM(stdz.mod_area_int)
r.squaredGLMM(stdz.mod_stemvol_int)

# Combine all model effects into one data frame
stdz.mod_all_fix <- rbind(stdz.mod_fvfm_fix, stdz.mod_spad_fix, stdz.mod_area_fix,
													stdz.mod_thick_fix, stdz.mod_hlratio_fix, stdz.mod_stemvol_fix)

# Plot the effect sizes for multi-predictor models
ggplot(stdz.mod_all_fix) +
	geom_errorbar(aes(x = factor, ymin = slope - se, ymax = slope + se, colour = factor)) + 
	geom_point(aes(x = factor, y = slope, colour = factor), size = 5) + 
	facet_wrap(~trait, scales = "fixed")





