# Final analyses for elevation seedling paper
# John Godlee (johngodlee@gmail.com)
# 2019_08_12

# Clear Global Env.
rm(list = ls())

# Packages 
library(dplyr)  # %>% etc.
library(scales)  # rescale()
library(lme4)  # lmer()
library(tibble)  # rownames_to_column()
library(MuMIn)  # Weights()
library(ggplot2)  # ggplot()
library(tidyr)  # gather()
library(stargazer)  # stargazer()
library(stringr)  # str_extract()

# Set working directory
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("~/google_drive/postgrad/extra_projects/elevation_tree_comp/analysis/")

# Load and clean data

## Genus Level migration rates (Feeley et al. 2011)
genus_mig_rates <- read.csv("data/genus_mig_rates.csv")

## Data for rank abundance curve
aberg_census <- read.csv("data/aberg.csv")

## Species ranges
ranges <- read.csv("data/species_range.csv")

## Site locations
camp_loc <- read.csv("data/camp_loc.csv")

## Elevation of measurements
species_site_elev <- read.csv("data/species_site_elev.csv")

## Site environmental characteristics
site_char <- read.csv("data/site_char.csv")

## Seedling measurements
seedlings <- read.csv("data/seedlings.csv")

## Clean seedling measurements
seedlings_clean <- seedlings %>%
  filter(Comp.Y.N. == "Y") %>%  # Competition seedlings only
  filter(!is.na(LAI.4.ring)) %>%  # Remove NAs
  filter(!is.na(Comp.seed.total)) %>%
  filter(!is.na(Comp.adult.log.metric)) %>%
  filter(!is.na(Elevation)) %>%
  filter(!is.na(Species)) %>%
  mutate(LMA = Leaf.mass.dry.g / Leaf.area,  # Create Leaf Mass per area
    Leaf.height.ratio = No.leaves / Height.cm,  # Create Height to leaf ratio
    leaf_chl = 0.53 * exp(1)^(0.0364 * SPAD.mean),  # Create Chl-a
    LAI.4.ring_scale = rescale(LAI.4.ring),  # Rescale fixed effects
    Comp.seed.total_scale = rescale(Comp.seed.total),
    Comp.adult.log.metric_scale = rescale(Comp.adult.log.metric),
    Elevation_scale = rescale(Elevation),
    num_leaves = No.leaves + No.cot) %>%  # Create total number leaves
  dplyr::select(-Adult.Seedling,  # Remove columns unused in analysis
    -Soil.temp.1, -Soil.temp.2, -Soil.temp.3,
    -Soil.mois.1, -Soil.mois.2, -Soil.mois.3,
    -Height.mm,
    -No.leaves,
    -No.cot,
    -Stem.volume.mm3,
    -L.time, -L.temp, -L.light,
    -D.time, -D.temp, -D.light,
    -Canopy.open.per,
    -SPAD.1, -SPAD.2, -SPAD.3,
    -Leaf.thick.1.mm, -Leaf.thick.2.mm, -Leaf.thick.3.mm,
    -Comp.Y.N.) %>%
  rename("site" = "Site",
    "elev_code" = "Elevation.code",
    "species" = "Species", 
    "ind_code" = "Individual.code",
    "lat" = "Lat.DD",
    "lon" = "Long.DD",
    "coll_date" = "Collection.date",
    "coll_day" = "Collection.day",
    "elev" = "Elevation", 
    "soil_temp_mean" = "Soil.temp.mean",
    "soil_temp_sd" = "Soil.temp.SD",
    "soil_mois_mean" = "Soil.mois.mean",
    "soil_mois_sd" = "Soil.mois.SD",
    "height_cm" = "Height.cm",
    "width_mm" = "Width.mm",
    "stem_vol_cm3" = "Stem.volume.cm3",
    "stem_thick_per_height" = "Stem.thick.per.height",
    "qpsii" = "QPSII",
    "d_fvfm" = "D.FvFm",
    "lai" = "LAI.4.ring",
    "spad_mean" = "SPAD.mean",
    "spad_sd" = "SPAD.SD",
    "leaf_thick_mean_mm" = "Leaf.thick.mean.mm",
    "leaf_thick_sd" = "Leaf.thick.SD",
    "leaf_mass_fresh_g" = "Leaf.mass.fresh.g",
    "leaf_mass_faa_g"  = "Leaf.mass.FAA.g",
    "leaf_mass_dry_g" = "Leaf.mass.dry.g",
    "leaf_area_cm2" = "Leaf.area",
    "lma" = "LMA",
    "comp_seed_same_sp" = "Comp.seed.same.sp",
    "comp_seed_diff_sp" = "Comp.seed.diff.sp",
    "comp_seed_total" = "Comp.seed.total",
    "comp_adult_total" = "Comp.adult.total",
    "comp_adult_metric" = "Comp.adult.metric",
    "comp_adult_metric_log" = "Comp.adult.log.metric",
    "leaf_height_ratio" = "Leaf.height.ratio",
    "leaf_chl" = "leaf_chl",
    "lai_scale" = "LAI.4.ring_scale",
    "comp_seed_total_scale" = "Comp.seed.total_scale",
    "comp_adult_metric_log_scale" = "Comp.adult.log.metric_scale",
    "elev_scale" = "Elevation_scale",
    "num_leaves" = "num_leaves")

# Create traits only dataframe, for plotting and analysis
seedlings_traits <- seedlings_clean %>%
  dplyr::select(species, elev_code, elev,
    leaf_area_cm2, leaf_mass_dry_g, stem_vol_cm3, lma, 
    leaf_height_ratio, leaf_thick_mean_mm, leaf_chl, d_fvfm) %>%
  group_by(species, elev_code) 

# Create labels for plotting trait data
traits_levels <- c("d_fvfm", "leaf_height_ratio", "leaf_area_cm2", "leaf_chl", 
  "leaf_mass_dry_g", "leaf_thick_mean_mm", "lma", "stem_vol_cm3")
  
traits_labels <- c(
  expression("Chlorophyll-"*alpha), 
  expression("D" ~ F[v] / F[m]), 
  expression("Leaf:height" ~ "ratio" ~ (n ~ cm^-1)), 
  expression("Leaf" ~ "area" ~ (cm^2)),
  expression("Leaf" ~ "biomass" ~ "(g)"),
  expression("LMA" ~ (g ~ cm^-2)),
  expression("Mean" ~ "leaf" ~ "thickness" ~ "(mm)"),
  expression("Stem" ~ "vol." ~ (cm^3))
  )

# Summarise data by taking mean trait values per site
seedlings_traits_elev_code_summ <- seedlings_traits %>%
  summarise_all(list(mean), na.rm = TRUE) %>%
  gather(key = "var", value = "value", -species, -elev, -elev_code) %>%
  mutate(var_exp = factor(var, 
    levels = traits_levels,
    labels = traits_labels))

# Create long format dataframe for facet wrapping trait data
seedlings_traits_gather <- seedlings_traits %>%
  gather(key = "var", value = "value", -species, -elev, -elev_code) %>%
  mutate(var_exp = factor(var, 
    levels = traits_levels,
    labels = traits_labels))

# Create labels for plotting environmental data
env_levels <- c("soil_temp_mean", "soil_mois_mean",
      "lai", "comp_seed_total", "comp_adult_total", "comp_adult_metric_log")

env_labels <- c(
  expression("Mean" ~ "soil" ~ "temp" ~ (degree * C)),
  expression("Mean" ~ "soil" ~ "mois" ~ ('%')),
  expression("LAI" (m^2 ~ m^{-2})),
  expression("Herb." ~ "plant" ~ "abundance" (n ~ 3.14 ~ m^{-2})),
  expression("Adult" ~ "tree" ~ "abundance" (n ~ 78.54 ~ m^{-2})),
  expression("Iterative" ~ "Seedling" ~ "Index"))

# Create long format dataframe for facet wrapping environmental data
seedlings_env_gather <- seedlings_clean %>%
  dplyr::select(species, elev,
    soil_temp_mean, soil_mois_mean, 
    lai, comp_seed_total, comp_adult_total, comp_adult_metric_log) %>%
  gather(key = "var", value = "value", -species, -elev) %>%
  mutate(var_exp = factor(var, 
    levels = env_levels,
    labels = env_labels))
  
  
# Interaction plots of plant traits over elevation codes
spaghetti <- ggplot(seedlings_traits_elev_code_summ, aes(x = elev_code, y = value, 
  group = species, colour = species)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~var_exp, scales = "free_y", labeller = label_parsed) + 
  theme_classic() +
  theme(panel.grid.major.x = element_line(colour = "grey"))

ggsave(file="../manuscript/img/spaghetti.pdf", plot=spaghetti, width=10, height=5)

# Boxplots of variation in plant traits across species

box <- ggplot(seedlings_traits_gather, aes(x=species, y=value)) + 
  geom_boxplot(aes(fill = species), colour = "black") + 
  facet_wrap(~var_exp, scales = "free_y", labeller = label_parsed) + 
  theme_classic()

ggsave(file="../manuscript/img/box.pdf", plot=box, width=10, height=5)

# Table of where species are sampled
species_site_summ <- species_site_elev %>%
  dplyr::select(-range) %>%
  spread(key = position, value = site)

fileConn <- file("../manuscript/include/species_sample_loc.tex")
writeLines(stargazer(species_site_summ, 
  summary = FALSE, 
  rownames = FALSE,
  label = "species_sample_loc", digit.separate = 0), fileConn)
close(fileConn)

# How many individuals sampled at each elevation 
species_elevcode_tally <- data.frame(table(seedlings_clean$species, seedlings_clean$elev_code) [,])

fileConn <- file("../manuscript/include/species_elevcode_tally.tex")
writeLines(stargazer(species_elevcode_tally, 
  summary = FALSE, rownames = FALSE,
  label = "species_elevcode_tally", digit.separate = 0), fileConn)
close(fileConn)

# Genus level migration rates plot
mig <- ggplot(genus_mig_rates, aes(x = genus, y = mig_rate_abund_m, fill = genus)) + 
  geom_bar(stat = "identity") + 
    theme_classic() + 
    labs(x = "Genus", y = expression("Migration" ~ "rate" ~ (m ~ y^-1))) + 
    theme(legend.position = "none")

ggsave(file="../manuscript/img/mig.pdf", plot=mig, width=10, height=5)

# Site level environmental variables (Whitaker?)
fileConn <- file("../manuscript/include/site_char.tex")
writeLines(stargazer(site_char, summary = FALSE, type = "latex", 
  rownames = F, label = "site_char", digit.separate = 0), fileConn)
close(fileConn)

# Plot ranges and sample locations
ranges_spread <- ranges %>%
  spread(max_min, range)

ranges_ggplot <- ggplot() + 
  geom_abline(aes(intercept = elev_mean, slope = 0), 
    linetype = 2, size = 0.5, 
    data = camp_loc) +
  geom_point(aes(x = species, y = range, colour = position), 
    size = 5, 
    data = species_site_elev) + 
  geom_point(aes(x = species, y = range),
    shape = 15, size = 4, 
    data = ranges) + 
  geom_segment(aes(x = species, y = min, xend = species, yend = max), 
    data = ranges_spread) + 
  xlab("Species") + 
  ylab("Elevation (m)") + 
  scale_colour_discrete(name = "Plot position") +
  theme_classic()

ggsave(file = "../manuscript/img/ranges.pdf", plot = ranges_ggplot, width = 10, height = 5)

# Trait variables scatterplots with linear model fits
traits_elev_scatter <- ggplot(seedlings_traits_gather, aes(x=elev, y = value, colour = species)) + 
    geom_point() + 
    geom_smooth(aes(fill = species, colour = species), method = lm, se = T) + 
    facet_wrap(~var_exp, scales = "free", labeller = label_parsed) + 
    theme_classic() + 
    labs(y = "")

ggsave(file = "../manuscript/img/traits_elev_scatter.pdf", plot = traits_elev_scatter, width = 10, height = 5)

# Trait variable slopes with error bars

## Split gathered df into list by trait, 
## then by species, 
## then unlist to flat list of dfs
seedlings_traits_gather_list <- split(seedlings_traits_gather, seedlings_traits_gather$var)

seedlings_traits_gather_list_species <- lapply(seedlings_traits_gather_list, function(x){
  split(x, x$species)
})

seedlings_traits_gather_list_species_unlist <- unlist(seedlings_traits_gather_list_species, recursive = FALSE)

## Create linear models
seedlings_traits_gather_list_species_unlist_mods <- lapply(seedlings_traits_gather_list_species_unlist, function(x){
  ifelse(all(is.na(x$value)) == TRUE, 
    return(NA), 
    return(lm(value ~ elev, data = x)))
})

## Extract slopes and standard errors
lm_slope <- unlist(lapply(seedlings_traits_gather_list_species_unlist_mods, function(x){
  ifelse(class(x) == "lm",
    return(summary(x)$coeff[2]),
    return(NA))
}))

lm_se <- unlist(lapply(seedlings_traits_gather_list_species_unlist_mods, function(x){
  ifelse(class(x) == "lm",
    return(summary(x)$coeff[, 2][2]),
    return(NA))
}))

## Combine into data frame
lm_df <- data.frame(mod = names(seedlings_traits_gather_list_species_unlist_mods), 
  slope = lm_slope, se = lm_se)

## Add columns for response and species
lm_df$response <- gsub("\\..*", "", lm_df$mod)

lm_df$species <- gsub(".*\\.", "", lm_df$mod)

## Plot
traits_elev_slopes <- ggplot() + 
  geom_point(data = lm_df, 
    aes(x = species, y = slope, colour = species),
    size = 2) + 
  geom_errorbar(data = lm_df, 
    aes(x = species, ymin = slope - se, ymax = slope + se, colour = species),
    width = 1) + 
  geom_hline(yintercept = 0, linetype = 5) + 
  facet_wrap(~response, scales = "free_y") + 
  theme_classic() + 
  theme(legend.position = "none")

ggsave(file = "../manuscript/img/traits_elev_slopes.pdf", plot = traits_elev_slopes, width = 10, height = 5)

# Competition radius calculations

## Calculate competition radius
k <- 2
comp_radius <- data.frame("site" = camp_loc$site, "elev_mean" = camp_loc$elev_mean, k)
comp_radius <- left_join(comp_radius, dplyr::select(site_char, Site_Code, trees_ha), by = c("site" = "Site_Code"))
comp_radius$c_r <- comp_radius$k * (sqrt(10000 / comp_radius$trees_ha))

## Linear regression of trees / ha vs elevation to get VC
comp_radius_vc <- filter(comp_radius,  
  site %in% c("PA400", "PA800", "SP1500", "SP1750", "TRU08",
    "TRU07", "TRU06", "TRU04", "TRU02"))

trees_ha_elev <- lm(trees_ha~elev_mean, data = comp_radius_vc)

## Make a plot of with a linear regression - filling in VC
comp_radius_fit <- ggplot(comp_radius, aes(x = elev_mean, y = trees_ha)) + 
  geom_point() +
  geom_smooth(method = lm) + 
  geom_text(aes(label = site), 
    colour = c(rep("black", 2),"red",rep("black", 7)),
    hjust =-0.2, size = 4) + 
  scale_x_continuous(limits = c(0, 3600)) + 
  ylab(expression(Trees~ha^"-1")) + xlab("Elevation (m)") + 
  theme_classic()

ggsave(file = "../manuscript/img/comp_radius_fit.pdf", plot = comp_radius_fit, width = 10, height = 5)

## Pretty table of competition radius per site
comp_radius$c_r <- signif(comp_radius$c_r, digits = 1)
comp_radius_k <- dplyr::select(comp_radius, -k, -elev_mean)

fileConn <- file("../manuscript/include/comp_radius.tex")
writeLines(stargazer(comp_radius_k, summary = FALSE, type = "latex",
  rownames = F, label = "comp_radius", digit.separate = 0), fileConn)
close(fileConn)

# Plot rank abundance curve with species measured highlighted
## Create summary data frame 
aberg_census <- aberg_census %>%
  mutate(genus_species = paste(genus, specie))

aberg_census_summ <- aberg_census %>%
  filter(cod %in% c("TRU-02", "TRU-04", "TRU-06", 
    "TRU-07", "TRU-08", "SPD-01", 
    "SPD-02", "PAN-01", "PAN-02")) %>%
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
      T, F),
    myrcia = case_when(genus_species == "Myrcia splendens" | 
        genus_species == "Myrcia fallax" | 
        genus_species == "Myrcia rostrata" ~ T))

rank_abund <- ggplot() + 
  geom_point(data = filter(aberg_census_summ, sampled == T), 
    aes(x = id, y = n), 
    colour = "#E03A50FE", alpha = 1, size = 5) +
  geom_point(data = filter(aberg_census_summ, myrcia == T), 
    aes(x = id, y = n), 
    colour = "#53C90E", alpha = 1, size = 5) +
  geom_point(data = filter(aberg_census_summ, sampled == F), 
    aes(x = id, y = n), 
    colour = "black", alpha = 0.5, size = 1) + 
  theme_classic() + 
  labs(x = "Rank", y = "N")

ggsave(file = "../manuscript/img/rank_abund.pdf", plot = rank_abund, width = 10, height = 5)

# Compare random effects structure across models for future analysis

## Create predictors
predictors <- seedlings_clean %>%
  dplyr::select(species, site, comp_adult_metric_log_scale, 
elev_scale, lai_scale, comp_seed_total_scale)

## Create responses
responses <- seedlings_clean %>%
  dplyr::select(d_fvfm, leaf_chl, leaf_height_ratio, 
    leaf_area_cm2, stem_vol_cm3, leaf_thick_mean_mm)

## Create models
mod_list_intercept <- lapply(responses, function(x){
  lapply(predictors[3:6], function(y){
    lmer(x ~ y + (1|predictors$species) + (1|predictors$site),
      REML = F)
  })
})

mod_list_slope <- lapply(responses, function(x){
  lapply(predictors[3:6], function(y){
    lmer(x ~ y + (y|predictors$species) + (1|predictors$site),
      REML = F)
  })
})

## Rename models
for(x in 1:length(mod_list_intercept)){
  names(mod_list_intercept[[x]]) <- paste0(names(mod_list_intercept[[x]]), "-", names(mod_list_intercept[x]), "-intercept")
  }

for(x in 1:length(mod_list_slope)){
  names(mod_list_slope[[x]]) <- paste0(names(mod_list_slope[[x]]), "-", names(mod_list_slope[x]), "-slope")
}


## Collapse list of lists and rename
mod_list_intercept_collapse <- unlist(mod_list_intercept)
mod_list_slope_collapse <- unlist(mod_list_slope)

## Random effects models for each of the response variables
rand_mod_list <- lapply(responses, function(x){
  lmer(x ~ (1|predictors$species) + (1|predictors$site),
    REML = F)
  })

names(rand_mod_list) <- paste0("NA-", names(rand_mod_list))

## Combine models
mod_list_all <- c(mod_list_intercept_collapse, mod_list_slope_collapse, rand_mod_list)

## Clean names
names(mod_list_all) <- gsub(".*\\.", "", names(mod_list_all))

## Add AIC and convert to data frame
mod_output <- as.data.frame(t(as.data.frame(lapply(mod_list_all, function(i) AIC(i)))))
names(mod_output)[names(mod_output) == 'V1'] <- "AIC"

## Convert rownames to columns
mod_output <- rownames_to_column(mod_output, "model")

## Create column of response variables, fixed effects, and slope/intercept for grouping
mod_output$fixed_eff <- gsub("\\..*", "", mod_output$model)

mod_output$response <- str_extract(mod_output$model, "(?<=\\.)(.*?)(?=\\.)|(?<=\\.)(.*?)$")

mod_output$rsri <- gsub(".*\\.", "", mod_output$model)
mod_output$rsri <- ifelse(mod_output$rsri %in% c("intercept", "slope"),
  mod_output$rsri, 
  "intercept")

## Calculate fit statistics
mod_output <- mod_output %>%
  group_by(response) %>%
  mutate(akaike_weight = Weights(AIC)) %>%
  ungroup() %>%
  mutate(r2m = unlist(unname(lapply(mod_list_all, function(i) r.squaredGLMM(i)[1]))),
    r2c = unlist(unname(lapply(mod_list_all, function(i) r.squaredGLMM(i)[2]))))

## Extract fixed effect slopes with standard error bars
mod_list_summ <- lapply(mod_list_all, function(i){
  summary(i)
  })

get_coeffs <- function(x){
  stopifnot(x$objClass == "lmerMod")
  slope = x$coefficients[2]
  se = x$coefficients[4]
  
  coeffs = c(slope, se)
  return(coeffs)
}

mod_list_coeff <- lapply(mod_list_summ, get_coeffs)

mod_output$model_slope <- unname(unlist(lapply(mod_list_coeff, "[[", 1)))

mod_output$model_se <- unname(unlist(lapply(mod_list_coeff, "[[", 2)))

## Mark random effects models
mod_output$is_rand <- ifelse(mod_output$fixed_eff == "NA", TRUE, FALSE)


# Choose slope or intercept, for each model
daic_intercept_slope <- mod_output %>%
  filter(is_rand == FALSE) %>%
  mutate(fixed_eff_response = paste0(fixed_eff, "-", response)) %>%
  dplyr::select(-akaike_weight, -r2m, -r2c, -model_slope, 
    -model_se, -is_rand, -fixed_eff, -response, -model) %>%
  spread(key = fixed_eff_response, value = AIC) %>%
  column_to_rownames(var = "rsri") %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "fixed_eff_response") %>%
  rename(aic_intercept = intercept, 
    aic_slope = slope) %>%
  mutate(
    fixed_eff = c(
      rep("comp_adult_metric_log_scale", times = 6),
      rep("comp_seed_total_scale", times = 6),
      rep("elev_scale", times = 6),
      rep("lai_scale", times = 6)
      ),
    response = c(
      rep(c("d_fvfm", "leaf_area_cm2", "leaf_chl", 
        "leaf_height_ratio", "leaf_thick_mean_mm", "stem_vol_cm3"), times = 4)
    ),
    daic_rsri = abs(aic_intercept) - abs(aic_slope)) %>%
  mutate(best_model = case_when(
    daic_rsri > 0 ~ "slope",
    daic_rsri < 0 ~ "intercept"
  )) %>%
  mutate(best_model_name = paste0(fixed_eff_response, "-", best_model)) %>%
  mutate(best_model_name = gsub("-", ".", best_model_name))


# Create new model list which chooses the optimal random effect structure for each combination of fixed and random effects

# Plot slopes of best single predictor models
## Rename factor levels for plotting
mod_output_best <- mod_output %>%
  filter(is_rand == FALSE) %>%
  filter(model %in% daic_intercept_slope$best_model_name) %>%
  mutate(response_exp = factor(response, 
    levels = c("d_fvfm", "leaf_chl", "leaf_height_ratio", 
      "leaf_area_cm2", "stem_vol_cm3", "leaf_thick_mean_mm"), 
    labels = c(
      expression("D" ~ F[v] / F[m]), 
      expression("Chlorophyll-"*alpha), 
      expression("Leaf:height" ~ "ratio" ~ (n ~ cm^-1)), 
      expression("Leaf" ~ "area" ~ (cm^2)),
      expression("Stem" ~ "vol." ~ (cm^3)),
      expression("Mean" ~ "leaf" ~ "thickness" ~ "(mm)")
    )),
    fixed_eff_exp = factor(fixed_eff, 
      levels = c("comp_adult_metric_log_scale", "elev_scale",
        "lai_scale", "comp_seed_total_scale"),
      labels = c("ISI", "Elev.", "LAI", "Herb."))) %>% 
  mutate(rsri_exp = factor(rsri, 
    levels = c("intercept", "slope"),
    labels = c("Intercept", "Slope")))

## Plot the slopes of each model
single_pred_slope <- ggplot(mod_output_best, 
  aes(x = fixed_eff_exp, y = model_slope)) + 
  geom_errorbar(
    aes(ymin = model_slope-model_se, ymax = model_slope+model_se, 
    colour = fixed_eff_exp)) + 
  geom_point(aes(fill = fixed_eff_exp, shape = rsri_exp),
    colour = "black", size = 2) + 
  geom_hline(aes(yintercept = 0), linetype = 2) + 
  facet_wrap(~response_exp, scales = "free_y", labeller = label_parsed) + 
  theme_classic() + 
  theme(legend.position = "right") +
  labs(x = "Fixed effect", y = "Slope") + 
  scale_shape_manual(name = "Random\neffects\nstructure", values = c(21, 23)) +
  scale_fill_discrete(guide = FALSE) + 
  scale_colour_discrete(guide = FALSE)

ggsave(file="../manuscript/img/single_pred_slope.pdf", plot=single_pred_slope, width=10, height=5)

# Plot the R2m for each single predictor model
single_pred_r2m <- ggplot(filter(mod_output_best), 
  aes(x = fixed_eff_exp, y = r2m)) + 
  geom_bar(stat = "identity", aes(fill = fixed_eff_exp), colour = "black") +
  facet_wrap(~response_exp, scales = "fixed", labeller = label_parsed) + 
  theme_classic() + 
  theme(legend.position = "none") +
  labs(x = "Fixed effect", y = expression(R[m]^2))

ggsave(file="../manuscript/img/single_pred_r2m.pdf", plot=single_pred_r2m, width=10, height=5)

# Calculate the dAIC against a null model for each best single predictor

rand_mod_aic <- mod_output %>%
  filter(is_rand == TRUE) %>%
  dplyr::select(response, AIC) %>%
  rename(AIC_rand = AIC)

mod_output_best <- mod_output_best %>%
  left_join(., rand_mod_aic, by = "response") %>%
  mutate(daic_rand = AIC_rand - AIC)

single_pred_daic <- ggplot(mod_output_best, 
  aes(x = fixed_eff_exp, y = daic_rand)) + 
  geom_bar(stat = "identity", aes(fill = fixed_eff_exp), 
    colour = "black") + 
  facet_wrap(~response_exp, scales = "free_y", labeller = label_parsed) + 
  theme_classic() + 
  theme(legend.position = "none") +
  labs(x = "Fixed effect", y = expression(delta*"AIC"[r]))  
  
ggsave(file = "../manuscript/img/single_pred_daic.pdf", plot = single_pred_daic, width = 10, height = 5)

# Find best multiple predictor model for each response variable
