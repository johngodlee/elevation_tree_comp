# SEECC presentation script

* 1
	* Introduction
	* I’m going to discuss a project where I investigated the potential role of competition and forest structure in limiting the elevational range shifts of tropical tree species.
 * 2
 	* We know that climate change is causing an increase in temperature, changes in precipitation regime, and an increase in the frequency of extreme weather events.
	* Populations of plants can adapt to these changes in climate by altering their morphology (i.e. their average leaf thickness), their phenology (the timing of life events like flowering and seeding), or their physiology (i.e. changing their photochemistry).
	* If species fail to adapt to climate change fast enough, this will result in higher mortality in some areas, and higher recruitment in others. Causing range shifts, either to higher elevations or higher latitudes.
	* A key area of ecological research is in predicting how species will move, whether their ranges will contract or expand, and whether ecosystem functioning may change as a result.
 * 3	
 	* The majority of studies that have aimed to do this in the past have used bioclimatic envelope models, which use current species distributions, which are then correlated with climate variables and projected into the future to predict how a species range with alter.
* 4 
	* But bioclimatic envelope models don’t always produce accurate results, especially when tracking range shifts across elevation.
	* These graphs  published as part of a data synthesis investigating global trends of species range shifts show that models predicting latitudinal range shifts are doing a good job at explaining the observed trends, but models over elevation predict a lot of upslope movement that we aren’t seeing in real life.
* 5
	* A potential reason for the elevational models not matching up with observed range shifts is because they don't include non-climatic variables that co-vary with elevation such as herbivory, competition interactions, and other aspects of the biotic environment.
 * 7
 	* So, in this study I focused specifically on whether changes in forest structure can explain variation in plant stress across elevation. 
	* We would expect that those individuals experiencing stress are less likely to survive to produce healthy seeds and so could present a barrier to upslope migration.
	* From this we can make a judgement as to whether proxies for forest structure could be included in these predictive range shift models to improve their accuracy.
* 8
	* To answer this question I worked on a transect of permanent survey plots on the Eastern flank of the Andes, where previous studies have identified the area as a potential upslope migration corridor.
	* Along this transect, which runs through lowland wet forest at about 400 m through cloud forest and up to the tree line at about 3200 m, I chose 10 plots in which to sample the seedlings of nine tree species.
* 9
	* The nine species had contrasting ranges and growth forms, some dicots some monocots, but all made up a decent proportion of the biomass in their respective plots.
	* For each of these species I sampled individuals at the top, middle, and bottom of their range.
* 10
	* I quantified seedling stress levels by measuring chlorophyll fluorescence and chlorophyll content.
	* Chlorophyll fluorescence can be used as a proxy for photosynthetic efficiency, i.e. how effectively light energy is used by the plant. 
	* When light hits a chloroplast, the majority of it is absorbed, and from there the energy can be used to drive the light reactions of photosynthesis, or be dissipated as heat, or as fluorescence.
	* In a healthy plant, we would expect a higher proportion to be used in photosynthesis...
* 11
	* But in a stressed plant, where some of the photosynthetic machinery has been shut down because of UV damage or something, we would expect less to be used in photosynthesis, and therefore more to be dissipated as fluorescence, which we can then measure.
	* A decrease in chlorophyll content is associated with various forms of plant stress including temperature and nutrient stress, and a decrease in chlorophyll content limits overall photosynthetic capacity.
	* Both of these measures can be used to infer whether a plant will be successful over its lifetime in terms of fitness.
* 12
	* Around each seedling, I also measured three proxies to describe the forest structure. 
	* I measured canopy density using the Leaf Area Index calculated from a hemispherical photo of the canopy above each seedling
	* I measured adult-seedling root competition, using  a metric that takes into account the distances and trunk diameters of competitor trees.
	* and I measured seedling abundance by counting the number of seedlings around the seedling I measured.
* 13
	* To determine whether forest structure affects either of the plant stress variables and in what way, I conducted a series of linear mixed effects models, each using one competition variable.
	* I included random effects of plot and species to account for the repeated measurements at each plot and to account for baseline differences in the physiology of each species I sampled.
* 14
	* To see what combination of forest structure variables best capture the variation in plant stress, I compared models with various combinations of fixed effects using AIC and pseudo-r-squared estimates.
* 15
	* The graph on the left shows the model quality of each of the single fixed effect models,  with models above the red line being better than null models, and the graph on the right shows how much of the variance in the stress variables is explained by each fixed effect.
	* For photosynthetic efficiency, elevation, root competition and canopy density accounted for a good deal of the variation in photosynthetic efficiency while seedling abundance accounted for very little variance.
	* Only the seedling abundance model for chlorophyll content was any better than a random effects model and this still only explained a very small amount of the variation in chlorophyll content, so chlorophyll content remains fairly unaffected by the forest structure variables.
* 16
	* This plot shows the fixed effect slopes for each of these models.
	* Firstly, it's easy to see that chlorophyll content wasn't really affected by any of these competition variables, and it might be that these species just don't seem to undergo chlorosis as a stress response, or that they are not stressed enough to have lost chlorophyll.
	* But for photosynthetic efficiency we have some interesting results.
	* As root competition increases, the photosynthetic efficiency decreases and this might be because adult trees can poach nutrients from nearby seedlings, preventing them from synthesising the necessary photosynthetic machinery.
	* Conversely, as the canopy becomes more dense, photosynthetic efficiency increases. And this might be because as the canopy becomes more dense, the micro-climate becomes more stable, at least on a diurnal scale, and there is less chance of rapid leaf temperature increases as a sun-fleck passes over. So in this way, a dense canopy facilitates the early growth of these seedlings.
* 17
	* The best model for photosynthetic efficiency, explaining the most variance while remaining as parsimonious as possible included elevation and both of the adult interactions, root competition and canopy density, but not seedling interactions, and it explained over 50% of the variance in photosynthetic efficiency.
	* While chlorophyll content was best explained by a random effects model, with none of the variables, not even elevation, having a real effect.
* 18	* In summary,	* We found from the model selection process that the adult tree forest structure variables do have an effect on seedling stress levels, and that their inclusion will likely improve the accuracy of range-shift models, given that our best model incorporated both elevation and these competition variables.
	* Considering the specific effects of forest structure on seedling stress in this particular system, we know that as you move into the cloud forest zone, at about 1500 m, the canopy becomes noticeably thinner and the trees are found much closer together owing to the slopes becoming steeper and a rapid turnover of the adult tree species composition, resulting in higher root competition, So, given the results we got from the single fixed effects models, it seems likely that this sudden transition could present a barrier to further upslope migration of lowland species, which become stressed as the canopy density decreases and root competition increases.
* 19
	* In the future I would like to see an experiment that transplants seedlings outside of their current ranges to provide more predictive data. We may have only sampled seedlings that were going to survive in the place we found them.	* In this experiment I would like to use a larger cohort of seedlings and measure them over time as they grow to see whether our variations in plant stress can be correlated with mortality of individuals	* Also I would like to conduct similar experiments on rare species as there is good evidence that rare species may be much more sensitive to variation in competitive environment but that they still perform important ecosystem functions.
	* Lastly, I am sure some people are thinking how are we supposed to get this data on forest structure, isn't it laborious and tedious to collect on the large scales we need, and I think we should be looking to what remote sensing and drone data can give us in that respect. We can now measure forest structure parameters over larger scales and much much quicker than the methods I used in this study, and I think this shows promise as a way to adequately model the effects of forest structure on climate change induced range shifts.
	* IF ANYONE WANTS TO TALK ABOUT THOSE THINGS AFTER, COME FIND ME
	*  THANKS 




