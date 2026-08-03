# bs_trait_dtrmnts_repro_success_mediation

# Manuscript Citation
Laubach, Z. M., Keller, K. P., Safran, R. J., Tsunekage, T., Levin, I. I.,  2026. Testing the mediating role of social interactions on the relationship between age and plumage traits on reproductive success. _Ethology_

----

### Author names, affiliations and contact information
* Zachary M. Laubach, 1, 2
* Email: zachary.laubach@colorado.edu

* Kayleigh P. Keller, 2
* Rebecca J. Safran, 1
* Rebecca J. Safran, 1
* Toshi Tsunekage, 3
* Iris I. Levin, 3

* 1 Department of Ecology and Evolutionary Biology, University of Colorado, Boulder, CO, USA
* 2 Department of Statistics, Colorado State University, Fort Collins, CO, USA
* 3 Department of Biology, Kenyon College, Gambier, OH, USA


----

### Paper summary:
Differential reproduction is a key driver of evolution that is determined by individual characteristics and mating opportunities, including mate choice. Social interactions between conspecifics are hypothesized to be important in facilitating mate choice and reproductive success but are difficult to measure. Using data from 52 adult barn swallows (Hirundo rustica erythrogaster), whose social interactions were measured via proximity tags, we tested the hypothesis that social interactions mediate both the relationship between age (a proxy for experience) as well as condition-dependent plumage traits, and their associations with reproductive success. We found that older female barn swallows had higher fecundity and that older males have higher paternity. Older males achieved higher paternity through extra pair copulations, not by greater paternity with their social mate. Longer tail streamers were associated with greater fecundity/paternity in both sexes, but this effect was independent of age only among females. Darker ventral plumage coloration was not associated with higher reproductive success in either sex. We observed that older males appear to be less social with conspecifics, as indicated by fewer numbers of social interactions, though these associations were only marginally significant in males. The size of the effect of female age on conspecific social strength was similar to males, but was not significant in females. Interestingly, females with fewer social interactions had higher fecundity. Finally, we found no evidence of mediation by the number of social interactions. Taken together, our results suggest that older, more experienced birds can produce more offspring while being less social. 

![Data collection](/cover_image.png "cover image")


----

### File organization and description

#### Purpose: This repository contains the scripts, data, and output necessary to reproduce the analysis for this paper. Scripts are modular and numbered in the order they are to be run. An overview of the repository organization and file descriptions are provided below.

#### There are 4 subdirectories in this repository.

* _data -_ this subdirectory contains the raw data (.csv) and the processed data (.RData) for all analyses. The .RData follow the modular scripts such that the processed data saved at the end of one script is loaded in the sequential script, thus allowing the user to enter the analysis pipeline at an stage without having to re-run all previous steps.

* _source_code -_ subdirectory contains any custom functions that are used in the scripts.

* _scripts -__ this subdirectory contains the sequentially numbered R scripts for load and cleaning the data and for all downstream analyses. Given that the scripts are modular, they can be run alone. For complete reproducibility run the scripts in sequential order starting with script 1.

* _output -__ this subdirectory contains summary tables and figures generated as part of the analysis.

#### Data description
* _chr15_attrib_df.R -_ This file has 52 rows and 50 columns that contain the tidy data used in the main analyses. Each row corresponds to an indidivual adult barn swallow. The columns include four data types, including 9 characater (chr) variables, 38 numeric variables (num), and 3 factor variables (fact). Below are the metadata for the variables that are used in the analyses presented in th is paper. Other variables in the dataframe are not used in this paper and are not explicitly defined
	-	Band.ID = chr, the individual bird's ID
	-	Tag = chr, the proximity tag ID
	-	Year = num, the sampling year
	-	Site = chr, the breeding site location
	-	Sex = fact, each bird's sex, female(f) or male(m)
	-	Age.category = fact, each bird's categorical age, sy - second year (1yr old, first time breeders) and asy - after second year (2 or more years old,  multiple year breeders
	-	Age = num, each birds's age in years
	-	Nest = chr, the nest ID
	-	Mean.TS = num, the average of 3 right tail streamer length measures (mm)
	-	B_avg.bright = num, the average of 3 measures of light reflectance from natural belly feathers; darker plumage has a lower spectral signal
	-	T_avg.bright = num, the average of 3 measures of light reflectance from natural throat feathers
	-	R_avg.bright = num, the average of 3 measures of light reflectance from natural breast feathers
	-	Post.R_avg.brigth = num, the average of 3 measures of light reflectance from breast feathers after the experimental ventral darkening treatment
	-	R.bright.treat.and.orig = num, concatenation of the variables R_avg.bright and Post.R_avg.brigth
	-	soc.mate.ID.nest.1 = chr, the ID of each bird's social mate at their first nest
	-	soc.mate.tag.nest.1 = chr, the proximity tag ID of each bird's social mate at their first nest
	-	total.fecundity = num, the total number of eggs a female laid during the breeding season (from first through third nesting attempts)
	-	nest.2.3.tot.pat = num, a male's total number of genetically determined offspring (from second through third nesting attempts)
	-	nest.2.3.epp = num, a male's extra-pair genetically determined offspring (from second through third nesting attempts)
	-	nest.2.3.spp = num, a male's within social-pair genetically determined offspring (from second through third nesting attempts)
	-	strength.fxm = num, node-level strength quantified as the the total number of social interactions from proximity based female-male social networks 
	-	degree.fxm = num, node-level degree quantified as the the total number of social partners from female-male proximity based social networks 
	-	strength.fxf = num, node-level strength quantified as the the total number of social interactions from proximity based female-female social networks 
	-	degree.fxf = num, node-level degree quantified as the the total number of social partners from female-female proximity based social networks 
	-	strength.mxm = num, node-level strength quantified as the the total number of social interactions from proximity based female-female social networks 
	-	degree.mxm = num, node-level degree quantified as the the total number of social partners from female-female proximity based social networks 


#### Script description

* _0_import_raw_data.R -_ Load barn swallow parental care and offspring physiology data. The 0_import_data script simply loads the raw data to be used in the analysis. If data are pulled from a remote location, skip this step and proceed to script 1.

* _1_tidy_attrib_data.R -_ Tidy and join phenotyic and demographic data. 

* _2_tidy_repro_data.R -_ Tidy and join reproductive data.

* _3_tidy_intx_data.R -_ Tidy and join social interaction data and build adjacency matrices.

* _4_quant_soc_intx.R -_ Quantify the social interaction data to create node and level network measures from pre- manipulation social networks

* _5_descrprtv_stats_node_data.R -_ Calculate descriptive statistics and explore data 

* _6_med_and_tot_effects_models.R -_ Run mediation and total effects models for traits and reproductive success


#### Output description

#### Plot used to create Figure 3
* _female_panel_plot.pdf -_ Dot and whisker plot showing estimates from the female barn swallow causal mediation models for a) the effect of age on fecundity and mediation by strength from female-male interaction networks and b) the effect of tail streamer length on fecundity and mediation by strength from female-male interaction networks. Estimates include the total effect estimated from each model and its decomposition into the natural direct effect (aka, average direct effect; ADE) and the natural indirect effect (aka, average causal mediation effect; ACME).


----

### Software and version information

#### This analysis was performed in R
*	R version 4.6.0 (2026-04-24)
*	Platform: x86_64-apple-darwin20
*	Running under: macOS Tahoe 26.5.2

#### R package versions
attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
assortnet_0.20       igraph_2.3.2        gridExtra_2.3.1     geomtextpath_0.2.0  lmPerm_2.1.6         performance_0.17.0  rempsyc_0.2.0       broom_1.0.13        broom.mixed_0.2.9.7  mediation_4.5.1     sandwich_3.1-1      mvtnorm_1.4-1       Matrix_1.7-5         MASS_7.3-65         patchwork_1.3.2     here_1.0.2          lubridate_1.9.5      forcats_1.0.1       stringr_1.6.0       dplyr_1.2.1        
purrr_1.2.2          readr_2.2.0         tidyr_1.3.2         tibble_3.3.1        ggplot2_4.0.3        tidyverse_2.0.0 	 naniar_1.1.0

