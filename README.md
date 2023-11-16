# Summary
These are MCMC samplers and processing and run scripts for an inhomogenous Poisson point process model to integrated detection/non-detection data and presence-only data to estimate the expected abundance of species. This model enables the user to intergrate these two datatypes by assuming they share the same underlying data generating process (an inhomogneous Poisson point process). These models are adapted from Koshkina et al. (2017) [Methods in Ecology and Evolution] but diverge from the original specification to accomodate _t_ sampling periods (survey years). This model fits temporally varying mean-centered random year effects of all three parameters in the model (detection probability, p, mean expected abundance, λ, and thinning, b). 


# The working directory

Below you fill find descriptions of each folder in this repository and files contained within them.

# data

This folder has x files.

**1. allNY_10kmgrid_alllandscapecovs.csv**

This file contains all of the summarized spatial covariate data used in the analysis. Each row is a different 10km 2 pixel in New York State, each column is a covariate, and each cell is a value.

Deciduous = Proportion of a 10 km2 grid cell made up of deciduous forest (NLCD, 2019)
Coniferous =  Proportion of a 10 km2 grid cell made up of coniferous forest (NLCD, 2019)
Mixed =  Proportion of a 10 km2 grid cell made up of mixed forest (NLCD, 2019)
Pasture =  Proportion of a 10 km2 grid cell made up of pasture  (NLCD, 2019)
Cultivated.Crops =  Proportion of a 10 km2 grid cell made up of cultivated crops  (NLCD, 2019)
road_density = Mean number of km of road per km 2 in each grid cell (calculated from primary and secondary roads raster provided by the NYSDEC, hosted on the github)
elevation = Mean elevation (m) of the 10 km2 grid cell (calculated from Digital Elevation Models of New York State provided by the NYSDEC, hosted on the github)
forest_edge = Edge density of combined class of all forest (coniferous, mixed, and deciduous). 
PO_den = density of presence-only points (all mammals from online repositories from 2013-2021).
PO_denbear = density of presence-only points (bear sightings from iSeeMammals for sampling years)
The models are coded up via the nimble package version 0.13.1 (de Valpine et al. 2022). 

**2. allNY_2016-2018_bearPAdata_pixelID_10kmgrid_remake.csv**

This file contains all the detection/non-detection data for black bears collected during summer surveys from 2016-2018 in the southern tier of New York State. Each row is a site, each column is an occasion, each cell is a detection/non-detection record.

**3. allNY_2016_2018_bearPOrecord_pixelID_10kmgrid_noweeklydups_noinaturalist.csv**

This file contains all the presence-only data used in this analysis for black bears  (maximum 1 per 10 km2 pixel per week) between the years 2016-2018 both from public online repositories and iSeeMammals, a citizen science bear monitoring project run by Cornell University

**4. NY_blackbears_ordinaldates.csv**

This file contains the sampling dates (in ordinal format) for each sampling occasion from the summer surveys from 2016-2018. Each row is a site, each column is an occasion.

**5. allNY_2013-2021_coyotePArecords_pixelID_10kmgrid.csv**

This file contains all the detection/non-detection data for coyotes collected during winter surveys from 2013-2021 in the southern tier of New York State. Each row is a site, each column is an occasion, each cell is a detection/non-detection record.

**6. allNY_2013_2021_coyotePOrecord_pixelID_10kmgrid_noweeklydups.csv**

This file contains all the presence-only data used in this analysis for coyotes (maximum 1 per 10 km2 pixel per week) between the years 2013-2021 both from public online repositories

**7. 'juliandays_allNY_2013-2021.csv'**

This file contains the sampling dates (in ordinal format) for each sampling occasion for the winter surveys from 2013-2021. Each row is a site, each column is an occasion.

**8. allNY_2013-2021_bobcatPArecords_pixelID_10kmgrid.csv**

This file contains all the detection/non-detection data for bobcats collected during winter surveys from 2013-2021 in the southern tier of New York State. Each row is a site, each column is an occasion, each cell is a detection/non-detection record.

**9. allNY_2013_2021_bobcatPOrecord_pixelID_10kmgrid_noweeklydups.csv**

This file contains all the presence-only data used in this analysis for bobcats (maximum 1 per 10 km2 pixel per week) between the years 2013-2021 both from public online repositories


# models

**1. nimble_ppp_integrated_model_bear_10km_yearindexed_randomeffects_centrered.R**


# scripts
