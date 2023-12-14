# Summary
These are MCMC samplers and processing and run scripts for an inhomogenous Poisson point process model to integrate detection/non-detection data and presence-only data to estimate the expected abundance of species. This model enables the user to intergrate these two datatypes by assuming they share the same underlying data generating process (an inhomogeneous Poisson point process). These models are adapted from Koshkina et al. (2017) [Methods in Ecology and Evolution] but diverge from the original specification to accomodate _t_ primary sampling periods (survey years). This model fits temporally varying mean-centered random year effects of all three parameters in the model (detection probability, p, mean expected abundance, λ, and thinning, b). 


# The working directory

Below you fill find descriptions of each folder in this repository and files contained within them.

# The data folder (./data)

This folder has x files.

**1. allNY_10kmgrid_alllandscapecovs_final.csv**

This file contains all of the summarized spatial covariate data used in the analysis. Each row is a different 10km 2 pixel in New York State, each column is a covariate, and each cell is a value.

| **Covariate**           | **Description**                                                                       | **Source**                                                                                              |
|-------------------------|---------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------|
| Deciduous               | Proportion of a 10 km2 grid cell made up of deciduous forest                          | NLCD, 2019 (https://www.mrlc.gov/data/nlcd-2019-land-cover-conus)                                       |
| Coniferous              | Proportion of a 10 km2 grid cell made up of coniferous forest                         | NLCD, 2019 (https://www.mrlc.gov/data/nlcd-2019-land-cover-conus)                                       |
| Mixed                   | Proportion of a 10 km2 grid cell made up of mixed forest                              | NLCD, 2019 (https://www.mrlc.gov/data/nlcd-2019-land-cover-conus)                                       |
| Pasture                 | Proportion of a 10 km2 grid cell made up of pasture                                   | NLCD, 2019 (https://www.mrlc.gov/data/nlcd-2019-land-cover-conus)                                       |
| Cultivated.Crops        | Proportion of a 10 km2 grid cell made up of cultivated crops                          | NLCD, 2019 (https://www.mrlc.gov/data/nlcd-2019-land-cover-conus)                                       |
| road_density            | Mean number of km of road per km 2 in each grid cell                                  | calculated from primary and secondary roads raster provided by the NYSDEC, hosted on the githu          |
| elevation               | Mean elevation (m) of the 10 km2 grid cell                                            | calculated from Digital Elevation Models of New York State provided by the NYSDEC, hosted on the github |
| forest_edge             | Edge density of combined class of all forest                                          | NLCD, 2019 (https://www.mrlc.gov/data/nlcd-2019-land-cover-conus)                                       |
| PO_join_count_truncated | density of presence-only points (all mammals from online repositories from 2013-2021) | calculated from all online repositories, hosted on github                                               |

**2. allNY_2016-2018_bearPAdata_pixelID_10kmgrid_final.csv**

This file contains all the detection/non-detection data for black bears collected during summer surveys from 2016-2018 in the southern tier of New York State. Each row is a site, each column is an occasion, each cell is a detection/non-detection record.

**3. allNY_2016_2018_bearPOrecord_pixelID_10kmgrid_noweeklydups_final.csv**

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


# The model folder (./models)

The models are coded up via the nimble package version 0.13.1 (de Valpine et al. 2022). 
Dec
**1. nimble_ppp_integrated_model_bear_10km_yearindexed_randomeffects_centrered.R**

This is the nimble IPPP model that is fit to the black bear data files above. The code is commented out to describe each part of the model.

**2. nimble_ppp_integrated_model_coyote_10km_yearindexed_randomeffects_centered**

This is the nimble IPPP model that is fit to the coyote data files above. The code is commented out to describe each part of the model.

**3. nimble_ppp_integrated_model__bobcat_10km_indexingoveryear_randomeffects_centered**

This is the nimble IPPP model that is fit to the bobcat data files above. The code is commented out to describe each part of the model.


# The scripts folder (./scripts)

**1. bear_PO_PA_model_10km_nimble_formattingandrunscript_indexoveryear**

This is the formatting and run script for the bear data and model above. This code is commented throughout.

**2. coyote_PO_PA_model_10km_nimble_formattingandrunscript_indexoveryear**

This is the formatting and run script for the coyote data and model above. This code is commented throughout.

**3. bobcat_PO_PA_model_10km_nimble_formattingandrunscript_indexoveryear** 

This is the formatting nad run script for the bobcat data and model above. This code is commented throughout. 
