# PO-PA-Integrated-SDM
An adaption of the Koshkina et al. 2017 IDSM 

These are MCMC samplers and processing and run scripts for an inhomogenous Poisson point process model to integrated detection/non-detection data and presence-only data to estimate the expected abundance of species.
These models are adapted from Koshkina et al. (2017) [Methods in Ecology and Evolution] but diverge from the original specification to accomodate _t_ sampling periods (survey years). This model fits temporally varying
mean-centered random year effects of all three parameters in the model (detection probability, p, mean expected abundance, λ, and thinning, b). 

The models are coded up via the nimble package version 0.13.1 (de Valpine et al. 2022). 
