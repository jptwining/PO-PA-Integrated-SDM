setwd('C:/Users/jpt93/Documents/R/PO_PA model')

#read in detection non-detection data
PAdata <- read.csv("allNY_2016-2018_bearPAdata_pixelID_10kmgrid_remake.csv")

#read in presence-only data
POdata <- read.csv("allNY_2016_2018_bearPOrecord_pixelID_10kmgrid_noweeklydups_noinaturalist.csv")

#read in pixel covariate data
pixelcovdata <- read.csv("allNY_10kmgrid_alllandscapecovs.csv")

#read in observation cov data
date <- read.csv('NY_blackbears_ordinaldates.csv')

#remove 2016
remove2016 <- which(POdata$Year == 2016)
POdata <- POdata[-remove2016,]
table(POdata$Year)

#remove 2016
remove2016 <- which(PAdata$year == 2016)
PAdata <- PAdata[-remove2016,]
table(PAdata$year)

sum(PAdata[,3:7], na.rm=TRUE)

table(PAdata$year)
    date <- date[-remove2016,]
table(date$year)

# this will load these packages
library(jagsUI)
library(AHMbook)
library(raster)
library(mvtnorm)
library(runjags)
library(vioplot)
library(scales)
library(coda)



########################## data formatting ########################################################
#split PO into years

POdatayear <- split(POdata$FID_allNY_10km, POdata$Year)

unique(PAdata$year)
sitenames=unique(PAdata$Station)
years=unique(PAdata$year)
n.year=length(years)
J=rep(NA,n.year)
year.sites=list(2)
for(i in 1:n.year){
  J[i]=length(PAdata$Station[PAdata$year==years[i]])
  year.sites[[i]]=unique(PAdata$Station[PAdata$year==years[i]])
}
maxJ=max(J)
totalJ=sum(J)
rep(c(0, 1, 2), each = 10)
K=5


y <- PAdata[,3:7]
#3D matrix (site, occ, year)
y.new =array(0,dim=c(maxJ,K,n.year)) #new data 
for(i in 1:nrow(PAdata)){ #loop through each row
  this.year=which(years==PAdata$year[i]) #get year for this row
  this.site=which(year.sites[[this.year]]==PAdata$Station[i]) #get site for this row
  y.new[this.site,1:5,this.year]=as.numeric(PAdata[i,3:7]) #force numeric
}



for(t in 1:n.year){
  print(c(sum(y.new[1:J[t],,t], na.rm=TRUE),sum(PAdata[PAdata$year==years[t],3:7], na.rm=TRUE)))
}


K2D <- 1*(!is.na(y.new))

# sum(K2D)*7
y.new[is.na(y.new)] <- 0

#if count data, set counts great than 1 to 1
y.new[y.new>=1]<-1


y[is.na(y)] <- 0

sum(y.new)
#if count data, set counts great than 1 to 1

#create occassion matrix
occ<- y
for(i in 1:nrow(y)){
  tmp<- y[i,]
  n.samp<- sum(!is.na(tmp))
  tmp[!is.na(tmp)]<- 1:n.samp
  occ[i,]<- tmp
}

# create previous capture (behavioral effect) covariate
Jtemp <- dim(y)[1]; K <- dim(y)[2]
prev <- matrix(0,nrow=Jtemp,ncol=K)
for(i in 1:nrow(y)){
  tmp<- y[i,]
  if(sum(y[i,],na.rm=T)>0){
    first <- min(which(y[i,]>0))
    if(first<K){
      prev[i,(first+1):K]<-1
    }}}


occ <- as.matrix(occ)
prev <- as.matrix(prev)
occ.c = standardize(occ)

occyear <- cbind(occ.c, PAdata)
prevyear <- cbind(prev, PAdata)

#format occ covs
occ.new=array(0,dim=c(maxJ,K,n.year)) #new data 
for(i in 1:nrow(occ)){ #loop through each row
  this.year=which(years==occyear$year[i]) #get year for this row
  this.site=which(year.sites[[this.year]]==occyear$Station[i]) #get site for this row
  occ.new[this.site,1:5,this.year]=as.numeric(occyear[i,1:5]) #force numeric
}

for(t in 1:n.year){
  print(c(sum(occ.new[1:J[t],,t]),sum(occyear[occyear$year==years[t],1:5])))
}


#format prev cov
prev.new=array(0,dim=c(maxJ,K,n.year)) #new data 
for(i in 1:nrow(occ)){ #loop through each row
  this.year=which(years==prevyear$year[i]) #get year for this row
  this.site=which(year.sites[[this.year]]==prevyear$Station[i]) #get site for this row
  prev.new[this.site,1:5,this.year]=as.numeric(prevyear[i,1:5]) #force numeric
}

for(t in 1:n.year){
  print(c(sum(prev.new[1:J[t],,t]),sum(prevyear[prevyear$year==years[t],1:5])))
}

#dates
#read in ordinal day data, remove NAs and format
date <- date[7:11]
date <- as.matrix(date)
range(date)
date.c = standardize(date)

dateyear <- cbind(date.c, PAdata)
# date[is.na(y)] <- NA # fix the NAs

#for week with year last
date.new=array(0,dim=c(maxJ,K,n.year)) #new data 
for(i in 1:nrow(PAdata)){ #loop through each row
  this.year=which(years==PAdata$year[i]) #get year for this row
  this.site=which(year.sites[[this.year]]==PAdata$Station[i]) #get site for this row
  date.new[this.site,1:5,this.year]=as.numeric(dateyear[i,1:5]) #force numeric
}

for(t in 1:n.year){
  print(c(sum(date.new[1:J[t],,t]),sum(dateyear[dateyear$year==years[t],1:5])))
}



#create aggregated site covariates
pixelcovdata$All_Forest <- pixelcovdata[,"Coniferous"] + pixelcovdata[,'Mixed'] + pixelcovdata[,'Deciduous']
pixelcovdata$ConifMixed <- pixelcovdata[,"Coniferous"] + pixelcovdata[,'Mixed']
pixelcovdata$Wetlands <- pixelcovdata[,"Woody.Wetlands"] + pixelcovdata[,'Herbaceus.Wetlands']
pixelcovdata$Developed <- pixelcovdata[,"Developed.Low"] + pixelcovdata[,'Developed.Mid']+ pixelcovdata[,'Developed.High']
pixelcovdata$PO_den[(is.na(pixelcovdata$PO_den))] <- 0
log_PO_count <- log(pixelcovdata$PO_join_count_truncated+1)


#create combined covariate classes
pixelcovdata$All_Forest <- pixelcovdata[,"Coniferous"] + pixelcovdata[,'Mixed'] + pixelcovdata[,'Deciduous']
pixelcovdata$ConifMixed <- pixelcovdata[,"Coniferous"] + pixelcovdata[,'Mixed']
pixelcovdata$Wetlands <- pixelcovdata[,"Woody.Wetlands"] + pixelcovdata[,'Herbaceus.Wetlands']
pixelcovdata$Developed <- pixelcovdata[,"Developed.Low"] + pixelcovdata[,'Developed.Mid']+ pixelcovdata[,'Developed.High']

#transform the PO density covariate 
# #log transform
log_PO_count <- log(pixelcovdata$PO_join_count_truncated+1)
table(is.na(log_PO_count))


#standarize and scale all covariates to have mean and unit variance of 0 (standardize is function in the AHMbook package, very useful)
elevation.new = standardize(pixelcovdata$elevation)
conifmixed.new = standardize(pixelcovdata$ConifMixed)
deciduous.new = standardize(pixelcovdata$Deciduous)
forest.new = standardize(pixelcovdata$All_Forest)
crops.new = standardize(pixelcovdata$Cultivated.Crops)
edge.new = standardize(pixelcovdata$forest_edge)
wetlands.new = standardize(pixelcovdata$Wetlands)
pixelcovdata$road_density[(is.na(pixelcovdata$road_density))] <- 0
road.new = standardize(pixelcovdata$road_density)
# poden.new = standardize(log_po_den)
POcount.new = standardize(log_PO_count)
developed.new = standardize(pixelcovdata$Developed)


## large constant needed for Bernoulli one's trick
constant <- 10000

####number of pixels in landscape
G <- nrow(pixelcovdata)


## what year each PO data point came from 
npo_count <- sapply(POdatayear, length)

npo_count<- as.vector(npo_count)

opp_year = rep(1:n.year, times = npo_count)

# The number of presence only data points per season
npo = as.numeric(npo_count)
# The total number of all presence only data points
all_npo = sum(npo_count)

#checking for NAs which will cause trouble later
table(is.na(conifmixed.new))
table(is.na(wetlands.new))
table(is.na(edge.new))
table(is.na(POcount.new))
table(is.na(occ.new))
table(is.na(date.new))
table(is.na(y))

#check for colinearity between covariates
library("car")
mydata <- pixelcovdata[,c(9 ,18, 19, 20,25, 30, 33, 15)]
head(mydata)
mymodel <- lm(random~., data=mydata)
vif(mymodel) #no conlinearity 
# eps.p = eps.p, sd.p = 1, sd.p = 1, eps.p = eps.new)
# #random year level effect on psi inits
beta0.it  <- rnorm(J)
# 
# #random year level effect on psi inits
alpha0.it <- rnorm(J)

gamma0.it <- rnorm(J)



#use this to tell the sampler to restart from last run. Set to NA if first run
restart.file=NA
#if restarting, pull out the states to start from
if(!is.na(restart.file)){
  load(restart.file)
  stateList=vector("list",chains)
  for(i in 1:n.chain){
    stateList[[i]]=out[[i]]$stateList
  }
  #discard
  rm(out)
  gc()
}

#############################nimble model##############################################
n.chains = 3 #number of chains 
chains = vector("list", n.chains) #for full model
chains2 = vector("list", n.chains) #for sampling directly from the posterior (will be more heavily thinned)
out = for (chain in 1:n.chains){
  source("nimble_ppp_integrated_model_bear_10km_yearindexed_randomeffects_centrered.R")
  source("Restart MCMC Functions.R")
  
  # Parameters to be monitored
  parameters <-  c("beta0", "beta", "alpha0", "alpha", "zsum", "gamma0", "gamma", "beta0.it", 'alpha0.it')
  
  #data and constants
  Nimdata <- list(y = y.new, ones = rep(1, nrow(POdata))) #data goes here
  constants=list(npixel = G, conifmixed = conifmixed.new, deciduous = deciduous.new, edge = edge.new, POden = POcount.new, pa_pixel = PAdata$FID_allNY_10km,
                 po_pixel = POdata$FID_allNY_10km, opp_year = opp_year, npo = as.numeric(npo_count), m = nrow(POdata), K2D = K2D,
                 CONSTANT = constant, crops = crops.new, date = date.new, elevation = elevation.new,
                 nsurveys = dim(y.new)[2], prev = prev.new, J = J, nyears = n.year, cell_area = log(10)) #constants go here
  
  str(K2D)
  #inits
  Niminits=list(alpha=c(0,0,0), alpha0 = 0, beta=c(0,0,0,0,0), beta0 = runif(1), gamma0 = 0, gamma = 0, beta0.it = beta0.it, sd.lam = 1, gamma0.it = gamma0.it, sd.b = 1, alpha0.it = alpha0.it, sd.p = 1) #setting initial vales. When working on log or logit scale, 0 is usually a half decent starting value.
  
  #thinning rate
  nt = 3
  nt2 = 240
  # Build the model, configure the mcmc, and compile
  start.time<-Sys.time()
  Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                        inits=Niminits)
  
  conf <- configureMCMC(Rmodel,monitors=parameters, monitors2 = parameters, thin=nt, thin2 = nt2, useConjugacy = TRUE, enableWAIC=TRUE)
  
  # Build and compile
  Rmcmc <- buildMCMC(conf)
  Cmodel <- compileNimble(Rmodel)
  Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
  
  #if we are restarting from past states
  if(!is.na(restart.file)){
    ##restore mcmc and sampler states
    modelState <- stateList[[chain]]$modelState
    mcmcState <- stateList[[chain]]$mcmcState
    ## restore the saved "state" into the new model and new MCMC
    setModelState(Cmodel, modelState)
    setMCMCstate(conf, Cmcmc, mcmcState)
  }
  
  
    # n.iter = 4000
  n.iter = 30000
  # Run the model.
  start.time2<-Sys.time()
  Cmcmc$run(n.iter,reset=FALSE) #can keep running this line to get more sample
  end.time<-Sys.time()
  time1=end.time-start.time  # total time for compilation, replacing samplers, and fitting
  time2=end.time-start.time2 # post-com  
  
  
  #get the chains
    # n.iter = 2000
    # n.burn = 1000
  mvSamples= as.matrix(Cmcmc$mvSamples)
  mvSamples2 = as.matrix(Cmcmc$mvSamples2)
        # plot(mcmc(mvSamples[n.burn:nrow(mvSamples),]))
  chains[[chain]]=mvSamples
  chains2[[chain]]=mvSamples2
  
  #save MCMC and sampler states
  stateList <- list(modelState = getModelState(Cmodel),
                    mcmcState = getMCMCstate(conf, Cmcmc))
  
 allthestuff = list(chains=chains,chains2=chains2,
              stateList=stateList)
}

# save(allthestuff, file = "bear_SDM_elevation_decid_confimixed_crops_saved_chainsandstatelist.RData")
load("bear_SDM_elevation_decid_confimixed_crops_saved_chainsandstatelist.RData")
chains <- allthestuff$chains
chains2 <- allthestuff$chains2
# save(chains, file = "bear_SDM_converged_chains_yearindex.RData")
n.iter <- 10000 #number of iterations in thinned chains
n.burn <- 5000 #per chain burn in
a=mcmc.list(mcmc(chains[[1]][n.burn:n.iter,]),
            mcmc(chains[[2]][n.burn:n.iter,]),
            mcmc(chains[[3]][n.burn:n.iter,])) #make a list including all three chains

plot(a) #plot the traceplots
gelman.diag(a) #get gelman-rubin diagnostics
gelman <- gelman.diag(a) #get gelman-rubin diagnostics

str(gelman)
write.csv(gelman$psrf, 'bear_ISDM_gelman.csv')
a=runjags::combine.mcmc(a)    #combine the three chains    

#extract point estimates and credible intervals

mean(a[,"beta[1]"])
quantile(a[,"beta[1]"], probs = c(2.5, 97.5)/100)

mean(a[,"beta[2]"])
quantile(a[,"beta[2]"], probs = c(2.5, 97.5)/100)

mean(a[,"beta[3]"])
quantile(a[,"beta[3]"], probs = c(2.5, 97.5)/100)

mean(a[,"beta[4]"])
quantile(a[,"beta[4]"], probs = c(2.5, 97.5)/100)

mean(a[,"beta[5]"])
quantile(a[,"beta[5]"], probs = c(2.5, 97.5)/100)


mean(a[,"gamma[1]"])
quantile(a[,"gamma[1]"], probs = c(2.5, 97.5)/100)

mean(a[,"alpha[1]"])
quantile(a[,"alpha[1]"], probs = c(2.5, 97.5)/100)

mean(a[,"alpha[2]"])
quantile(a[,"alpha[2]"], probs = c(2.5, 97.5)/100)

mean(a[,"alpha[3]"])
quantile(a[,"alpha[3]"], probs = c(2.5, 97.5)/100)


# save(chains, file = "bear_SDM_converged_chains_yearindexrandomeff.RData")
# save(chains2, file = "bear_SDM_converged_chains_yearindexrandomeff2.RData")

