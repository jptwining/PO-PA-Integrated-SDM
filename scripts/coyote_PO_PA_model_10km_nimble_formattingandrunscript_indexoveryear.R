setwd('./data')
#read in detection non-detection data
PAdata <- read.csv("allNY_2013-2021_coyotePArecords_pixelID_10kmgrid.csv")

#read in presence-only data
POdata <- read.csv("allNY_2013_2021_coyotePOrecord_pixelID_10kmgrid_noweeklydups.csv")

#read in pixel covariate data
pixelcovdata <- read.csv("allNY_10kmgrid_alllandscapecovs_final.csv")

#read in observation cov data
date <- read.csv('juliandays_allNY_2013-2021.csv')


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
unique(POdata$Year)

unique(PAdata$year)
sitenames=unique(PAdata$sitename)
years=unique(PAdata$year)
n.year=length(years)
J=rep(NA,n.year)
year.sites=list(9)
for(i in 1:n.year){
  J[i]=length(PAdata$sitename[PAdata$year==years[i]])
  year.sites[[i]]=unique(PAdata$sitename[PAdata$year==years[i]])
}
maxJ=max(J)
totalJ=sum(J)
K=3


y <- PAdata[,(5:7)]

sum(y)

#3D matrix (site, occ, year)
y.new =array(0,dim=c(maxJ,K,n.year)) #new data 
for(i in 1:nrow(PAdata)){ #loop through each row
  this.year=which(years==PAdata$year[i]) #get year for this row
  this.site=which(year.sites[[this.year]]==PAdata$sitename[i]) #get site for this row
  y.new[this.site,1:3,this.year]=as.numeric(PAdata[i,5:7]) #force numeric
}


for(t in 1:n.year){
  print(c(sum(y.new[1:J[t],,t], na.rm=TRUE),sum(PAdata[PAdata$year==years[t],5:7], na.rm=TRUE)))
}

#create operational matrix
K2D <- 1*(!is.na(y.new))

#make NAs in observed data to 0
y.new[is.na(y.new)] <- 0

#if count data, set counts great than 1 to 1
y.new[y.new>=1]<-1



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


#matrix and standardize occ and behav covariate
occ <- as.matrix(occ)
prev <- as.matrix(prev)
occ.c = standardize(occ)

occyear <- cbind(occ.c, PAdata)
prevyear <- cbind(prev, PAdata)
#format occ covs
occ.new=array(0,dim=c(maxJ,K,n.year)) #new data 
for(i in 1:nrow(occ)){ #loop through each row
  this.year=which(years==occyear$year[i]) #get year for this row
  this.site=which(year.sites[[this.year]]==occyear$sitename[i]) #get site for this row
  occ.new[this.site,1:3,this.year]=as.numeric(occyear[i,1:3]) #force numeric
}

for(t in 1:n.year){
  print(c(sum(occ.new[1:J[t],,t]),sum(occyear[occyear$year==years[t],1:3])))
}


#format prev cov
prev.new=array(0,dim=c(maxJ,K,n.year)) #new data 
for(i in 1:nrow(occ)){ #loop through each row
  this.year=which(years==prevyear$year[i]) #get year for this row
  this.site=which(year.sites[[this.year]]==prevyear$sitename[i]) #get site for this row
  prev.new[this.site,1:3,this.year]=as.numeric(prevyear[i,1:3]) #force numeric
}

for(t in 1:n.year){
  print(c(sum(prev.new[1:J[t],,t]),sum(prevyear[prevyear$year==years[t],1:3])))
}

#dates
#read in ordinal day data, remove NAs and format
date <- read.csv('juliandays_allNY_2013-2021.csv')
date <- date[3:5]
range(date)
date <- as.matrix(date)
date.c = standardize(date)
dateyear <- cbind(date.c, PAdata)
# date[is.na(y)] <- NA # fix the NAs

#for week with year last
date.new=array(0,dim=c(maxJ,K,n.year)) #new data 
for(i in 1:nrow(PAdata)){ #loop through each row
  this.year=which(years==PAdata$year[i]) #get year for this row
  this.site=which(year.sites[[this.year]]==PAdata$sitename[i]) #get site for this row
  date.new[this.site,1:3,this.year]=as.numeric(dateyear[i,1:3]) #force numeric
}

for(t in 1:n.year){
  print(c(sum(date.new[1:J[t],,t]),sum(dateyear[dateyear$year==years[t],1:3])))
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
mydata <- pixelcovdata[,c(9 ,18, 20, 24,25, 30, 15)]
head(mydata)
mymodel <- lm(random~., data=mydata)
vif(mymodel) #no conlinearity 


# #random year level effect on lam inits
beta0.it  <- rnorm(J)

# #random year level effect on p and b inits
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
  source("./models/nimble_ppp_integrated_model_coyote_10km_yearindexed_randomeffects_centered.R")
  source("./models/Restart MCMC Functions.R")
  # Parameters monitored
  parameters <-  c("beta0", "beta", "alpha0", "alpha", "zsum", "gamma0", "gamma", "beta0.it", "alpha0.it", "gamma0.it")
  
  #data and constants
  Nimdata <- list(y = y.new, ones = rep(1, nrow(POdata))) #data goes here
  constants=list(npixel = G, conifmixed = conifmixed.new, deciduous = deciduous.new, edge = edge.new, POden = POcount.new, pa_pixel = PAdata$FID_allNY_10km,
                 po_pixel = POdata$FID_allNY_10km, opp_year = opp_year, npo = as.numeric(npo_count), m = nrow(POdata), K2D = K2D,
                 CONSTANT = constant, road = road.new, date = date.new, crops = crops.new,
                 nsurveys = dim(y.new)[2], occ = occ.new, J = J, nyears = n.year, cell_area = log(10)) #constants go here
  
  #inits
  Niminits=list(alpha=c(0,0,0), alpha0 = 0, beta=c(0,0,0,0,0), beta0 = runif(1), gamma0 = 0, gamma = 0, beta0.it = beta0.it, sd.lam = 1, gamma0.it = gamma0.it, sd.b = 1, alpha0.it = alpha0.it, sd.p = 1) #setting initial vales. When working on log or logit scale, 0 is usually a half decent starting value.
  
  #thinning rate
  nt = 3
  nt2 = 15
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

#set iter for thinned chains and burn in
n.iter <- 10000
n.burn < 5000
#combine the chains and burn
a=mcmc.list(mcmc(chains[[1]][n.burn:n.iter,]),
            mcmc(chains[[2]][n.burn:n.iter,]),
            mcmc(chains[[3]][n.burn:n.iter,]))

#plot the chains
plot(a)

#gelman-rubin diagnostics
gelman <- gelman.diag(a)

save(allthestuff, file = "coyote_SDM_saved_chainsandstatelist_randomsmonitorerd.RData")

