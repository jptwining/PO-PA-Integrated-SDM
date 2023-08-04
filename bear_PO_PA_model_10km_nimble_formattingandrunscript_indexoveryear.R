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
#remove NAs
# pixelcovdata$PO_den[(is.na(pixelcovdata$PO_den))] <- 0
# pixelcovdata$PO_den[pixelcovdata$PO_den>=10] <- 10
# pixelcovdata$Podenbear[(is.na(pixelcovdata$Podenbear))] <- 0
# pixelcovdata$Podenbear_trunc[(is.na(pixelcovdata$Podenbear_trunc))] <- 0

range(pixelcovdata$PO_join_count)
range(pixelcovdata$PO_join_count_truncated)

# #log transform
# log_po_den <- (log(pixelcovdata$PO_den+1))
# log_po_bearden <- log(pixelcovdata$Podenbear+1)
log_PO_count <- log(pixelcovdata$PO_join_count_truncated+1)
table(is.na(log_PO_count))
# plot(pixelcovdata$PO_join_count)
# plot(log_po_bearden)
# 
# hist(pixelcovdata$PO_join_count_truncated)

hist(log_PO_count)


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
setwd('C:/Users/jpt93/Documents/R/PO_PA model')
load("bear_SDM_elevation_decid_confimixed_crops_saved_chainsandstatelist.RData")
chains <- allthestuff$chains
chains2 <- allthestuff$chains2
##this is for sampling from the posterior all working withs chains2 (otherwise same as above)
n.iter <- 10000 #number of iterations in thinned chains
n.burn <- 9500#per chain burn in
b=mcmc.list(mcmc(chains[[1]][n.burn:n.iter,]),
            mcmc(chains[[2]][n.burn:n.iter,]),
            mcmc(chains[[3]][n.burn:n.iter,])) #make a list including all three chains

b=runjags::combine.mcmc(b)   #combine the three chains         

 plot(b)

 gelman.diag(b)
 # 
# plot(x = a[,"T1obs"], y = a[,"T1sim"])
# 
# Tobs <- as.data.frame(a[,'T1obs'])
# Tsim <- as.data.frame(a[,"T1sim"])
# 
# fittest <- cbind(Tobs, Tsim)
# colnames(fittest) <- c("Tobs", "Tsim")
# 
# range(fittest$Tobs)
# 
# reg <- lm(fittest$Tobs ~ fittest$Tsim)
# plot(x=fittest$Tobs, y = fittest$Tsim) + abline(0,1)
# 
# ggplot(fittest, aes(x = scale(Tobs), y = scale(Tsim))) + geom_point() +
#   theme_classic() + geom_abline(intercept = 0, slope = 1) + ylab("simulated") + xlab("observed")
#sampling from posterior to predict abundance and occupancy across the landscape

#colnames(mvSamples)           
cell_area = log(10) #set the cell area

n.samples = nrow(b) #numer of samples 
loglambda = matrix(NA, n.samples, nrow(pixelcovdata))   #sample directly from the posterior of the thinned chains
for (i in 1:n.samples){
  for (j in 1:nrow(pixelcovdata)){
    loglambda[i,j] = b[,'beta0'][[i]] + b[,'beta[1]'][[i]] * conifmixed.new[j] + b[,'beta[2]'][[i]] * deciduous.new[j] +
      b[,'beta[3]'][[i]] * edge.new[j] + b[,'beta[4]'][[i]] * elevation.new[j]+ b[,'beta[5]'][[i]] * crops.new[j] + cell_area
  }
}

lambda = matrix(NA, n.samples, nrow(pixelcovdata))   #create new matrix for lambda estimates

lambda=exp(loglambda) #exponent of log lambda to get lambda - THIS IS WHAT YOU WANT FOR EXPECTED ABUNDANCE

psi = matrix(NA, n.samples, nrow(pixelcovdata))   #create new matrix for psi estimates

psi= (1 - exp(-lambda)) #convert lambda to psi 

head(psi) #check it looks good

psi <- as.data.frame(psi) #make dataframe

lambda <- as.data.frame(lambda)
meanpsi <- colMeans(psi) #take a mean of all the iterations so you have a single point value
meanpsi <- as.data.frame(meanpsi) #make dataframe

lambda <- as.data.frame(lambda)

meanlambda <- colMeans(lambda)
meanlambda <- as.data.frame(meanlambda) #make dataframe

write.csv(meanpsi, "bear_integrated_final_SDM_mean_psi_converged.csv") #save it. 

write.csv(meanlambda, "bear_integrated_final_SDM_mean_lambda_converged.csv") #save it. 


hist(meanlambda)
########I THEN MAKE THE MAPS IN ARCMAP ###################


#marginal estimates of thinning parameter and PO covariate
range(POcount.new)
POden = seq(-0.5, 4, 0.02)

n.samples = nrow(b)
logit_b = matrix(NA, n.samples, length(POden))   
for (i in 1:n.samples){
  for (j in 1:length(POden)){
    logit_b[i,j] = b[,'gamma0'][[i]] + b[,'gamma[1]'][[i]] * POden[j] 
  }
}

b_POden <- plogis(logit_b)

POden_b_mean= colMeans(b_POden)
POdenCIs=apply(b_POden,2,quantile, c(0.025,0.975), na.rm=TRUE) #90% percentile intervals

POden_preds <- data.frame(mean = POden_b_mean, 
                          lower = POdenCIs[1,],
                          upper = POdenCIs[2,],
                          covariate = POden)

ggplot(POden_preds, aes(x = POden, y = mean)) + geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, linetype = "dashed")+ xlab("presence-only record density")+ylab("conditional probability of public detection")+
  geom_path(size=1)+
  theme_classic()

#marginal estimates of detection probability #######################################
range(date.new)
range(prev.new)
date = seq(-1.6, 2.5, 0.02)
prev0 = 0

# plot(date.new)
# hist(date)
n.samples = nrow(b)
logit_p_day_prev0 = matrix(NA, n.samples, length(date))   
for (i in 1:n.samples){
  for (j in 1:length(date)){
    logit_p_day_prev0[i,j] = b[,'alpha0'][[i]] + b[,'alpha[1]'][[i]] *  date[j] + b[,'alpha[2]'][[i]] *  pow(date[j], 2) + b[,'alpha[3]'][[i]] * prev0
  }
}



p_day_prev0 <- plogis(logit_p_day_prev0)

p_day_prev0_mean= colMeans(p_day_prev0)
p_day_prev0CIs=apply(p_day_prev0,2,quantile, c(0.025,0.975), na.rm=TRUE) #90% percentile intervals

p_day_prev0_preds <- data.frame(mean = p_day_prev0_mean, 
                               lower = p_day_prev0CIs[1,],
                               upper = p_day_prev0CIs[2,],
                               date = date,
                               prev = prev0)

ggplot(p_day_prev0_preds, aes(x = date, y = mean)) + geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, linetype = "dashed")+ xlab("presence-only record density")+ylab("thinning rate")+
  geom_path(size=1)+
  theme_classic()


prev1 = 1

# plot(date.new)
# hist(date)
n.samples = nrow(b)
logit_p_day_prev1 = matrix(NA, n.samples, length(date))   
for (i in 1:n.samples){
  for (j in 1:length(date)){
    logit_p_day_prev1[i,j] = b[,'alpha0'][[i]] + b[,'alpha[1]'][[i]] *  date[j] + b[,'alpha[2]'][[i]] *  pow(date[j], 2) + b[,'alpha[3]'][[i]] * prev1
  }
}



p_day_prev1 <- plogis(logit_p_day_prev1)

p_day_prev1_mean= colMeans(p_day_prev1)
p_day_prev1CIs=apply(p_day_prev1,2,quantile, c(0.025,0.975), na.rm=TRUE) #90% percentile intervals

p_day_prev1_preds <- data.frame(mean = p_day_prev1_mean, 
                                lower = p_day_prev1CIs[1,],
                                upper = p_day_prev1CIs[2,],
                                date = date,
                                prev = prev1)

ggplot(p_day_prev1_preds, aes(x = date, y = mean)) + geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, linetype = "dashed")+ xlab("presence-only record density")+ylab("thinning rate")+
  geom_path(size=1)+
  theme_classic()


p_day_prev_preds <- rbind(p_day_prev0_preds, p_day_prev1_preds)

mean(p_day_prev_preds$mean)
mean(p_day_prev_preds$lower)
mean(p_day_prev_preds$upper)
p_day_prev_preds$prev <- as.factor(p_day_prev_preds$prev)

ggplot(p_day_prev_preds, aes(x = date, y = mean, colour = prev)) + geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0, size = 0.75, linetype = "dashed")+ xlab("ordinal date")+ylab("detection probability")+
  geom_path(size=1.5)+
  theme_classic() + scale_colour_grey(name="learnt response", labels = c("0", "1")) + ylim(0,1)


#marginal occupancy elevation
#fix non focal parameters at their mean
conifmixed = 0
deciduous = 0 
edge = 0
elevation = seq(-1.5, 4.45, 0.02)
crops = 0
range(elevation.new)

length(elevation)


n.samples = nrow(b)
loglambda_elevation = matrix(NA, n.samples, length(elevation))   
for (i in 1:n.samples){
  for (j in 1:length(elevation)){
    loglambda_elevation[i,j] = b[,'beta0'][[i]] + b[,'beta[1]'][[i]] * conifmixed + b[,'beta[2]'][[i]] * deciduous +
      b[,'beta[3]'][[i]] * edge + b[,'beta[4]'][[i]] * elevation[j] +  b[,'beta[5]'][[i]] * crops + cell_area
  }
}



lambda_elevation <- exp(loglambda_elevation)

elevationlambdamean= colMeans(lambda_elevation)
elevationCIs=apply(lambda_elevation,2,quantile, c(0.05,0.95), na.rm=TRUE) #90% percentile intervals

elevation_preds <- data.frame(mean = elevationlambdamean, 
                              lower = elevationCIs[1,],
                              upper = elevationCIs[2,],
                              covariate = elevation)

ggplot(elevation_preds, aes(x = elevation, y = mean)) + geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, linetype = "dashed")+ xlab("elevation")+ylab("expected abundance")+
         geom_path(size=1)+
         theme_classic()
       
    ###############################################################################   
#marginal occupancy conifmixed
range(conifmixed.new)
plot(conifmixed.new)
#fix non focal parameters at their mean
conifmixed = seq(-1, 3, 0.02)
deciduous = 0 
edge = 0
elevation = 0
crops = 0

length(conifmixed)


n.samples = nrow(b)
loglambda_conifmixed = matrix(NA, n.samples, length(conifmixed))   
for (i in 1:n.samples){
  for (j in 1:length(conifmixed)){
    loglambda_conifmixed[i,j] = b[,'beta0'][[i]] + b[,'beta[1]'][[i]] * conifmixed[j] + b[,'beta[2]'][[i]] * deciduous +
      b[,'beta[3]'][[i]] * edge + b[,'beta[4]'][[i]] * elevation +  b[,'beta[5]'][[i]] * crops + cell_area
  }
}



lambda_conifmixed <- exp(loglambda_conifmixed)

conifmixedlambdamean= colMeans(lambda_conifmixed)
conifmixedCIs=apply(lambda_conifmixed,2,quantile, c(0.05,0.95), na.rm=TRUE) #90% percentile intervals

conifmixed_preds <- data.frame(mean = conifmixedlambdamean, 
                              lower = conifmixedCIs[1,],
                              upper = conifmixedCIs[2,],
                              covariate = conifmixed)

ggplot(conifmixed_preds, aes(x = conifmixed, y = mean)) + geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, linetype = "dashed")+ xlab("conifmixed")+ylab("expected abundance")+
  geom_path(size=1)+
  theme_classic()       
       
       
###############################################################################   
#marginal occupancy deciduous
range(deciduous.new)
plot(deciduous.new)
#fix non focal parameters at their mean
conifmixed = 0
deciduous = seq(-1, 1.5, 0.02)
edge = 0
elevation = 0
crops = 0

length(deciduous)


n.samples = nrow(b)
loglambda_deciduous = matrix(NA, n.samples, length(deciduous))   
for (i in 1:n.samples){
  for (j in 1:length(deciduous)){
    loglambda_deciduous[i,j] = b[,'beta0'][[i]] + b[,'beta[1]'][[i]] * conifmixed + b[,'beta[2]'][[i]] * deciduous[j] +
      b[,'beta[3]'][[i]] * edge + b[,'beta[4]'][[i]] * elevation +  b[,'beta[5]'][[i]] * crops + cell_area
  }
}



lambda_deciduous <- exp(loglambda_deciduous)

deciduouslambdamean= colMeans(lambda_deciduous)
deciduousCIs=apply(lambda_deciduous,2,quantile, c(0.05,0.95), na.rm=TRUE) #90% percentile intervals

deciduous_preds <- data.frame(mean = deciduouslambdamean, 
                               lower = deciduousCIs[1,],
                               upper = deciduousCIs[2,],
                              covariate = deciduous)

ggplot(deciduous_preds, aes(x = deciduous, y = mean)) + geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, linetype = "dashed")+ xlab("deciduous")+ylab("expected abundance")+
  geom_path(size=1)+
  theme_classic()             
       
       
       
###############################################################################   
#marginal occupancy edge
range(edge.new)
plot(edge.new)
#fix non focal parameters at their mean
conifmixed = 0
deciduous = 0
edge = seq(-1.9, 4, 0.02)
elevation = 0
crops = 0

length(edge)


n.samples = nrow(b)
loglambda_edge = matrix(NA, n.samples, length(edge))   
for (i in 1:n.samples){
  for (j in 1:length(edge)){
    loglambda_edge[i,j] = b[,'beta0'][[i]] + b[,'beta[1]'][[i]] * conifmixed + b[,'beta[2]'][[i]] * deciduous +
      b[,'beta[3]'][[i]] * edge[j] + b[,'beta[4]'][[i]] * elevation +  b[,'beta[5]'][[i]] * crops + cell_area
  }
}



lambda_edge <- exp(loglambda_edge)

edgelambdamean= colMeans(lambda_edge)
edgeCIs=apply(lambda_edge,2,quantile, c(0.05,0.95), na.rm=TRUE) #90% percentile intervals

edge_preds <- data.frame(mean = edgelambdamean, 
                              lower = edgeCIs[1,],
                              upper = edgeCIs[2,],
                         covariate = edge)

ggplot(edge_preds, aes(x = edge, y = mean)) + geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, linetype = "dashed")+ xlab("edge")+ylab("expected abundance")+
  geom_path(size=1)+
  theme_classic()             

###############################################################################   
#marginal occupancy crops
range(crops.new)
plot(crops.new)
#fix non focal parameters at their mean
conifmixed = 0
deciduous = 0
edge = 0
elevation = 0
crops = seq(-0.5, 4.5, 0.02)

length(crops)


n.samples = nrow(b)
loglambda_crops = matrix(NA, n.samples, length(crops))   
for (i in 1:n.samples){
  for (j in 1:length(crops)){
    loglambda_crops[i,j] = b[,'beta0'][[i]] + b[,'beta[1]'][[i]] * conifmixed + b[,'beta[2]'][[i]] * deciduous +
      b[,'beta[3]'][[i]] * edge + b[,'beta[4]'][[i]] * elevation +  b[,'beta[5]'][[i]] * crops[j] + cell_area
  }
}



lambda_crops <- exp(loglambda_crops)

cropslambdamean= colMeans(lambda_crops)
cropsCIs=apply(lambda_crops,2,quantile, c(0.05,0.95), na.rm=TRUE) #90% percentile intervals

crops_preds <- data.frame(mean = cropslambdamean, 
                         lower = cropsCIs[1,],
                         upper = cropsCIs[2,],
                         covariate = crops)

ggplot(crops_preds, aes(x = crops, y = mean)) + geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, linetype = "dashed")+ xlab("crops")+ylab("expected abundance")+
  geom_path(size=1)+
  theme_classic()             


x <- 4

all_cond_preds_df <- rbind(cbind(conifmixed_preds, Covariate = "coniferous and mixed forest", Value = scale(conifmixed_preds[,x])),
                           cbind(deciduous_preds, Covariate = "deciduous forest", Value = scale(deciduous_preds[,x])),
                           cbind(edge_preds, Covariate = "forest edge density", Value = scale(edge_preds[,x])),
                           cbind(elevation_preds, Covariate = "elevation", Value = scale(elevation_preds[,x])),
                           cbind(crops_preds, Covariate = "croplands", Value = scale(crops_preds[,x])))

all_cond_preds_df$Covariate <- as.factor(all_cond_preds_df$Covariate)

all_cond_preds_df$mean <- all_cond_preds_df$mean/5
all_cond_preds_df$lower <- all_cond_preds_df$lower/5
all_cond_preds_df$upper <- all_cond_preds_df$upper/5

#plot
ggplot(all_cond_preds_df, aes(x = Value, y = mean)) + geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, colour = "blue", fill = "blue", linetype = "dashed")+
  geom_path(size=1, colour="blue")+
  theme_classic()+
  facet_grid(all_cond_preds_df$Covariate)+
  theme(legend.position="none") + ylab(expression("marginal expected abundance ("~lambda~")")) + xlab("covariate value (scaled)")

       


#posterior predictive checks

dataNodes <- Rmodel$getNodeNames(dataOnly=TRUE)
parentNodes <- Rmodel$getParents(dataNodes, stochOnly = TRUE)  # `getParents` is new in nimble 0.11.0
simNodes <- Rmodel$getDependencies(parentNodes, self = FALSE)

# Build the model, configure the mcmc, and compile
model <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
            inits=Niminits)

cmodel  <- compileNimble(model)
mcmc    <- buildMCMC(cmodel, monitors = parentNodes)

cmcmc   <- compileNimble(mcmc, project = model)

samples <- runMCMC(cmcmc, niter = 1000, nburnin = 500)

nSamp <- nrow(samples)
n <- length(y)
ppSamples <- matrix(0, nSamp, n)

ppSamples <- matrix(0, nrow = nSamp, ncol =
                      length(Rmodel$expandNodeNames(dataNodes, returnScalarComponents = TRUE)))
postNames <- colnames(samples)

identical(colnames(samples),Rmodel$expandNodeNames(mcmc$mvSamples$getVarNames()))


ppSamplerNF <- nimbleFunction(
  setup = function(model, mcmc) {
    dataNodes <- model$getNodeNames(dataOnly = TRUE)
    parentNodes <- model$getParents(dataNodes, stochOnly = TRUE)
    cat("Stochastic parents of data are:", paste(parentNodes, collapse = ','), ".\n")
    simNodes <- model$getDependencies(parentNodes, self = FALSE)
    vars <- mcmc$mvSamples$getVarNames()  # need ordering of variables in mvSamples / samples matrix
    cat("Using posterior samples of:", paste(vars, collapse = ','), ".\n")
    n <- length(model$expandNodeNames(dataNodes, returnScalarComponents = TRUE))
  },
  run = function(samples = double(2)) {
    nSamp <- dim(samples)[1]
    ppSamples <- matrix(nrow = nSamp, ncol = n)
    for(i in 1:nSamp) {
      values(model, vars) <<- samples[i, ]
      model$simulate(simNodes, includeData = TRUE)
      ppSamples[i, ] <- values(model, dataNodes)
    }
    returnType(double(2))
    return(ppSamples)
  })

ppSampler <- ppSamplerNF(model, mcmc)

cppSampler <- compileNimble(ppSampler, project = model)

colnames(samples)

identical(colnames(samples), model$expandNodeNames(mcmc$mvSamples$getVarNames()))

system.time(ppSamples_via_nf <- cppSampler$run(samples))

# set.seed(1)
# system.time({
#   for(i in seq_len(nSamp)) {
#     values(cmodel, postNames) <- samples[i, ]  # assign 'flattened' values
#     cmodel$simulate(simNodes, includeData = TRUE)
#     ppSamples[i, ] <- values(cmodel, dataNodes)
#   }
# })
# 
# for(i in 1:nSamp){
#   cmodel[["beta0"]] <- samples[i, "beta0"]             
#   cmodel[["beta[1]"]] <- samples[i, "beta[1]"]
#   cmodel[["beta[2]"]] <- samples[i, "beta[2]"]
#   cmodel[["beta[3]"]] <- samples[i, "beta[3]"]
#   cmodel[["beta[4]"]] <- samples[i, "beta[4]"]
#   cmodel[["beta[5]"]] <- samples[i, "beta[5]"]
#   cmodel$simulate(simNodes, includeData = TRUE)
#   ppSamples[i, ] <- cmodel[["y"]]
# }
# 
#   cmodel[["alpha0"]] <- samples[i, "alpha0"]  
#   cmodel[["alpha[1]"]] <- samples[i, "alpha[1]"]  
#   cmodel[["alpha[2]"]] <- samples[i, "alpha[2]"]  
#   cmodel[["alpha[3]"]] <- samples[i, "alpha[3]"]  
#   cmodel[["zsum"]] <- samples[i, "zsum"]  
#   cmodel[["gamma0"]] <- samples[i, "gamma0"]  
#   cmodel[["gamma[1]"]] <- samples[i, "gamma[1]"]  
#   cmodel$simulate(simNodes, includeData = TRUE)
#   ppSamples[i, ] <- cmodel[["y"]]
# }


start.time<-Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)
conf <- configureMCMC(Rmodel,monitors=parameters, monitors2 = parameters, thin=nt, thin2 = nt2, useConjugacy = TRUE, enableWAIC=TRUE)

# Build and compile
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
