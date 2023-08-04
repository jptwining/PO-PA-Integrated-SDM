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


