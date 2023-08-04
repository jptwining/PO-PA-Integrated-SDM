setwd('C:/Users/jpt93/Documents/R/PO_PA model')
load("coyote_SDM_saved_chainsandstatelist.RData")
chains <- allthestuff$chains
chains2 <- allthestuff$chains2
n.iter <- 10000
n.burn <- 5000
a=mcmc.list(mcmc(chains[[1]][n.burn:n.iter,]),
            mcmc(chains[[2]][n.burn:n.iter,]),
            mcmc(chains[[3]][n.burn:n.iter,]))

plot(a)

gelman.diag(a)
gelman <- gelman.diag(a) #get gelman-rubin diagnostics

str(gelman)
write.csv(gelman$psrf, 'coyote_ISDM_gelman.csv')
a=runjags::combine.mcmc(a)           

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



n.iter <- 2000
n.burn <- 1500
b=mcmc.list(mcmc(chains2[[1]][n.burn:n.iter,]),
            mcmc(chains2[[2]][n.burn:n.iter,]),
            mcmc(chains2[[3]][n.burn:n.iter,]))

gelman.diag(b)

plot(b)

b=runjags::combine.mcmc(b)           

#PPC fit

plot(x = a[,"T1obs"], y = a[,"T1sim"])

Tobs <- as.data.frame(a[,'T1obs'])
Tsim <- as.data.frame(a[,"T1sim"])

fittest <- cbind(Tobs, Tsim)
colnames(fittest) <- c("Tobs", "Tsim")

range(fittest$Tobs)

reg <- lm(fittest$Tobs ~ fittest$Tsim)
plot(x=fittest$Tobs, y = fittest$Tsim) + abline(0,1)

ggplot(fittest, aes(x = scale(Tobs), y = scale(Tsim))) + geom_point() +
  theme_classic() + geom_abline(intercept = 0, slope = 1) + ylab("simulated") + xlab("observed")



######sampling from posterior to predict abundance and occupancy across the landscape#########


#colnames(mvSamples)           
cell_area = log(10)
matrix(NA, 3, 4)

n.samples = nrow(b)
loglambda = matrix(NA, n.samples, nrow(pixelcovdata))   
for (i in 1:n.samples){
  for (j in 1:nrow(pixelcovdata)){
    loglambda[i,j] = b[,'beta0'][[i]] + b[,'beta[1]'][[i]] * conifmixed.new[j] + b[,'beta[2]'][[i]] * deciduous.new[j] +
      b[,'beta[3]'][[i]] * edge.new[j] + b[,'beta[4]'][[i]] * road.new[j]+ b[,'beta[5]'][[i]] * crops.new[j] + cell_area
  }
}

lambda = matrix(NA, n.samples, nrow(pixelcovdata))   

lambda=exp(loglambda)


lambda <- as.data.frame(lambda)

meanlambda <- colMeans(lambda)

write.csv(meanlambda, "coyote_integrated__final_SDM_mean_lambda.csv")

psi = matrix(NA, n.samples, nrow(pixelcovdata))   

psi= (1 - exp(-lambda))

head(psi)

psi <- as.data.frame(psi)

meanpsi <- colMeans(psi)
meanpsi <- as.data.frame(meanpsi)

write.csv(meanpsi, "coyote_integrated_final_SDM_mean_psi_02.csv")


#marginal occupancy for ROADS
#fix non focal parameters at their mean
conifmixed = 0
deciduous = 0 
edge = 0
road = seq(-1, 7.1, 0.02)
crops = 0
range(road.new)

length(road)


n.samples = nrow(b)
loglambda_road = matrix(NA, n.samples, length(road))   
for (i in 1:n.samples){
  for (j in 1:length(road)){
    loglambda_road[i,j] = b[,'beta0'][[i]] + b[,'beta[1]'][[i]] * conifmixed + b[,'beta[2]'][[i]] * deciduous +
      b[,'beta[3]'][[i]] * edge + b[,'beta[4]'][[i]] * road[j] +  b[,'beta[5]'][[i]] * crops + cell_area
  }
}


#marginal estimates of thinning parameter and PO covariate
range(POcount.new)
POden = seq(-0.5, 4.7, 0.02)

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
  theme_classic()+ylim(0,1)

#marginal estimates of detection probability #######################################
range(date.new)
range(occ.new)
date = seq(-2.2, 2.1, 0.02)
occ = -1.1224692

# plot(date.new)
# hist(date)
n.samples = nrow(b)
logit_p_day_occ1 = matrix(NA, n.samples, length(date))   
for (i in 1:n.samples){
  for (j in 1:length(date)){
    logit_p_day_occ1[i,j] = b[,'alpha0'][[i]] + b[,'alpha[1]'][[i]] * occ + b[,'alpha[2]'][[i]] * date[j] + b[,'alpha[3]'][[i]] * pow(date[j], 2)
  }
}



p_day_occ1 <- plogis(logit_p_day_occ1)

p_day_occ1_mean= colMeans(p_day_occ1)
p_day_occ1CIs=apply(p_day_occ1,2,quantile, c(0.025,0.975), na.rm=TRUE) #90% percentile intervals

p_day_occ1_preds <- data.frame(mean = p_day_occ1_mean, 
                               lower = p_day_occ1CIs[1,],
                               upper = p_day_occ1CIs[2,],
                               date = date,
                               occ = occ)

ggplot(p_day_occ1_preds, aes(x = date, y = mean)) + geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, linetype = "dashed")+ xlab("presence-only record density")+ylab("thinning rate")+
  geom_path(size=1)+
  theme_classic()

occ = 0


logit_p_day_occ2 = matrix(NA, n.samples, length(date))   
for (i in 1:n.samples){
  for (j in 1:length(date)){
    logit_p_day_occ2[i,j] = b[,'alpha0'][[i]] + b[,'alpha[1]'][[i]] * occ + b[,'alpha[2]'][[i]] * date[j] + b[,'alpha[3]'][[i]] * pow(date[j], 2)
  }
}


p_day_occ2 <- plogis(logit_p_day_occ2)

p_day_occ2_mean= colMeans(p_day_occ2)
p_day_occ2CIs=apply(p_day_occ2,2,quantile, c(0.025,0.975), na.rm=TRUE) #90% percentile intervals

p_day_occ2_preds <- data.frame(mean = p_day_occ2_mean, 
                               lower = p_day_occ2CIs[1,],
                               upper = p_day_occ2CIs[2,],
                               date = date,
                               occ = occ)


occ = 1.1224692

logit_p_day_occ3 = matrix(NA, n.samples, length(date))   
for (i in 1:n.samples){
  for (j in 1:length(date)){
    logit_p_day_occ3[i,j] = b[,'alpha0'][[i]] + b[,'alpha[1]'][[i]] * occ + b[,'alpha[2]'][[i]] * date[j] + b[,'alpha[3]'][[i]] * pow(date[j], 2)
  }
}


p_day_occ3 <- plogis(logit_p_day_occ3)

p_day_occ3_mean= colMeans(p_day_occ3)
p_day_occ3CIs=apply(p_day_occ3,2,quantile, c(0.025,0.975), na.rm=TRUE) #90% percentile intervals

p_day_occ3_preds <- data.frame(mean = p_day_occ3_mean, 
                               lower = p_day_occ3CIs[1,],
                               upper = p_day_occ3CIs[2,],
                               date = date,
                               occ = occ)


p_day_occ_preds <- rbind(p_day_occ1_preds, p_day_occ2_preds, p_day_occ3_preds)

p_day_occ_preds$occ <- as.factor(p_day_occ_preds$occ)

ggplot(p_day_occ_preds, aes(x = date, y = mean, colour = occ, fill = occ)) + geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0, size = 0.75, linetype = "dashed")+ xlab("ordinal date")+ylab("detection probability")+
  geom_path(size=1.5)+ scale_colour_grey(name="occasion (j)", labels = c("1", "2", "3"))+scale_fill_grey(name="occasion (j)", labels = c("1", "2", "3"))+
  theme_classic()

logit(p[j,k,t]) <- alpha0.it[t] + alpha[1] * occ[j,k,t] + alpha[2] * date[j,k,t] + alpha[3] *  pow(date[j,k,t], 2)


b_POden <- plogis(logit_b)

POden_b_mean= colMeans(b_POden)
POdenCIs=apply(b_POden,2,quantile, c(0.025,0.975), na.rm=TRUE) #90% percentile intervals

POden_preds <- data.frame(mean = POden_b_mean, 
                          lower = POdenCIs[1,],
                          upper = POdenCIs[2,],
                          covariate = POden)

ggplot(POden_preds, aes(x = POden, y = mean)) + geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, linetype = "dashed")+ xlab("presence-only record density")+ylab("thinning rate")+
  geom_path(size=1)+
  theme_classic()

range(POden_b_mean)


logit(b[pixel,t]) <-  gamma0.it[t] + gamma[1] * POden[pixel]



lambda_road <- exp(loglambda_road)

roadlambdamean= colMeans(lambda_road)
roadCIs=apply(lambda_road,2,quantile, c(0.05,0.95), na.rm=TRUE) #90% percentile intervals

road_preds <- data.frame(mean = roadlambdamean, 
                         lower = roadCIs[1,],
                         upper = roadCIs[2,],
                         covariate = road)

ggplot(road_preds, aes(x = road, y = mean)) + geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, linetype = "dashed")+ xlab("road")+ylab("expected abundance")+
  geom_path(size=1)+
  theme_classic()


#marginal occupancy for crops
range(crops.new)
#fix non focal parameters at their mean
conifmixed = 0
deciduous = 0 
edge = 0
crops = seq(-0.5, 5, 0.02)
road = 0
range(crops.new)

length(crops)


n.samples = nrow(b)
loglambda_crops = matrix(NA, n.samples, length(crops))   
for (i in 1:n.samples){
  for (j in 1:length(crops)){
    loglambda_crops[i,j] = b[,'beta0'][[i]] + b[,'beta[1]'][[i]] * conifmixed + b[,'beta[2]'][[i]] * deciduous +
      b[,'beta[3]'][[i]] * edge + b[,'beta[4]'][[i]] * road +  b[,'beta[5]'][[i]] * crops[j] + cell_area
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



#marginal occupancy for crops
range(conifmixed.new)
plot(conifmixed.new)
#fix non focal parameters at their mean
conifmixed = seq(0, 4, 0.02)
deciduous = 0 
edge = 0
crops = 0
road = 0

length(conifmixed)


n.samples = nrow(b)
loglambda_conifmixed = matrix(NA, n.samples, length(conifmixed))   
for (i in 1:n.samples){
  for (j in 1:length(conifmixed)){
    loglambda_conifmixed[i,j] = b[,'beta0'][[i]] + b[,'beta[1]'][[i]] * conifmixed[j] + b[,'beta[2]'][[i]] * deciduous +
      b[,'beta[3]'][[i]] * edge + b[,'beta[4]'][[i]] * road +  b[,'beta[5]'][[i]] * crops + cell_area
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



#marginal occupancy for decidous
range(deciduous.new)
plot(deciduous.new)
#fix non focal parameters at their mean
conifmixed = 0
deciduous = seq(0, 2.5, 0.02)
edge = 0
crops = 0
road = 0

length(deciduous)


n.samples = nrow(b)
loglambda_deciduous = matrix(NA, n.samples, length(deciduous))   
for (i in 1:n.samples){
  for (j in 1:length(deciduous)){
    loglambda_deciduous[i,j] = b[,'beta0'][[i]] + b[,'beta[1]'][[i]] * conifmixed + b[,'beta[2]'][[i]] * deciduous[j] +
      b[,'beta[3]'][[i]] * edge + b[,'beta[4]'][[i]] * road +  b[,'beta[5]'][[i]] * crops + cell_area
  }
}



lambda_deciduous <- exp(loglambda_deciduous)

deciduouslambdamean= colMeans(lambda_deciduous)
deciduousCIs=apply(lambda_deciduous,2,quantile, c(0.05,0.95), na.rm=TRUE) #90% percentile intervals

deciduous_preds <- data.frame(mean = deciduouslambdamean, 
                              lower = deciduousCIs[1,],
                              upper = deciduousCIs[2,],
                              covariate = deciduous)

ggplot(deciduous_preds, aes(x = deciduous, y = mean)) + geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, linetype = "dashed")+ xlab("deciduous forest")+ylab("expected abundance")+
  geom_path(size=1)+
  theme_classic()



#marginal occupancy for edge
range(edge.new)
plot(edge.new)
#fix non focal parameters at their mean
conifmixed = 0
deciduous = 0
edge = seq(0, 4, 0.02)
crops = 0
road = 0

length(edge)


n.samples = nrow(b)
loglambda_edge = matrix(NA, n.samples, length(edge))   
for (i in 1:n.samples){
  for (j in 1:length(edge)){
    loglambda_edge[i,j] = b[,'beta0'][[i]] + b[,'beta[1]'][[i]] * conifmixed + b[,'beta[2]'][[i]] * deciduous +
      b[,'beta[3]'][[i]] * edge[j] + b[,'beta[4]'][[i]] * road +  b[,'beta[5]'][[i]] * crops + cell_area
  }
}



lambda_edge <- exp(loglambda_edge)

edgelambdamean= colMeans(lambda_edge)
edgeCIs=apply(lambda_edge,2,quantile, c(0.05,0.95), na.rm=TRUE) #90% percentile intervals

edge_preds <- data.frame(mean = edgelambdamean, 
                         lower = edgeCIs[1,],
                         upper = edgeCIs[2,],
                         covariate = edge)

ggplot(edge_preds, aes(x = edge, y = mean)) + geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, linetype = "dashed")+ xlab("edge forest")+ylab("expected abundance")+
  geom_path(size=1)+
  theme_classic()



x <- 4

all_cond_preds_df <- rbind(cbind(conifmixed_preds, Covariate = "coniferous and mixed forest", Value = scale(conifmixed_preds[,x])),
                           cbind(deciduous_preds, Covariate = "deciduous forest", Value = scale(deciduous_preds[,x])),
                           cbind(edge_preds, Covariate = "forest edge density", Value = scale(edge_preds[,x])),
                           cbind(road_preds, Covariate = "road density", Value = scale(road_preds[,x])),
                           cbind(crops_preds, Covariate = "croplands", Value = scale(crops_preds[,x])))

all_cond_preds_df$Covariate <- as.factor(all_cond_preds_df$Covariate)


#plot
ggplot(all_cond_preds_df, aes(x = Value, y = mean, colour = "blue", fill = "blue")) + geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, linetype = "dashed")+
  geom_path(size=1)+
  theme_classic()+
  facet_grid(all_cond_preds_df$Covariate)+
  theme(legend.position="none") + ylab(expression("marginal expected abundance ("~lambda~")")) + xlab("covariate value (scaled)")




coefs <- read.csv('IDSM_allsp_coefs.csv')


coefs$overlap<- as.factor(coefs$overlap)

thinningcoefs <- subset(coefs, submodel == "thinning")
statecoefs <- subset(coefs, submodel == "state")
observationcoefs <- subset(coefs, submodel == "observation")

ggplot(data = statecoefs, aes(x = mean, y = covariate, group = species, colour = species, fill = species))+geom_point(aes(shape = species, colour = species, fill = species), size = 2.5, position=position_dodge(width=1)) +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0, size = 0.75, position=position_dodge(width=1))+ theme_classic()+ scale_colour_manual(values = c("7", '4','3')) + scale_fill_manual(values = c('7','4','2'))+
  ylab("landscape covariates") + xlab("coefficient estimates")  + xlim(-1,1.5)+ geom_vline(xintercept=0, linetype="dashed")# Make plot black and white with no background grid

ggplot(data = observationcoefs, aes(x = mean, y = covariate, group = species, colour = species, fill = species))+geom_point(aes(shape = species, colour = species, fill = species), size = 2.5, position=position_dodge(width=1)) +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0, size = 0.75, position=position_dodge(width=1))+ theme_classic()+ scale_colour_manual(values = c("7", '4','3')) + scale_fill_manual(values = c('7','4','2'))+
  ylab("observation covariates") + xlab("coefficient estimates")  + xlim(-1,2)+ geom_vline(xintercept=0, linetype="dashed")# Make plot black and white with no background grid

ggplot(data = thinningcoefs, aes(x = mean, y = covariate, group = species, colour = species, fill = species))+geom_point(aes(shape = species, colour = species, fill = species), size = 2.5, position=position_dodge(width=1)) +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0, size = 0.75, position=position_dodge(width=1))+ theme_classic()+ scale_colour_manual(values = c("7", '4','3')) + scale_fill_manual(values = c('7','4','2'))+
  ylab("thinning covariates") + xlab("coefficient estimates")  + xlim(-1,3.5)+ geom_vline(xintercept=0, linetype="dashed")# Make plot black and white with no background grid



