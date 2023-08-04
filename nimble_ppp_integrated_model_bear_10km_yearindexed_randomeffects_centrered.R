library(nimble)

NimModel <- nimbleCode({
  
  # Nimble version of the Fidino (2021) / Koshkina et al. 2017 integrated species distribution model.
  #
  # --- Priors ---
  #Priors and hyperpriors for state model
  beta0 ~ dnorm(0, 0.01)
  for (k in 1:5){
    beta[k] ~ dnorm(0, 0.01)  
  }
  tau.lam <- pow(sd.lam, -2)
  sd.lam ~ dunif(0, 2)
  #
  #Priors and hyperpriors  for observation model
  alpha0 <- logit(mean.p)
  mean.p ~ dunif(0, 1)
  for (k in 1:3){
    alpha[k] ~ dnorm(0, 0.1)
  }
  tau.p <- pow(sd.p, -2)
  sd.p ~ dunif(0, 2)
  gamma0 ~ dnorm(0, 0.01)
  gamma[1] ~ dnorm(0, 0.01)  
  tau.b <- pow(sd.b, -2)
  sd.b ~ dunif(0, 2)
  
  # # #year level random effect on PO thinning
  for (t in 1:nyears){ #loop over years
    gamma0.it[t] ~ dnorm(gamma0, tau.b)
  }
  # # #year level random effect on lambda
  for (t in 1:nyears){ #loop over years
    beta0.it[t] ~ dnorm(beta0, tau.lam)
  }
  # # #year level random effect on PA detection
  for (t in 1:nyears){ #loop over years
    alpha0.it[t] ~ dnorm(alpha0, tau.p)
  }
  
  # The latent-state model
  for(pixel in 1:npixel){
    for (t in 1:nyears){    #
    # beta = latent state model regression coefficients, abundance covariates go on line 33
    # cell_area = log area of grid cell 
    #eps.lam = year level random effect on lambda
      
    # Species presence in a gridcell as a Bernoulli trial
    z[pixel,t] ~ dbern(1 - exp(-lambda[pixel,t]))
    log(lambda[pixel,t]) <- beta0.it[t] + beta[1] * conifmixed[pixel] + beta[2] * deciduous[pixel] + beta[3] * edge[pixel] + beta[4] * elevation[pixel] + beta[5] * crops[pixel] + cell_area
    
    #Thinning poisson process for observing a species in PO data
    logit(b[pixel,t]) <-  gamma0.it[t] + gamma[1] * POden[pixel]
    }
  }
  for (t in 1:nyears){
  # The presence_only data model denominator, which is the thinned poisson process across the whole region
  po_denominator[t] <- inprod(lambda[1:npixel, t], b[1:npixel, t] ) / npo[t]
  }
  #
  # Loop through each presence-only data point  using Bernoulli one's trick
  for(po in 1:m){
    ones[po] ~ dbern(
      exp(log(lambda[po_pixel[po], opp_year[po]]*b[po_pixel[po], opp_year[po]]) - log(po_denominator[opp_year[po]])) / CONSTANT)
  } 
  #
  #Observation model for replicated counts, observation covariates go on line 52
      for (t in 1:nyears) {
        for (j in 1:J[t]){
          for (k in 1:nsurveys) {
      y[j, k, t] ~ dbern((p[j,k,t] * K2D[j,k,t]) * z[pa_pixel[j], t])
      logit(p[j,k,t]) <- alpha0.it[t] + alpha[1] * date[j,k,t] + alpha[2] * pow(date[j,k,t], 2)+ alpha[3] * prev[j,k,t]
    }
    }
  }
  for (t in 1:nyears) {
  # Derived parameter, the number of cells occupied
  zsum[t] <- sum(z[1:npixel,t])
  }
}
)



