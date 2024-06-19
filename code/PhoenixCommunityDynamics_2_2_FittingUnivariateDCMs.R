rm(list=ls())
gc()


library(jagsUI)
library(rjags)
#library(runjags)
library(abind)

library(dplyr)
library(beepr)


setwd("~/GitHub/metrophxwildlife-communitydynamics")
set.seed(54321)

#### Import Dataset ####
  load("./data/model_inputs/ModelInputData_DCM.RData")
  
  
  # Three seasons of data from the Metro Phoenix Wildlife Study
  str(y)        # species observations
  str(ysum)     # species observations, summarized across survey occasions
  n.site        # number of sites = dim(y)[1]
  n.survey      # number of survey occasions = dim(y)[2]
  n.season      # number of seasons = dim(y)[3]
  n.spp         # number of species = dim(y)[4]
  str(K.survey) # number of occasions in which each site was sampled in each season
  names.common  # names of the species sampled
  
  # Site occupancy and detection covariates
  str(data.site)
  
    # standardized versions of seasonally-varying covariates
    # these were the only ones that were separate objects, because they were standardized by season
    str(K)            # number of days in which each site was surveyed
    str(K.std)        # K, but standardized by season
    str(ndvi.1km.std)
    # str(ndvihet.1km.std)
  
  # Species covariates
  str(data.spp)

  # Code for Calculating CPO #
  # this can be calculated in a line-efficient manner using a function, such as this one adapted from code provided by Mason Fidino
  calc_cpo <- function(posterior = NULL ,data = NULL){
    lik <- as.matrix(posterior$samples, chains = TRUE)         # $samples is for outputs from jagsUI, instead of $mcmc
    lik <- lik[,grep("lik", colnames(lik))]              # subset to only the parameters that have "lik" in their column names
    cpo_vec <- apply(1/lik, 2, sum)                      # sum the inverse of likelihood estimates across chains*samples values
    cpo_vec <- nrow(lik) / cpo_vec                       # harmonic mean of likelihood estimates across chains*samples values
    cpo <- -sum(log(cpo_vec), na.rm = TRUE)
    return(cpo)
  }
  # but, it's also only a few lines to calculate without the function
  
  # we will store these CPO summary statistic values in a dataframe
    # Note that, in earlier iterations of this analysis, we had two global
    # models, one including quadratic urbanization as a potential predictor term
    # In the manuscript, we presented the global model results without this 
    # quadratic, but we have kept it here, out of posterity and convenience.
  scales <- c(100, 1000, 5000)
  preds <- c("urbanization","urbanization quad", "heterogeneity", "vegetation greenness")
  n.scale <- length(scales)
  n.pred <- 4
  CPO <- data.frame(
    "predictor" = rep(preds, each = n.scale),
    "scale" = rep(scales, times = n.pred),
    "sumCPO" = rep(NA, times = n.scale*n.pred)
  )
  CPO
  
  
#### Setup for MCMC ####
  # Across all models, we will use these same settings and monitor the same parameters
  # MCMC settings FOR A SCALED DOWN VERSION OF THE MODEL
  na <- 100	   # pre-burnin
  nb <- 370	 # burn-in
  ni <- 400  # iterations (including burn-in with jagsUI)
  nt <- 3
  nc <- 3
  str(ysum)   # check that this is integer data (will not run otherwise)
  
  # Parameters monitored (could also add "z")
  params <- c(
    # community-level intercepts
      "mu.lpsi",
      "mu.lphi",
      "mu.lgamma",
    # community-level covariate effects
      "mu.beta1.psi", "mu.beta1.phi", "mu.beta1.gamma",
      "mu.beta2.psi", "mu.beta2.phi", "mu.beta2.gamma", 
      "mu.beta3.psi", "mu.beta3.phi", "mu.beta3.gamma",
      "mu.beta4.psi", "mu.beta4.phi", "mu.beta4.gamma",
      "mu.alpha1",
    # species trait covariate effects
      "dint.psi.trait1", "dint.psi.trait2",
      "dint.phi.trait1", "dint.phi.trait2",
      "dint.gamma.trait1", "dint.gamma.trait2",
      "deff1.psi.trait1", "deff1.psi.trait2",
      "deff1.phi.trait1", "deff1.phi.trait2",
      "deff1.gamma.trait1", "deff1.gamma.trait2",
      "deff2.psi.trait1", "deff2.psi.trait2",
      "deff2.phi.trait1", "deff2.phi.trait2",
      "deff2.gamma.trait1", "deff2.gamma.trait2",
      "deff3.psi.trait1", "deff3.psi.trait2",
      "deff3.phi.trait1", "deff3.phi.trait2",
      "deff3.gamma.trait1", "deff3.gamma.trait2",
      "deff4.psi.trait1", "deff4.psi.trait2",
      "deff4.phi.trait1", "deff4.phi.trait2",
      "deff4.gamma.trait1","deff4.gamma.trait2",
    # species-level covariate effects
      "beta1.psi", "beta1.phi", "beta1.gamma",
      "beta2.psi", "beta2.phi", "beta2.gamma",
      "beta3.psi", "beta3.phi", "beta3.gamma",
      "beta4.psi", "beta4.phi", "beta4.gamma",
      "alpha1",
    # Community Composition Metrics
      "rich",
      #"sum.psi",
      #"relative.psi",
      "hill1",
      #"H",
      #"E",
    # Species-level intercepts
      "lp",
      "lp.year",
      "lpsi1",
      "lphi",
      "lgamma",
    # Species-level occurrences
      #"n.occ",
      #"z",
      # Likelihood
      "lik"#,
    # Species-level Real Parameters
    #"p",
    #"gamma",
    #"psi1", 
    #"phi",
    #"psi"
  )
  
#### Fit Univariate Models ####
##### Urbanization at Different Scales #####
  ##### 100 m #####
    str(bdata <- list(
      y = ysum, 
      K = K.survey, 
      nsite = n.site, 
      nsurvey = n.survey,
      nseason = n.season,
      nspec = n.spp,
      cov.occ1 = data.site$imp100m_std,
      #cov.occ2 = (data.site$imp1km_std)^2,
      #cov.occ3 = data.site$sdhi1km_std,
      #cov.occ4 = as.matrix(ndvi.1km.std),
      trait1 = data.spp$logmass_std,
      trait2 = data.spp$dietdiv_std,
      det1 = as.matrix(K.std)
    )
    )
    
    # Initial values 
    zst <- apply(y, c(1,3,4), function(x) max(x, na.rm = T))  # Observed occurrence as inits for z
    inits <- function(){ list(
      z = zst
    )}
    
    # Call JAGS, check convergence, summarize posteriors
    # n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (start.time <- Sys.time())
    out <- jags(bdata, inits, params, 
                "./code/jags/model_DCM_urbonly.R", 
                n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (end.time <- Sys.time())
    elapsed.time <- difftime(end.time, start.time, units='mins')
    cat(paste(paste('Posterior computed in ', elapsed.time, sep=''), ' minutes', sep=''))
    
    #out$summary ; jags.View(out) ; print(out, 3) # not shown
    #par(mfrow = c(3,3)); traceplot(out, c("rich"))
    (CPO$sumCPO[1] <- calc_cpo(out))
    
    msum <- out$summary
    head(msum)
    
    beep(sound = "coin")
    
    # Export model outputs 
    #saveRDS(out, "./data/model_outputs/DCM_univariate/univariateDCM_urb100m.rds")
    #write.csv(msum, "./data/model_outputs/DCM_univariate/univariateDCM_urb100m_summary.csv")
    
  
  ##### 1 km #####
    bdata <- list(
      y = ysum, 
      K = K.survey, 
      nsite = n.site, 
      nsurvey = n.survey,
      nseason = n.season,
      nspec = n.spp,
      cov.occ1 = data.site$imp1km_std,   # this is all we have to change between urbanization models
      #cov.occ2 = (data.site$imp1km_std)^2,
      #cov.occ3 = data.site$sdhi1km_std,
      #cov.occ4 = as.matrix(ndvi.1km.std),
      trait1 = data.spp$logmass_std,
      trait2 = data.spp$dietdiv_std,
      det1 = as.matrix(K.std)
    )
    
    
    # Initial values 
    zst <- apply(y, c(1,3,4), function(x) max(x, na.rm = T))  # Observed occurrence as inits for z
    inits <- function(){ list(
      z = zst
    )}
    
    # Call JAGS, check convergence, summarize posteriors
    # n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (start.time <- Sys.time())
    out <- jags(bdata, inits, params, 
                "./code/jags/model_DCM_urbonly.R", 
                n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (end.time <- Sys.time())
    elapsed.time <- difftime(end.time, start.time, units='mins')
    cat(paste(paste('Posterior computed in ', elapsed.time, sep=''), ' minutes', sep=''))
    
    #out$summary ; jags.View(out) ; print(out, 3) # not shown
    #par(mfrow = c(3,3)); traceplot(out, c("rich"))
    
    # Add the calculated summary statistic to the appropriate row
    (CPO$sumCPO[2] <- calc_cpo(out))
    
    msum <- out$summary
    head(msum)
    
    beep(sound = "coin")
    
    # Export model outputs 
    #saveRDS(out, "./data/model_outputs/DCM_univariate/univariateDCM_urb1km.rds")
    #write.csv(msum, "./data/model_outputs/DCM_univariate/univariateDCM_urb1km_summary.csv")
    
    
  ##### 5 km #####
    bdata <- list(
      y = ysum, 
      K = K.survey, 
      nsite = n.site, 
      nsurvey = n.survey,
      nseason = n.season,
      nspec = n.spp,
      cov.occ1 = data.site$imp5km_std,
      #cov.occ2 = (data.site$imp1km_std)^2,
      #cov.occ3 = data.site$sdhi1km_std,
      #cov.occ4 = as.matrix(ndvi.1km.std),
      trait1 = data.spp$logmass_std,
      trait2 = data.spp$dietdiv_std,
      det1 = as.matrix(K.std)
    )
    
    
    # Initial values 
    zst <- apply(y, c(1,3,4), function(x) max(x, na.rm = T))  # Observed occurrence as inits for z
    inits <- function(){ list(
      z = zst
    )}
    
    # Call JAGS, check convergence, summarize posteriors
    # n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (start.time <- Sys.time())
    out <- jags(bdata, inits, params, 
                "./code/jags/model_DCM_urbonly.R", 
                n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (end.time <- Sys.time())
    elapsed.time <- difftime(end.time, start.time, units='mins')
    cat(paste(paste('Posterior computed in ', elapsed.time, sep=''), ' minutes', sep=''))
    
    #out$summary ; jags.View(out) ; print(out, 3) # not shown
    #par(mfrow = c(3,3)); traceplot(out, c("rich"))
    
    # Add the calculated summary statistic to the appropriate row
    (CPO$sumCPO[3] <- calc_cpo(out))
    
    msum <- out$summary
    head(msum)
    
    beep(sound = "coin")
    
    # Export model outputs 
    #saveRDS(out, "./data/model_outputs/DCM_univariate/univariateDCM_urb5km.rds")
    #write.csv(msum, "./data/model_outputs/DCM_univariate/univariateDCM_urb5km_summary.csv")
    
    
    
    
##### Urbanization with Quadratic at Different Scales ##### 
    ##### 100 m #####
    bdata <- list(
      y = ysum, 
      K = K.survey, 
      nsite = n.site, 
      nsurvey = n.survey,
      nseason = n.season,
      nspec = n.spp,
      cov.occ1 = data.site$imp1km_std,
      cov.occ2 = (data.site$imp100m_std)^2,
      #cov.occ3 = data.site$sdhi1km_std,
      #cov.occ4 = as.matrix(ndvi.1km.std),
      trait1 = data.spp$logmass_std,
      trait2 = data.spp$dietdiv_std,
      det1 = as.matrix(K.std)
    )
    
    
    # Initial values 
    zst <- apply(y, c(1,3,4), function(x) max(x, na.rm = T))  # Observed occurrence as inits for z
    inits <- function(){ list(
      z = zst
    )}
    
    # Call JAGS, check convergence, summarize posteriors
    # n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (start.time <- Sys.time())
    out <- jags(bdata, inits, params, 
                "./code/jags/model_DCM_urbquad.R", 
                n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (end.time <- Sys.time())
    elapsed.time <- difftime(end.time, start.time, units='mins')
    cat(paste(paste('Posterior computed in ', elapsed.time, sep=''), ' minutes', sep=''))
    
    #out$summary ; jags.View(out) ; print(out, 3) # not shown
    #par(mfrow = c(3,3)); traceplot(out, c("rich"))
    
    # Add the calculated summary statistic to the appropriate row
    (CPO$sumCPO[4] <- calc_cpo(out))
    
    msum <- out$summary
    head(msum)
    
    beep(sound = "coin")
    
    # Export model outputs 
    #saveRDS(out, "./data/model_outputs/DCM_univariate/univariateDCM_het1km.rds")
    #write.csv(msum, "./data/model_outputs/DCM_univariate/univariateDCM_impquad100m_summary.csv")
    
    
    ##### 1 km #####
    bdata <- list(
      y = ysum, 
      K = K.survey, 
      nsite = n.site, 
      nsurvey = n.survey,
      nseason = n.season,
      nspec = n.spp,
      cov.occ1 = data.site$imp1km_std,
      cov.occ2 = (data.site$imp1km_std)^2,
      #cov.occ3 = data.site$sdhi1km_std,
      #cov.occ4 = as.matrix(ndvi.1km.std),
      trait1 = data.spp$logmass_std,
      trait2 = data.spp$dietdiv_std,
      det1 = as.matrix(K.std)
    )
    
    
    # Initial values 
    zst <- apply(y, c(1,3,4), function(x) max(x, na.rm = T))  # Observed occurrence as inits for z
    inits <- function(){ list(
      z = zst
    )}
    
    # Call JAGS, check convergence, summarize posteriors
    # n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (start.time <- Sys.time())
    out <- jags(bdata, inits, params, 
                "./code/jags/model_DCM_urbquad.R", 
                n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (end.time <- Sys.time())
    elapsed.time <- difftime(end.time, start.time, units='mins')
    cat(paste(paste('Posterior computed in ', elapsed.time, sep=''), ' minutes', sep=''))
    
    #out$summary ; jags.View(out) ; print(out, 3) # not shown
    #par(mfrow = c(3,3)); traceplot(out, c("rich"))
    
    # Add the calculated summary statistic to the appropriate row
    (CPO$sumCPO[5] <- calc_cpo(out))
    
    msum <- out$summary
    head(msum)
    
    beep(sound = "coin")
    
    # Export model outputs 
    #saveRDS(out, "./data/model_outputs/DCM_univariate/univariateDCM_het1km.rds")
    #write.csv(msum, "./data/model_outputs/DCM_univariate/univariateDCM_impquad1km_summary.csv")
    
    ##### 5 km #####
    bdata <- list(
      y = ysum, 
      K = K.survey, 
      nsite = n.site, 
      nsurvey = n.survey,
      nseason = n.season,
      nspec = n.spp,
      cov.occ1 = data.site$imp1km_std,
      cov.occ2 = (data.site$imp5km_std)^2,
      #cov.occ3 = data.site$sdhi1km_std,
      #cov.occ4 = as.matrix(ndvi.1km.std),
      trait1 = data.spp$logmass_std,
      trait2 = data.spp$dietdiv_std,
      det1 = as.matrix(K.std)
    )
    
    
    # Initial values 
    zst <- apply(y, c(1,3,4), function(x) max(x, na.rm = T))  # Observed occurrence as inits for z
    inits <- function(){ list(
      z = zst
    )}
    
    # Call JAGS, check convergence, summarize posteriors
    # n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (start.time <- Sys.time())
    out <- jags(bdata, inits, params, 
                "./code/jags/model_DCM_urbquad.R", 
                n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (end.time <- Sys.time())
    elapsed.time <- difftime(end.time, start.time, units='mins')
    cat(paste(paste('Posterior computed in ', elapsed.time, sep=''), ' minutes', sep=''))
    
    #out$summary ; jags.View(out) ; print(out, 3) # not shown
    #par(mfrow = c(3,3)); traceplot(out, c("rich"))
    
    # Add the calculated summary statistic to the appropriate row
    (CPO$sumCPO[6] <- calc_cpo(out))
    
    msum <- out$summary
    head(msum)
    
    beep(sound = "coin")
    
    # Export model outputs 
    #saveRDS(out, "./data/model_outputs/DCM_univariate/univariateDCM_het1km.rds")
    #write.csv(msum, "./data/model_outputs/DCM_univariate/univariateDCM_impquad5km_summary.csv")
    
    
    
##### Heterogeneity at Different Scales #####
    ##### 100 m #####
    str(bdata <- list(
      y = ysum, 
      K = K.survey, 
      nsite = n.site, 
      nsurvey = n.survey,
      nseason = n.season,
      nspec = n.spp,
      #cov.occ1 = data.site$imp100m_std,    # now this doesn't matter
      #cov.occ2 = (data.site$imp1km_std)^2,
      cov.occ3 = data.site$sdhi100m_std,
      #cov.occ3 = data.site$pd100m_std,
      #cov.occ4 = as.matrix(ndvi.1km.std),
      trait1 = data.spp$logmass_std,
      trait2 = data.spp$dietdiv_std,
      det1 = as.matrix(K.std)
    )
    )
    
    # Initial values 
    zst <- apply(y, c(1,3,4), function(x) max(x, na.rm = T))  # Observed occurrence as inits for z
    inits <- function(){ list(
      z = zst
    )}
    
    # Call JAGS, check convergence, summarize posteriors
    # n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (start.time <- Sys.time())
    out <- jags(bdata, inits, params, 
                "./code/jags/model_DCM_hetonly.R", 
                n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (end.time <- Sys.time())
    elapsed.time <- difftime(end.time, start.time, units='mins')
    cat(paste(paste('Posterior computed in ', elapsed.time, sep=''), ' minutes', sep=''))
    
    #out$summary ; jags.View(out) ; print(out, 3) # not shown
    #par(mfrow = c(3,3)); traceplot(out, c("rich"))
    (CPO$sumCPO[7] <- calc_cpo(out))
    #2871.189 using SDHI
    #2876.209 using PD
    
    msum <- out$summary
    head(msum)
    
    beep(sound = "coin")
    
    # Export model outputs 
    #saveRDS(out, "./data/model_outputs/DCM_univariate/univariateDCM_het100m.rds")
    #write.csv(msum, "./data/model_outputs/DCM_univariate/univariateDCM_het100m_summary.csv")
    
    
    ##### 1 km #####
    bdata <- list(
      y = ysum, 
      K = K.survey, 
      nsite = n.site, 
      nsurvey = n.survey,
      nseason = n.season,
      nspec = n.spp,
      #cov.occ1 = data.site$imp1km_std,
      #cov.occ2 = (data.site$imp1km_std)^2,
      cov.occ3 = data.site$sdhi1km_std,
      #cov.occ3 = data.site$pd1km_std,
      #cov.occ4 = as.matrix(ndvi.1km.std),
      trait1 = data.spp$logmass_std,
      trait2 = data.spp$dietdiv_std,
      det1 = as.matrix(K.std)
    )
    
    
    # Initial values 
    zst <- apply(y, c(1,3,4), function(x) max(x, na.rm = T))  # Observed occurrence as inits for z
    inits <- function(){ list(
      z = zst
    )}
    
    # Call JAGS, check convergence, summarize posteriors
    # n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (start.time <- Sys.time())
    out <- jags(bdata, inits, params, 
                "./code/jags/model_DCM_hetonly.R", 
                n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (end.time <- Sys.time())
    elapsed.time <- difftime(end.time, start.time, units='mins')
    cat(paste(paste('Posterior computed in ', elapsed.time, sep=''), ' minutes', sep=''))
    
    #out$summary ; jags.View(out) ; print(out, 3) # not shown
    #par(mfrow = c(3,3)); traceplot(out, c("rich"))
    
    # Add the calculated summary statistic to the appropriate row
    (CPO$sumCPO[8] <- calc_cpo(out))
    # 2878.79 for SDHI
    # 2883.207 for PD
    
    msum <- out$summary
    head(msum)
    
    beep(sound = "coin")
    
    # Export model outputs 
    #saveRDS(out, "./data/model_outputs/DCM_univariate/univariateDCM_het1km.rds")
    #write.csv(msum, "./data/model_outputs/DCM_univariate/univariateDCM_het1km_summary.csv")
    
    ##### 5 km #####
    bdata <- list(
      y = ysum, 
      K = K.survey, 
      nsite = n.site, 
      nsurvey = n.survey,
      nseason = n.season,
      nspec = n.spp,
      #cov.occ1 = data.site$imp5km_std,  
      #cov.occ2 = (data.site$imp1km_std)^2,
      cov.occ3 = data.site$sdhi5km_std,
      #cov.occ3 = data.site$pd5km_std,
      #cov.occ4 = as.matrix(ndvi.1km.std),
      trait1 = data.spp$logmass_std,
      trait2 = data.spp$dietdiv_std,
      det1 = as.matrix(K.std)
    )
    
    
    # Initial values 
    zst <- apply(y, c(1,3,4), function(x) max(x, na.rm = T))  # Observed occurrence as inits for z
    inits <- function(){ list(
      z = zst
    )}
    
    # Call JAGS, check convergence, summarize posteriors
    # n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (start.time <- Sys.time())
    out <- jags(bdata, inits, params, 
                "./code/jags/model_DCM_hetonly.R", 
                n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (end.time <- Sys.time())
    elapsed.time <- difftime(end.time, start.time, units='mins')
    cat(paste(paste('Posterior computed in ', elapsed.time, sep=''), ' minutes', sep=''))
    
    #out$summary ; jags.View(out) ; print(out, 3) # not shown
    #par(mfrow = c(3,3)); traceplot(out, c("rich"))
    
    # Add the calculated summary statistic to the appropriate row
    (CPO$sumCPO[9] <- calc_cpo(out))
    # 2867.71 for SDHI
    #
    
    msum <- out$summary
    head(msum)
    
    beep(sound = "coin")
    
    # Export model outputs 
    #saveRDS(out, "./data/model_outputs/DCM_univariate/univariateDCM_het5km.rds")
    #write.csv(msum, "./data/model_outputs/DCM_univariate/univariateDCM_het5km_summary.csv")
  
  
    
##### Vegetation Greenness at Different Scales #####
    ##### 100 m #####
    bdata <- list(
      y = ysum, 
      K = K.survey, 
      nsite = n.site, 
      nsurvey = n.survey,
      nseason = n.season,
      nspec = n.spp,
      #cov.occ1 = data.site$imp1km_std,
      #cov.occ2 = (data.site$imp5km_std)^2,
      #cov.occ3 = data.site$sdhi1km_std,
      cov.occ4 = as.matrix(ndvi.100m.std),
      trait1 = data.spp$logmass_std,
      trait2 = data.spp$dietdiv_std,
      det1 = as.matrix(K.std)
    )
    
    
    # Initial values 
    zst <- apply(y, c(1,3,4), function(x) max(x, na.rm = T))  # Observed occurrence as inits for z
    inits <- function(){ list(
      z = zst
    )}
    
    # Call JAGS, check convergence, summarize posteriors
    # n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (start.time <- Sys.time())
    out <- jags(bdata, inits, params, 
                "./code/jags/model_DCM_vegonly.R", 
                n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (end.time <- Sys.time())
    elapsed.time <- difftime(end.time, start.time, units='mins')
    cat(paste(paste('Posterior computed in ', elapsed.time, sep=''), ' minutes', sep=''))
    
    #out$summary ; jags.View(out) ; print(out, 3) # not shown
    #par(mfrow = c(3,3)); traceplot(out, c("rich"))
    (CPO$sumCPO[10] <- calc_cpo(out))
    
    msum <- out$summary
    head(msum)
    
    beep(sound = "coin")
    
    # Export model outputs 
    #saveRDS(out, "./data/model_outputs/DCM_univariate/univariateDCM_veg100m.rds")
    #write.csv(msum, "./data/model_outputs/DCM_univariate/univariateDCM_veg100m_summary.csv")
    
    
    ##### 1 km #####
    bdata <- list(
      y = ysum, 
      K = K.survey, 
      nsite = n.site, 
      nsurvey = n.survey,
      nseason = n.season,
      nspec = n.spp,
      #cov.occ1 = data.site$imp1km_std,
      #cov.occ2 = (data.site$imp5km_std)^2,
      #cov.occ3 = data.site$sdhi1km_std,
      cov.occ4 = as.matrix(ndvi.1km.std),
      trait1 = data.spp$logmass_std,
      trait2 = data.spp$dietdiv_std,
      det1 = as.matrix(K.std)
    )
    
    
    # Initial values 
    zst <- apply(y, c(1,3,4), function(x) max(x, na.rm = T))  # Observed occurrence as inits for z
    inits <- function(){ list(
      z = zst
    )}
    
    # Call JAGS, check convergence, summarize posteriors
    # n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (start.time <- Sys.time())
    out <- jags(bdata, inits, params, 
                "./code/jags/model_DCM_vegonly.R", 
                n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (end.time <- Sys.time())
    elapsed.time <- difftime(end.time, start.time, units='mins')
    cat(paste(paste('Posterior computed in ', elapsed.time, sep=''), ' minutes', sep=''))
    
    #out$summary ; jags.View(out) ; print(out, 3) # not shown
    #par(mfrow = c(3,3)); traceplot(out, c("rich"))
    
    # Add the calculated summary statistic to the appropriate row
    (CPO$sumCPO[11] <- calc_cpo(out))
    
    msum <- out$summary
    head(msum)
    
    beep(sound = "coin")
    
    # Export model outputs 
    #saveRDS(out, "./data/model_outputs/DCM_univariate/univariateDCM_veg1km.rds")
    #write.csv(msum, "./data/model_outputs/DCM_univariate/univariateDCM_veg1km_summary.csv")
    
    ##### 5 km #####
    bdata <- list(
      y = ysum, 
      K = K.survey, 
      nsite = n.site, 
      nsurvey = n.survey,
      nseason = n.season,
      nspec = n.spp,
      #cov.occ1 = data.site$imp1km_std,
      #cov.occ2 = (data.site$imp5km_std)^2,
      #cov.occ3 = data.site$sdhi1km_std,
      cov.occ4 = as.matrix(ndvi.5km.std),
      trait1 = data.spp$logmass_std,
      trait2 = data.spp$dietdiv_std,
      det1 = as.matrix(K.std)
    )
    
    
    # Initial values 
    zst <- apply(y, c(1,3,4), function(x) max(x, na.rm = T))  # Observed occurrence as inits for z
    inits <- function(){ list(
      z = zst
    )}
    
    # Call JAGS, check convergence, summarize posteriors
    # n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (start.time <- Sys.time())
    out <- jags(bdata, inits, params, 
                "./code/jags/model_DCM_vegonly.R", 
                n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (end.time <- Sys.time())
    elapsed.time <- difftime(end.time, start.time, units='mins')
    cat(paste(paste('Posterior computed in ', elapsed.time, sep=''), ' minutes', sep=''))
    
    #out$summary ; jags.View(out) ; print(out, 3) # not shown
    #par(mfrow = c(3,3)); traceplot(out, c("rich"))
    
    # Add the calculated summary statistic to the appropriate row
    (CPO$sumCPO[12] <- calc_cpo(out))
    
    msum <- out$summary
    head(msum)
    
    beep(sound = "coin")
    
    # Export model outputs 
    #saveRDS(out, "./data/model_outputs/DCM_univariate/univariateDCM_veg5km.rds")
    #write.csv(msum, "./data/model_outputs/DCM_univariate/univariateDCM_veg5km_summary.csv")
    
    
    
  
    
#### Global Models ####
    CPO
    # top models
    # urbanization: 1 km
    # heterogeneity: 100 m 
    # veg: 1 km
    # quad urbanization: 5 km
    # MCMC settings
    na <- 1000	   # pre-burnin
    nb <- 3700	 # burn-in
    ni <- 4000  # iterations (including burn-in with jagsUI)
    nt <- 3
    nc <- 3
    
  ##### global model (1 km urb, 5 km het, 100 m veg #####
    cor(data.site$imp1km, (data.site$imp1km_std)^2) 
    cor(data.site$imp1km, data.site$sdhi5km)    # r = -0.51. That is rather high
    cor(data.site$sdhi5km, (data.site$imp1km_std)^2)    # -0.42 is okay
    cor(data.site$imp1km, data.site$sdhi100m)   # Let's try the next best model for heterogeneity: r = 0.03 is much better. 
    cor(data.site$sdhi100m, (data.site$imp1km_std)^2)   # but, now heterogeneity is more correlated with quadratic urb (-0.59)
    cor(data.site$sdhi1km, (data.site$imp1km_std)^2)
    
    cor(data.site$imp1km, data.site$ndvi_wd1km)
    cor(data.site$imp1km, data.site$ndvi_ww1km)
    cor(data.site$imp1km, data.site$ndvi_cw1km)   # only correlated in one season
    
    
    bdata <- list(
      y = ysum, 
      K = K.survey, 
      nsite = n.site, 
      nsurvey = n.survey,
      nseason = n.season,
      nspec = n.spp,
      cov.occ1 = data.site$imp1km_std,
      #cov.occ2 = (data.site$imp5km_std)^2,
      cov.occ3 = data.site$sdhi5km_std,
      cov.occ4 = as.matrix(ndvi.1km.std),
      trait1 = data.spp$logmass_std,
      trait2 = data.spp$dietdiv_std,
      det1 = as.matrix(K.std)
    )
    
    
    # Initial values 
    zst <- apply(y, c(1,3,4), function(x) max(x, na.rm = T))  # Observed occurrence as inits for z
    inits <- function(){ list(
      z = zst
    )}
    
    # Call JAGS, check convergence, summarize posteriors
    # n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (start.time <- Sys.time())
    out <- jags(bdata, inits, params, 
                "./code/jags/model_DCM_global.R", 
                n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (end.time <- Sys.time())
    elapsed.time <- difftime(end.time, start.time, units='mins')
    cat(paste(paste('Posterior computed in ', elapsed.time, sep=''), ' minutes', sep=''))
    
    #out$summary ; jags.View(out) ; print(out, 3) # not shown
    #par(mfrow = c(3,3)); traceplot(out, c("rich"))
    
    # Add the calculated summary statistic to the appropriate row
    #(CPO$sumCPO[11] <- calc_cpo(out))
    calc_cpo(out)
    # 2592.797
    
    msum <- out$summary
    head(msum)
    
    beep(sound = "coin")
    
    # Export model outputs 
    #saveRDS(out, "./data/model_outputs/DCM_univariate/univariateDCM_het1km.rds")
    #write.csv(msum, "./data/model_outputs/DCM_univariate/univariateDCM_impquad1km_summary.csv")
    
    
  ##### global model plus quadratic urbanization #####
    bdata <- list(
      y = ysum, 
      K = K.survey, 
      nsite = n.site, 
      nsurvey = n.survey,
      nseason = n.season,
      nspec = n.spp,
      cov.occ1 = data.site$imp1km_std,
      cov.occ2 = (data.site$imp1km_std)^2,
      cov.occ3 = data.site$sdhi5km_std,
      cov.occ4 = as.matrix(ndvi.1km.std),
      trait1 = data.spp$logmass_std,
      trait2 = data.spp$dietdiv_std,
      det1 = as.matrix(K.std)
    )
    
    
    # Initial values 
    zst <- apply(y, c(1,3,4), function(x) max(x, na.rm = T))  # Observed occurrence as inits for z
    inits <- function(){ list(
      z = zst
    )}
    
    # Call JAGS, check convergence, summarize posteriors
    # n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (start.time <- Sys.time())
    out <- jags(bdata, inits, params, 
                "./code/jags/model_DCM_global_plusquadurb.R", 
                n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (end.time <- Sys.time())
    elapsed.time <- difftime(end.time, start.time, units='mins')
    cat(paste(paste('Posterior computed in ', elapsed.time, sep=''), ' minutes', sep=''))
    
    #out$summary ; jags.View(out) ; print(out, 3) # not shown
    #par(mfrow = c(3,3)); traceplot(out, c("rich"))
    
    # Add the calculated summary statistic to the appropriate row
    #(CPO$sumCPO[11] <- calc_cpo(out))
    calc_cpo(out)
    # 2573.54  with 1 km quadratic urbanization (the same as the optimal scale for regular urb)
    # 2576.403 with 5 km quad urb (the optimal scale for the univariate model)
    
    msum <- out$summary
    head(msum)
    
    beep(sound = "coin")
    
    # Export model outputs 
    #saveRDS(out, "./data/model_outputs/DCM_univariate/univariateDCM_het1km.rds")
    #write.csv(msum, "./data/model_outputs/DCM_univariate/univariateDCM_impquad1km_summary.csv")  
    

    
#### Models to Examine Urban^2 vs. Heterogeneity
  ##### Urban + Urban^2 #####
    bdata <- list(
      y = ysum, 
      K = K.survey, 
      nsite = n.site, 
      nsurvey = n.survey,
      nseason = n.season,
      nspec = n.spp,
      cov.occ1 = data.site$imp1km_std,
      cov.occ2 = (data.site$imp1km_std)^2,
      #cov.occ3 = data.site$sdhi5km_std,
      #cov.occ4 = as.matrix(ndvi.1km.std),
      trait1 = data.spp$logmass_std,
      trait2 = data.spp$dietdiv_std,
      det1 = as.matrix(K.std)
    )
    
    
    # Initial values 
    zst <- apply(y, c(1,3,4), function(x) max(x, na.rm = T))  # Observed occurrence as inits for z
    inits <- function(){ list(
      z = zst
    )}
    
    # Call JAGS, check convergence, summarize posteriors
    # n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (start.time <- Sys.time())
    out <- jags(bdata, inits, params, 
                "./code/jags/model_DCM_urbhet.R", 
                n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (end.time <- Sys.time())
    elapsed.time <- difftime(end.time, start.time, units='mins')
    cat(paste(paste('Posterior computed in ', elapsed.time, sep=''), ' minutes', sep=''))
    
    #out$summary ; jags.View(out) ; print(out, 3) # not shown
    #par(mfrow = c(3,3)); traceplot(out, c("rich"))
    
    # Add the calculated summary statistic to the appropriate row
    #(CPO$sumCPO[11] <- calc_cpo(out))
    calc_cpo(out)
    # 2570.79
    
    msum <- out$summary
    head(msum, 10)
    # negative effects of urbanization and ALL urban^2 terms
    
    
    beep(sound = "coin")
    
  ##### Urban + Het #####
    bdata <- list(
      y = ysum, 
      K = K.survey, 
      nsite = n.site, 
      nsurvey = n.survey,
      nseason = n.season,
      nspec = n.spp,
      cov.occ1 = data.site$imp1km_std,
      cov.occ2 = data.site$sdhi5km_std,
      #cov.occ3 = data.site$sdhi5km_std,
      #cov.occ4 = as.matrix(ndvi.1km.std),
      trait1 = data.spp$logmass_std,
      trait2 = data.spp$dietdiv_std,
      det1 = as.matrix(K.std)
    )
    
    
    # Initial values 
    zst <- apply(y, c(1,3,4), function(x) max(x, na.rm = T))  # Observed occurrence as inits for z
    inits <- function(){ list(
      z = zst
    )}
    
    # Call JAGS, check convergence, summarize posteriors
    # n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (start.time <- Sys.time())
    out <- jags(bdata, inits, params, 
                "./code/jags/model_DCM_urbhet.R", 
                n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (end.time <- Sys.time())
    elapsed.time <- difftime(end.time, start.time, units='mins')
    cat(paste(paste('Posterior computed in ', elapsed.time, sep=''), ' minutes', sep=''))
    
    #out$summary ; jags.View(out) ; print(out, 3) # not shown
    #par(mfrow = c(3,3)); traceplot(out, c("rich"))
    
    # Add the calculated summary statistic to the appropriate row
    #(CPO$sumCPO[11] <- calc_cpo(out))
    calc_cpo(out)
    # 2585.533
    
    msum <- out$summary
    head(msum, 10)
    # negative effects of urbanization
    # negative effects of heterogeneity
    
    
    beep(sound = "coin")
    
    
  ##### Het only #####
    # From the 5 km SDHI model above
    # CPO stat = 2863.911
    # *positive* effects of 
    
  ##### Urban + Urban^2 + Het #####
    bdata <- list(
      y = ysum, 
      K = K.survey, 
      nsite = n.site, 
      nsurvey = n.survey,
      nseason = n.season,
      nspec = n.spp,
      cov.occ1 = data.site$imp1km_std,
      cov.occ2 = (data.site$imp1km_std)^2,
      cov.occ3 = data.site$sdhi5km_std,
      #cov.occ4 = as.matrix(ndvi.1km.std),
      trait1 = data.spp$logmass_std,
      trait2 = data.spp$dietdiv_std,
      det1 = as.matrix(K.std)
    )
    
    
    # Initial values 
    zst <- apply(y, c(1,3,4), function(x) max(x, na.rm = T))  # Observed occurrence as inits for z
    inits <- function(){ list(
      z = zst
    )}
    
    # Call JAGS, check convergence, summarize posteriors
    # n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (start.time <- Sys.time())
    out <- jags(bdata, inits, params, 
                "./code/jags/model_DCM_urbquadhet.R", 
                n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (end.time <- Sys.time())
    elapsed.time <- difftime(end.time, start.time, units='mins')
    cat(paste(paste('Posterior computed in ', elapsed.time, sep=''), ' minutes', sep=''))
    
    #out$summary ; jags.View(out) ; print(out, 3) # not shown
    #par(mfrow = c(3,3)); traceplot(out, c("rich"))
    
    # Add the calculated summary statistic to the appropriate row
    #(CPO$sumCPO[11] <- calc_cpo(out))
    calc_cpo(out)
    # 2585.533
    
    msum <- out$summary
    head(msum, 11)
    # negative effects of urbanization
    # negative effects of heterogeneity
    
    
    beep(sound = "coin")
#### Export the model selection criteria #### 
    CPO
    #write.csv(CPO, "./data/model_outputs/DCM_univariate/univariate_modelselectionCPO.csv", row.names = FALSE)
    
    