rm(list=ls())
gc()


library(jagsUI)
library(rjags)
#library(runjags)
library(abind)

library(dplyr)
library(beepr)
library(ggplot2)
library(ggcorrplot)


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
  str(data.site.long) # the same data in a long format
  
    # standardized versions of seasonally-varying covariates
    # these were the only ones that were separate objects, because they were standardized by season
    str(K)            # number of days in which each site was surveyed
    str(K.std)        # K, but standardized by season
    str(ndvi.1km.std)
    # str(ndvihet.1km.std)
  
  # Species covariates
  str(data.spp)
  dimnames(ysum)[[3]]  # the data is sorted by 'name_short'
  data.spp <- data.spp %>% arrange(name_short) # sort the species info to match

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

  
  # select species traits to be used as covariates
  traits.response <- data.spp %>% arrange(name_short) %>%
    # select which covariates to model
    select(
      logmass_std,
      dietdiv_std
    ) %>% as.matrix(nrow = nrow())  
  
  # vectors of taxonomic group membership, for controlling for phylogenetic effects
  order_vec <- data.spp %>% arrange(name_short) %>% pull(order_vec)
  class_vec <- data.spp %>% arrange(name_short) %>% pull(class_vec)
  
  
#### Setup for MCMC ####
  # Across all models, we will use these same settings and monitor the same parameters
  # MCMC settings FOR A SCALED DOWN VERSION OF THE MODEL
  na <- 1000	   # pre-burnin
  nb <- 3700	 # burn-in
  ni <- 4000  # iterations (including burn-in with jagsUI)
  nt <- 3
  nc <- 3
  str(ysum)   # check that this is integer data (will not run otherwise)
  
  # Parameters monitored (could also add "z")
  params <- c(
    # community-level intercepts and slopes
      "beta.comm.psi1",
      "beta.comm.phi",
      "beta.comm.gamma",
      "rho.comm",
    # species trait covariate effects
      "eff.traits.psi",
      "eff.traits.phi",
      "eff.traits.gamma",
      "eff.class.psi",
      "eff.order.psi",
    # species-level intercepts and slopes
      "beta.sp.psi1",
      "beta.sp.phi", 
      "beta.sp.gamma",
      "rho.sp.tot",
      # "rho.sp",
    # Community Composition Metrics
      # "rich",
      #"sum.psi",
      #"relative.psi",
      # "hill1",
      #"H",
      #"E",
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
  # select the covariates that will be included in the model here
  scale.model <- "100m"
  covariates.occ <- c(
    "intercept",
    "imp_std"#,
    # "impquad",
    # "sdhi_std",
    # "ndvi_std"
  )
  covariates.det <- c(
    "intercept",
    "effort"
  )
  
  # Organize the covariates as matrices/arrays that can be easily subset as vectors within the model
  # these vectors will be multiplied by the slope parameters within the regression term 
  # the last dimension is the number of covariates, including the intercept
  a.cov.occ <- array(NA, dim = c(n.site, n.season, length(covariates.occ))) 
  a.cov.det <- array(NA, dim = c(n.site, n.season, length(covariates.det)))
  # The first "covariate" will need to be multiplied the intercept, so let's start with '1' for that
  a.cov.occ[,,1] <- 1
  a.cov.det[,,1] <- 1
  
  # Add the covariates that are being included  
  for(s in 1:n.season){
    a.cov.occ[,s,-1] <- data.site.long %>% # for each covariate past the first (the intercept)...
      filter(scale == scale.model) %>%     # select the scale specified
      filter(season == unique(data.site.long$season)[s]) %>%  # select the season
      select(covariates.occ[-1]) %>% as.matrix()          # select the covariates, then add to the array of covariate values
  }
  
  a.cov.det[,,2] <- as.matrix(K.std)
  
    # str(
      bdata <- list(
      y = ysum, 
      K = K.survey, 
      nsite = n.site, 
      nsurvey = n.survey,
      nseason = n.season,
      nspec = n.spp,
      psi.cov = a.cov.occ,
      ncov.occ = dim(a.cov.occ)[[3]],
      rho.cov = a.cov.det,
      ncov.det = dim(a.cov.det)[[3]],
      traits = traits.response,
      n.traits = dim(traits.response)[[2]],
      class_vec = class_vec,
      n.classes = length(unique(class_vec)),
      order_vec = order_vec,
      n.orders = length(unique(order_vec)))
    # )
    
  
  
    # Initial values 
    zst <- apply(y, c(1,3,4), function(x) max(x, na.rm = T))  # Observed occurrence as inits for z
    inits <- function(){ list(
      z = zst
    )}
    
    # Call JAGS, check convergence, summarize posteriors
    # n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (start.time <- Sys.time())
    out <- jags(bdata, inits, params, 
                "./code/jags/model_DCM.R", 
                n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (end.time <- Sys.time())
    elapsed.time <- difftime(end.time, start.time, units='mins')
    cat(paste(paste('Posterior computed in ', elapsed.time, sep=''), ' minutes', sep=''))
    
    #out$summary ; jags.View(out) ; print(out, 3) # not shown
    #par(mfrow = c(3,3)); traceplot(out, c("rich"))
    (CPO$sumCPO[1] <- calc_cpo(out))
    
    # drop 'lik' from the samples, to reduce model object size before exporting 
    # now that the summary stat has been calculated, that parameter isn't necessary and takes up a lot of memory
    for(i in 1:nc){
      lik <- as.matrix(out$samples[i])
      not.lik <- lik[,-c(grep("lik", colnames(lik)))]  # drop the samples of the 'lik' parameter
      out$samples[[i]] <- not.lik
    }
    
    msum <- out$summary
    head(msum)
    
    beep(sound = "coin")
    
    # Export model outputs 
    #saveRDS(out, "./data/model_outputs/DCM_univariate/univariateDCM_urb100m.rds")
    # write.csv(msum, "./data/model_summaries/DCM_univariate/univariateDCM_urb100m_summary.csv")
    
  
  ##### 1 km #####
    # repeat from above
    # select the covariates that will be included in the model here
    scale.model <- "1km"
    covariates.occ <- c(
      "intercept",
      "imp_std"#,
      # "impquad",
      # "sdhi_std",
      # "ndvi_std"
    )
    covariates.det <- c(
      "intercept",
      "effort"
    )
    
    a.cov.occ <- array(NA, dim = c(n.site, n.season, length(covariates.occ))) # organize the covariates into arrays
    a.cov.det <- array(NA, dim = c(n.site, n.season, length(covariates.det)))
    a.cov.occ[,,1] <- 1   # set the intercepts to 1
    a.cov.det[,,1] <- 1
      
    for(s in 1:n.season){    # set the occupancy covariate values from each season
      a.cov.occ[,s,-1] <- data.site.long %>% # for each covariate past the first (the intercept)...
        filter(scale == scale.model) %>%     # select the scale specified
        filter(season == unique(data.site.long$season)[s]) %>%  # select the season
        select(covariates.occ[-1]) %>% as.matrix()          # select the covariates, then add to the array of covariate values
    }
    
    a.cov.det[,,2] <- as.matrix(K.std)   # set the detection covariates
    
    # str(
      bdata <- list(
      y = ysum, 
      K = K.survey, 
      nsite = n.site, 
      nsurvey = n.survey,
      nseason = n.season,
      nspec = n.spp,
      psi.cov = a.cov.occ,
      ncov.occ = dim(a.cov.occ)[[3]],
      rho.cov = a.cov.det,
      ncov.det = dim(a.cov.det)[[3]],
      traits = traits.response,
      n.traits = dim(traits.response)[[2]],
      class_vec = class_vec,
      n.classes = length(unique(class_vec)),
      order_vec = order_vec,
      n.orders = length(unique(order_vec)))
    # )
    
    
    # Initial values 
    zst <- apply(y, c(1,3,4), function(x) max(x, na.rm = T))  # Observed occurrence as inits for z
    inits <- function(){ list(
      z = zst
    )}
    
    # Call JAGS, check convergence, summarize posteriors
    # n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (start.time <- Sys.time())
    out <- jags(bdata, inits, params, 
                "./code/jags/model_DCM.R", 
                n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (end.time <- Sys.time())
    elapsed.time <- difftime(end.time, start.time, units='mins')
    cat(paste(paste('Posterior computed in ', elapsed.time, sep=''), ' minutes', sep=''))
    
    #out$summary ; jags.View(out) ; print(out, 3) # not shown
    #par(mfrow = c(3,3)); traceplot(out, c("rich"))
    
    # Add the calculated summary statistic to the appropriate row
    (CPO$sumCPO[2] <- calc_cpo(out))
    
    for(i in 1:nc){
      lik <- as.matrix(out$samples[i])
      not.lik <- lik[,-c(grep("lik", colnames(lik)))]  # drop the samples of the 'lik' parameter
      out$samples[[i]] <- not.lik
    }
    
    msum <- out$summary
    head(msum)
    
    beep(sound = "coin")
    
    # Export model outputs 
    #saveRDS(out, "./data/model_outputs/DCM_univariate/univariateDCM_urb1km.rds")
    # write.csv(msum, "./data/model_summaries/DCM_univariate/univariateDCM_urb1km_summary.csv")
    
    
  ##### 5 km #####
    # repeat from above
    # select the covariates that will be included in the model here
    scale.model <- "5km"
    covariates.occ <- c(
      "intercept",
      "imp_std"#,
      # "impquad",
      # "sdhi_std",
      # "ndvi_std"
    )
    covariates.det <- c(
      "intercept",
      "effort"
    )
    
    a.cov.occ <- array(NA, dim = c(n.site, n.season, length(covariates.occ))) # organize the covariates into arrays
    a.cov.det <- array(NA, dim = c(n.site, n.season, length(covariates.det)))
    a.cov.occ[,,1] <- 1   # set the intercepts to 1
    a.cov.det[,,1] <- 1
    
    for(s in 1:n.season){    # set the occupancy covariate values from each season
      a.cov.occ[,s,-1] <- data.site.long %>% # for each covariate past the first (the intercept)...
        filter(scale == scale.model) %>%     # select the scale specified
        filter(season == unique(data.site.long$season)[s]) %>%  # select the season
        select(covariates.occ[-1]) %>% as.matrix()          # select the covariates, then add to the array of covariate values
    }
    
    a.cov.det[,,2] <- as.matrix(K.std)   # set the detection covariates
    
    # str(
      bdata <- list(
      y = ysum, 
      K = K.survey, 
      nsite = n.site, 
      nsurvey = n.survey,
      nseason = n.season,
      nspec = n.spp,
      psi.cov = a.cov.occ,
      ncov.occ = dim(a.cov.occ)[[3]],
      rho.cov = a.cov.det,
      ncov.det = dim(a.cov.det)[[3]],
      traits = traits.response,
      n.traits = dim(traits.response)[[2]],
      class_vec = class_vec,
      n.classes = length(unique(class_vec)),
      order_vec = order_vec,
      n.orders = length(unique(order_vec)))
    # )
    
    
    # Initial values 
    zst <- apply(y, c(1,3,4), function(x) max(x, na.rm = T))  # Observed occurrence as inits for z
    inits <- function(){ list(
      z = zst
    )}
    
    # Call JAGS, check convergence, summarize posteriors
    # n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (start.time <- Sys.time())
    out <- jags(bdata, inits, params, 
                "./code/jags/model_DCM.R", 
                n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (end.time <- Sys.time())
    elapsed.time <- difftime(end.time, start.time, units='mins')
    cat(paste(paste('Posterior computed in ', elapsed.time, sep=''), ' minutes', sep=''))
    
    #out$summary ; jags.View(out) ; print(out, 3) # not shown
    #par(mfrow = c(3,3)); traceplot(out, c("rich"))
    
    # Add the calculated summary statistic to the appropriate row
    (CPO$sumCPO[3] <- calc_cpo(out))
    
    for(i in 1:nc){
      lik <- as.matrix(out$samples[i])
      not.lik <- lik[,-c(grep("lik", colnames(lik)))]  # drop the samples of the 'lik' parameter
      out$samples[[i]] <- not.lik
    }
    
    msum <- out$summary
    head(msum)
    
    beep(sound = "coin")
    
    # Export model outputs 
    #saveRDS(out, "./data/model_outputs/DCM_univariate/univariateDCM_urb5km.rds")
    # write.csv(msum, "./data/model_summaries/DCM_univariate/univariateDCM_urb5km_summary.csv")
    
    
    
    
##### Urbanization with Quadratic at Different Scales ##### 
    ##### 100 m #####
    # repeat from above
    # select the covariates that will be included in the model here
    scale.model <- "100m"
    covariates.occ <- c(
      "intercept",
      "imp_std",
      "impquad"#,
      # "sdhi_std",
      # "ndvi_std"
    )
    covariates.det <- c(
      "intercept",
      "effort"
    )
    
    a.cov.occ <- array(NA, dim = c(n.site, n.season, length(covariates.occ))) # organize the covariates into arrays
    a.cov.det <- array(NA, dim = c(n.site, n.season, length(covariates.det)))
    a.cov.occ[,,1] <- 1   # set the intercepts to 1
    a.cov.det[,,1] <- 1
    
    for(s in 1:n.season){    # set the occupancy covariate values from each season
      a.cov.occ[,s,-1] <- data.site.long %>% # for each covariate past the first (the intercept)...
        filter(scale == scale.model) %>%     # select the scale specified
        filter(season == unique(data.site.long$season)[s]) %>%  # select the season
        select(covariates.occ[-1]) %>% as.matrix()          # select the covariates, then add to the array of covariate values
    }
    
    a.cov.det[,,2] <- as.matrix(K.std)   # set the detection covariates
    
    # str(
      bdata <- list(
      y = ysum, 
      K = K.survey, 
      nsite = n.site, 
      nsurvey = n.survey,
      nseason = n.season,
      nspec = n.spp,
      psi.cov = a.cov.occ,
      ncov.occ = dim(a.cov.occ)[[3]],
      rho.cov = a.cov.det,
      ncov.det = dim(a.cov.det)[[3]],
      traits = traits.response,
      n.traits = dim(traits.response)[[2]],
      class_vec = class_vec,
      n.classes = length(unique(class_vec)),
      order_vec = order_vec,
      n.orders = length(unique(order_vec)))
    # )
    
    # Initial values 
    zst <- apply(y, c(1,3,4), function(x) max(x, na.rm = T))  # Observed occurrence as inits for z
    inits <- function(){ list(
      z = zst
    )}
    
    # Call JAGS, check convergence, summarize posteriors
    # n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (start.time <- Sys.time())
    out <- jags(bdata, inits, params, 
                "./code/jags/model_DCM.R", 
                n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (end.time <- Sys.time())
    elapsed.time <- difftime(end.time, start.time, units='mins')
    cat(paste(paste('Posterior computed in ', elapsed.time, sep=''), ' minutes', sep=''))
    
    #out$summary ; jags.View(out) ; print(out, 3) # not shown
    #par(mfrow = c(3,3)); traceplot(out, c("rich"))
    
    # Add the calculated summary statistic to the appropriate row
    (CPO$sumCPO[4] <- calc_cpo(out))
    
    for(i in 1:nc){
      lik <- as.matrix(out$samples[i])
      not.lik <- lik[,-c(grep("lik", colnames(lik)))]  # drop the samples of the 'lik' parameter
      out$samples[[i]] <- not.lik
    }
    
    msum <- out$summary
    head(msum)
    
    beep(sound = "coin")
    
    # Export model outputs 
    #saveRDS(out, "./data/model_outputs/DCM_univariate/univariateDCM_het1km.rds")
    # write.csv(msum, "./data/model_summaries/DCM_univariate/univariateDCM_impquad100m_summary.csv")
    
    
    ##### 1 km #####
    # repeat from above
    # select the covariates that will be included in the model here
    scale.model <- "1km"
    covariates.occ <- c(
      "intercept",
      "imp_std",
      "impquad"#,
      # "sdhi_std",
      # "ndvi_std"
    )
    covariates.det <- c(
      "intercept",
      "effort"
    )
    
    a.cov.occ <- array(NA, dim = c(n.site, n.season, length(covariates.occ))) # organize the covariates into arrays
    a.cov.det <- array(NA, dim = c(n.site, n.season, length(covariates.det)))
    a.cov.occ[,,1] <- 1   # set the intercepts to 1
    a.cov.det[,,1] <- 1
    
    for(s in 1:n.season){    # set the occupancy covariate values from each season
      a.cov.occ[,s,-1] <- data.site.long %>% # for each covariate past the first (the intercept)...
        filter(scale == scale.model) %>%     # select the scale specified
        filter(season == unique(data.site.long$season)[s]) %>%  # select the season
        select(covariates.occ[-1]) %>% as.matrix()          # select the covariates, then add to the array of covariate values
    }
    
    a.cov.det[,,2] <- as.matrix(K.std)   # set the detection covariates
    
    # str(
      bdata <- list(
      y = ysum, 
      K = K.survey, 
      nsite = n.site, 
      nsurvey = n.survey,
      nseason = n.season,
      nspec = n.spp,
      psi.cov = a.cov.occ,
      ncov.occ = dim(a.cov.occ)[[3]],
      rho.cov = a.cov.det,
      ncov.det = dim(a.cov.det)[[3]],
      traits = traits.response,
      n.traits = dim(traits.response)[[2]],
      class_vec = class_vec,
      n.classes = length(unique(class_vec)),
      order_vec = order_vec,
      n.orders = length(unique(order_vec)))
    # )
    
    
    # Initial values 
    zst <- apply(y, c(1,3,4), function(x) max(x, na.rm = T))  # Observed occurrence as inits for z
    inits <- function(){ list(
      z = zst
    )}
    
    # Call JAGS, check convergence, summarize posteriors
    # n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (start.time <- Sys.time())
    out <- jags(bdata, inits, params, 
                "./code/jags/model_DCM.R", 
                n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (end.time <- Sys.time())
    elapsed.time <- difftime(end.time, start.time, units='mins')
    cat(paste(paste('Posterior computed in ', elapsed.time, sep=''), ' minutes', sep=''))
    
    #out$summary ; jags.View(out) ; print(out, 3) # not shown
    #par(mfrow = c(3,3)); traceplot(out, c("rich"))
    
    # Add the calculated summary statistic to the appropriate row
    (CPO$sumCPO[5] <- calc_cpo(out))
    
    for(i in 1:nc){
      lik <- as.matrix(out$samples[i])
      not.lik <- lik[,-c(grep("lik", colnames(lik)))]  # drop the samples of the 'lik' parameter
      out$samples[[i]] <- not.lik
    }
    
    msum <- out$summary
    head(msum)
    
    beep(sound = "coin")
    
    # Export model outputs 
    #saveRDS(out, "./data/model_outputs/DCM_univariate/univariateDCM_het1km.rds")
    # write.csv(msum, "./data/model_summaries/DCM_univariate/univariateDCM_impquad1km_summary.csv")
    
    ##### 5 km #####
    # repeat from above
    # select the covariates that will be included in the model here
    scale.model <- "5km"
    covariates.occ <- c(
      "intercept",
      "imp_std",
      "impquad"#,
      # "sdhi_std",
      # "ndvi_std"
    )
    covariates.det <- c(
      "intercept",
      "effort"
    )
    
    a.cov.occ <- array(NA, dim = c(n.site, n.season, length(covariates.occ))) # organize the covariates into arrays
    a.cov.det <- array(NA, dim = c(n.site, n.season, length(covariates.det)))
    a.cov.occ[,,1] <- 1   # set the intercepts to 1
    a.cov.det[,,1] <- 1
    
    for(s in 1:n.season){    # set the occupancy covariate values from each season
      a.cov.occ[,s,-1] <- data.site.long %>% # for each covariate past the first (the intercept)...
        filter(scale == scale.model) %>%     # select the scale specified
        filter(season == unique(data.site.long$season)[s]) %>%  # select the season
        select(covariates.occ[-1]) %>% as.matrix()          # select the covariates, then add to the array of covariate values
    }
    
    a.cov.det[,,2] <- as.matrix(K.std)   # set the detection covariates
    
    # str(
      bdata <- list(
      y = ysum, 
      K = K.survey, 
      nsite = n.site, 
      nsurvey = n.survey,
      nseason = n.season,
      nspec = n.spp,
      psi.cov = a.cov.occ,
      ncov.occ = dim(a.cov.occ)[[3]],
      rho.cov = a.cov.det,
      ncov.det = dim(a.cov.det)[[3]],
      traits = traits.response,
      n.traits = dim(traits.response)[[2]],
      class_vec = class_vec,
      n.classes = length(unique(class_vec)),
      order_vec = order_vec,
      n.orders = length(unique(order_vec)))
    # )
    
    
    # Initial values 
    zst <- apply(y, c(1,3,4), function(x) max(x, na.rm = T))  # Observed occurrence as inits for z
    inits <- function(){ list(
      z = zst
    )}
    
    # Call JAGS, check convergence, summarize posteriors
    # n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (start.time <- Sys.time())
    out <- jags(bdata, inits, params, 
                "./code/jags/model_DCM.R", 
                n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (end.time <- Sys.time())
    elapsed.time <- difftime(end.time, start.time, units='mins')
    cat(paste(paste('Posterior computed in ', elapsed.time, sep=''), ' minutes', sep=''))
    
    #out$summary ; jags.View(out) ; print(out, 3) # not shown
    #par(mfrow = c(3,3)); traceplot(out, c("rich"))
    
    # Add the calculated summary statistic to the appropriate row
    (CPO$sumCPO[6] <- calc_cpo(out))
    
    for(i in 1:nc){
      lik <- as.matrix(out$samples[i])
      not.lik <- lik[,-c(grep("lik", colnames(lik)))]  # drop the samples of the 'lik' parameter
      out$samples[[i]] <- not.lik
    }
    
    msum <- out$summary
    head(msum)
    
    beep(sound = "coin")
    
    # Export model outputs 
    #saveRDS(out, "./data/model_outputs/DCM_univariate/univariateDCM_het1km.rds")
    # write.csv(msum, "./data/model_summaries/DCM_univariate/univariateDCM_impquad5km_summary.csv")
    
    
    
##### Heterogeneity at Different Scales #####
    ##### 100 m #####
    scale.model <- "100m"
    covariates.occ <- c(
      "intercept",
      # "imp_std",
      # "impquad"#,
      "sdhi_std"#,
      # "ndvi_std"
    )
    covariates.det <- c(
      "intercept",
      "effort"
    )
    
    a.cov.occ <- array(NA, dim = c(n.site, n.season, length(covariates.occ))) # organize the covariates into arrays
    a.cov.det <- array(NA, dim = c(n.site, n.season, length(covariates.det)))
    a.cov.occ[,,1] <- 1   # set the intercepts to 1
    a.cov.det[,,1] <- 1
    
    for(s in 1:n.season){    # set the occupancy covariate values from each season
      a.cov.occ[,s,-1] <- data.site.long %>% # for each covariate past the first (the intercept)...
        filter(scale == scale.model) %>%     # select the scale specified
        filter(season == unique(data.site.long$season)[s]) %>%  # select the season
        select(covariates.occ[-1]) %>% as.matrix()          # select the covariates, then add to the array of covariate values
    }
    
    a.cov.det[,,2] <- as.matrix(K.std)   # set the detection covariates
    
    # str(
      bdata <- list(
      y = ysum, 
      K = K.survey, 
      nsite = n.site, 
      nsurvey = n.survey,
      nseason = n.season,
      nspec = n.spp,
      psi.cov = a.cov.occ,
      ncov.occ = dim(a.cov.occ)[[3]],
      rho.cov = a.cov.det,
      ncov.det = dim(a.cov.det)[[3]],
      traits = traits.response,
      n.traits = dim(traits.response)[[2]],
      class_vec = class_vec,
      n.classes = length(unique(class_vec)),
      order_vec = order_vec,
      n.orders = length(unique(order_vec)))
    # )
    
    # Initial values 
    zst <- apply(y, c(1,3,4), function(x) max(x, na.rm = T))  # Observed occurrence as inits for z
    inits <- function(){ list(
      z = zst
    )}
    
    # Call JAGS, check convergence, summarize posteriors
    # n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (start.time <- Sys.time())
    out <- jags(bdata, inits, params, 
                "./code/jags/model_DCM.R", 
                n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (end.time <- Sys.time())
    elapsed.time <- difftime(end.time, start.time, units='mins')
    cat(paste(paste('Posterior computed in ', elapsed.time, sep=''), ' minutes', sep=''))
    
    #out$summary ; jags.View(out) ; print(out, 3) # not shown
    #par(mfrow = c(3,3)); traceplot(out, c("rich"))
    (CPO$sumCPO[7] <- calc_cpo(out))
    
    
    for(i in 1:nc){
      lik <- as.matrix(out$samples[i])
      not.lik <- lik[,-c(grep("lik", colnames(lik)))]  # drop the samples of the 'lik' parameter
      out$samples[[i]] <- not.lik
    }
    
    msum <- out$summary
    head(msum)
    
    beep(sound = "coin")
    
    # Export model outputs 
    #saveRDS(out, "./data/model_outputs/DCM_univariate/univariateDCM_het100m.rds")
    # write.csv(msum, "./data/model_summaries/DCM_univariate/univariateDCM_het100m_summary.csv")
    
    
    ##### 1 km #####
    scale.model <- "1km"
    covariates.occ <- c(
      "intercept",
      # "imp_std",
      # "impquad"#,
      "sdhi_std"#,
      # "ndvi_std"
    )
    covariates.det <- c(
      "intercept",
      "effort"
    )
    
    a.cov.occ <- array(NA, dim = c(n.site, n.season, length(covariates.occ))) # organize the covariates into arrays
    a.cov.det <- array(NA, dim = c(n.site, n.season, length(covariates.det)))
    a.cov.occ[,,1] <- 1   # set the intercepts to 1
    a.cov.det[,,1] <- 1
    
    for(s in 1:n.season){    # set the occupancy covariate values from each season
      a.cov.occ[,s,-1] <- data.site.long %>% # for each covariate past the first (the intercept)...
        filter(scale == scale.model) %>%     # select the scale specified
        filter(season == unique(data.site.long$season)[s]) %>%  # select the season
        select(covariates.occ[-1]) %>% as.matrix()          # select the covariates, then add to the array of covariate values
    }
    
    a.cov.det[,,2] <- as.matrix(K.std)   # set the detection covariates
    
    # str(
    bdata <- list(
      y = ysum, 
      K = K.survey, 
      nsite = n.site, 
      nsurvey = n.survey,
      nseason = n.season,
      nspec = n.spp,
      psi.cov = a.cov.occ,
      ncov.occ = dim(a.cov.occ)[[3]],
      rho.cov = a.cov.det,
      ncov.det = dim(a.cov.det)[[3]],
      traits = traits.response,
      n.traits = dim(traits.response)[[2]],
      class_vec = class_vec,
      n.classes = length(unique(class_vec)),
      order_vec = order_vec,
      n.orders = length(unique(order_vec)))
    # )
    
    
    # Initial values 
    zst <- apply(y, c(1,3,4), function(x) max(x, na.rm = T))  # Observed occurrence as inits for z
    inits <- function(){ list(
      z = zst
    )}
    
    # Call JAGS, check convergence, summarize posteriors
    # n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (start.time <- Sys.time())
    out <- jags(bdata, inits, params, 
                "./code/jags/model_DCM.R", 
                n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (end.time <- Sys.time())
    elapsed.time <- difftime(end.time, start.time, units='mins')
    cat(paste(paste('Posterior computed in ', elapsed.time, sep=''), ' minutes', sep=''))
    
    #out$summary ; jags.View(out) ; print(out, 3) # not shown
    #par(mfrow = c(3,3)); traceplot(out, c("rich"))
    
    # Add the calculated summary statistic to the appropriate row
    (CPO$sumCPO[8] <- calc_cpo(out))
    
    
    for(i in 1:nc){
      lik <- as.matrix(out$samples[i])
      not.lik <- lik[,-c(grep("lik", colnames(lik)))]  # drop the samples of the 'lik' parameter
      out$samples[[i]] <- not.lik
    }
    
    msum <- out$summary
    head(msum)
    
    beep(sound = "coin")
    
    # Export model outputs 
    #saveRDS(out, "./data/model_outputs/DCM_univariate/univariateDCM_het1km.rds")
    # write.csv(msum, "./data/model_summaries/DCM_univariate/univariateDCM_het1km_summary.csv")
    
    ##### 5 km #####
    scale.model <- "5km"
    covariates.occ <- c(
      "intercept",
      # "imp_std",
      # "impquad"#,
      "sdhi_std"#,
      # "ndvi_std"
    )
    covariates.det <- c(
      "intercept",
      "effort"
    )
    
    a.cov.occ <- array(NA, dim = c(n.site, n.season, length(covariates.occ))) # organize the covariates into arrays
    a.cov.det <- array(NA, dim = c(n.site, n.season, length(covariates.det)))
    a.cov.occ[,,1] <- 1   # set the intercepts to 1
    a.cov.det[,,1] <- 1
    
    for(s in 1:n.season){    # set the occupancy covariate values from each season
      a.cov.occ[,s,-1] <- data.site.long %>% # for each covariate past the first (the intercept)...
        filter(scale == scale.model) %>%     # select the scale specified
        filter(season == unique(data.site.long$season)[s]) %>%  # select the season
        select(covariates.occ[-1]) %>% as.matrix()          # select the covariates, then add to the array of covariate values
    }
    
    a.cov.det[,,2] <- as.matrix(K.std)   # set the detection covariates
    
    # str(
    bdata <- list(
      y = ysum, 
      K = K.survey, 
      nsite = n.site, 
      nsurvey = n.survey,
      nseason = n.season,
      nspec = n.spp,
      psi.cov = a.cov.occ,
      ncov.occ = dim(a.cov.occ)[[3]],
      rho.cov = a.cov.det,
      ncov.det = dim(a.cov.det)[[3]],
      traits = traits.response,
      n.traits = dim(traits.response)[[2]],
      class_vec = class_vec,
      n.classes = length(unique(class_vec)),
      order_vec = order_vec,
      n.orders = length(unique(order_vec)))
    # )
    
    
    # Initial values 
    zst <- apply(y, c(1,3,4), function(x) max(x, na.rm = T))  # Observed occurrence as inits for z
    inits <- function(){ list(
      z = zst
    )}
    
    # Call JAGS, check convergence, summarize posteriors
    # n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (start.time <- Sys.time())
    out <- jags(bdata, inits, params, 
                "./code/jags/model_DCM.R", 
                n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (end.time <- Sys.time())
    elapsed.time <- difftime(end.time, start.time, units='mins')
    cat(paste(paste('Posterior computed in ', elapsed.time, sep=''), ' minutes', sep=''))
    
    #out$summary ; jags.View(out) ; print(out, 3) # not shown
    #par(mfrow = c(3,3)); traceplot(out, c("rich"))
    
    # Add the calculated summary statistic to the appropriate row
    (CPO$sumCPO[9] <- calc_cpo(out))
    
    for(i in 1:nc){
      lik <- as.matrix(out$samples[i])
      not.lik <- lik[,-c(grep("lik", colnames(lik)))]  # drop the samples of the 'lik' parameter
      out$samples[[i]] <- not.lik
    }
    
    msum <- out$summary
    head(msum)
    
    beep(sound = "coin")
    
    # Export model outputs 
    #saveRDS(out, "./data/model_outputs/DCM_univariate/univariateDCM_het5km.rds")
    # write.csv(msum, "./data/model_summaries/DCM_univariate/univariateDCM_het5km_summary.csv")
  
  
    
##### Vegetation Greenness at Different Scales #####
    ##### 100 m #####
    scale.model <- "100m"
    covariates.occ <- c(
      "intercept",
      # "imp_std",
      # "impquad"#,
      # "sdhi_std"#,
      "ndvi_std"
    )
    covariates.det <- c(
      "intercept",
      "effort"
    )
    
    a.cov.occ <- array(NA, dim = c(n.site, n.season, length(covariates.occ))) # organize the covariates into arrays
    a.cov.det <- array(NA, dim = c(n.site, n.season, length(covariates.det)))
    a.cov.occ[,,1] <- 1   # set the intercepts to 1
    a.cov.det[,,1] <- 1
    
    for(s in 1:n.season){    # set the occupancy covariate values from each season
      a.cov.occ[,s,-1] <- data.site.long %>% # for each covariate past the first (the intercept)...
        filter(scale == scale.model) %>%     # select the scale specified
        filter(season == unique(data.site.long$season)[s]) %>%  # select the season
        select(covariates.occ[-1]) %>% as.matrix()          # select the covariates, then add to the array of covariate values
    }
    
    a.cov.det[,,2] <- as.matrix(K.std)   # set the detection covariates
    
    # str(
    bdata <- list(
      y = ysum, 
      K = K.survey, 
      nsite = n.site, 
      nsurvey = n.survey,
      nseason = n.season,
      nspec = n.spp,
      psi.cov = a.cov.occ,
      ncov.occ = dim(a.cov.occ)[[3]],
      rho.cov = a.cov.det,
      ncov.det = dim(a.cov.det)[[3]],
      traits = traits.response,
      n.traits = dim(traits.response)[[2]],
      class_vec = class_vec,
      n.classes = length(unique(class_vec)),
      order_vec = order_vec,
      n.orders = length(unique(order_vec)))
    # )
    
    
    # Initial values 
    zst <- apply(y, c(1,3,4), function(x) max(x, na.rm = T))  # Observed occurrence as inits for z
    inits <- function(){ list(
      z = zst
    )}
    
    # Call JAGS, check convergence, summarize posteriors
    # n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (start.time <- Sys.time())
    out <- jags(bdata, inits, params, 
                "./code/jags/model_DCM.R", 
                n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (end.time <- Sys.time())
    elapsed.time <- difftime(end.time, start.time, units='mins')
    cat(paste(paste('Posterior computed in ', elapsed.time, sep=''), ' minutes', sep=''))
    
    #out$summary ; jags.View(out) ; print(out, 3) # not shown
    #par(mfrow = c(3,3)); traceplot(out, c("rich"))
    (CPO$sumCPO[10] <- calc_cpo(out))
    
    for(i in 1:nc){
      lik <- as.matrix(out$samples[i])
      not.lik <- lik[,-c(grep("lik", colnames(lik)))]  # drop the samples of the 'lik' parameter
      out$samples[[i]] <- not.lik
    }
    
    msum <- out$summary
    head(msum)
    
    beep(sound = "coin")
    
    # Export model outputs 
    #saveRDS(out, "./data/model_outputs/DCM_univariate/univariateDCM_veg100m.rds")
    # write.csv(msum, "./data/model_summaries/DCM_univariate/univariateDCM_veg100m_summary.csv")
    
    
    ##### 1 km #####
    scale.model <- "1km"
    covariates.occ <- c(
      "intercept",
      # "imp_std",
      # "impquad"#,
      # "sdhi_std"#,
      "ndvi_std"
    )
    covariates.det <- c(
      "intercept",
      "effort"
    )
    
    a.cov.occ <- array(NA, dim = c(n.site, n.season, length(covariates.occ))) # organize the covariates into arrays
    a.cov.det <- array(NA, dim = c(n.site, n.season, length(covariates.det)))
    a.cov.occ[,,1] <- 1   # set the intercepts to 1
    a.cov.det[,,1] <- 1
    
    for(s in 1:n.season){    # set the occupancy covariate values from each season
      a.cov.occ[,s,-1] <- data.site.long %>% # for each covariate past the first (the intercept)...
        filter(scale == scale.model) %>%     # select the scale specified
        filter(season == unique(data.site.long$season)[s]) %>%  # select the season
        select(covariates.occ[-1]) %>% as.matrix()          # select the covariates, then add to the array of covariate values
    }
    
    a.cov.det[,,2] <- as.matrix(K.std)   # set the detection covariates
    
    # str(
    bdata <- list(
      y = ysum, 
      K = K.survey, 
      nsite = n.site, 
      nsurvey = n.survey,
      nseason = n.season,
      nspec = n.spp,
      psi.cov = a.cov.occ,
      ncov.occ = dim(a.cov.occ)[[3]],
      rho.cov = a.cov.det,
      ncov.det = dim(a.cov.det)[[3]],
      traits = traits.response,
      n.traits = dim(traits.response)[[2]],
      class_vec = class_vec,
      n.classes = length(unique(class_vec)),
      order_vec = order_vec,
      n.orders = length(unique(order_vec)))
    # )
    
    # Initial values 
    zst <- apply(y, c(1,3,4), function(x) max(x, na.rm = T))  # Observed occurrence as inits for z
    inits <- function(){ list(
      z = zst
    )}
    
    # Call JAGS, check convergence, summarize posteriors
    # n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (start.time <- Sys.time())
    out <- jags(bdata, inits, params, 
                "./code/jags/model_DCM.R", 
                n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (end.time <- Sys.time())
    elapsed.time <- difftime(end.time, start.time, units='mins')
    cat(paste(paste('Posterior computed in ', elapsed.time, sep=''), ' minutes', sep=''))
    
    #out$summary ; jags.View(out) ; print(out, 3) # not shown
    #par(mfrow = c(3,3)); traceplot(out, c("rich"))
    
    # Add the calculated summary statistic to the appropriate row
    (CPO$sumCPO[11] <- calc_cpo(out))
    
    for(i in 1:nc){
      lik <- as.matrix(out$samples[i])
      not.lik <- lik[,-c(grep("lik", colnames(lik)))]  # drop the samples of the 'lik' parameter
      out$samples[[i]] <- not.lik
    }
    
    msum <- out$summary
    head(msum)
    
    beep(sound = "coin")
    
    # Export model outputs 
    #saveRDS(out, "./data/model_outputs/DCM_univariate/univariateDCM_veg1km.rds")
    # write.csv(msum, "./data/model_summaries/DCM_univariate/univariateDCM_veg1km_summary.csv")
    
    ##### 5 km #####
    scale.model <- "5km"
    covariates.occ <- c(
      "intercept",
      # "imp_std",
      # "impquad"#,
      # "sdhi_std"#,
      "ndvi_std"
    )
    covariates.det <- c(
      "intercept",
      "effort"
    )
    
    a.cov.occ <- array(NA, dim = c(n.site, n.season, length(covariates.occ))) # organize the covariates into arrays
    a.cov.det <- array(NA, dim = c(n.site, n.season, length(covariates.det)))
    a.cov.occ[,,1] <- 1   # set the intercepts to 1
    a.cov.det[,,1] <- 1
    
    for(s in 1:n.season){    # set the occupancy covariate values from each season
      a.cov.occ[,s,-1] <- data.site.long %>% # for each covariate past the first (the intercept)...
        filter(scale == scale.model) %>%     # select the scale specified
        filter(season == unique(data.site.long$season)[s]) %>%  # select the season
        select(covariates.occ[-1]) %>% as.matrix()          # select the covariates, then add to the array of covariate values
    }
    
    a.cov.det[,,2] <- as.matrix(K.std)   # set the detection covariates
    
    # str(
    bdata <- list(
      y = ysum, 
      K = K.survey, 
      nsite = n.site, 
      nsurvey = n.survey,
      nseason = n.season,
      nspec = n.spp,
      psi.cov = a.cov.occ,
      ncov.occ = dim(a.cov.occ)[[3]],
      rho.cov = a.cov.det,
      ncov.det = dim(a.cov.det)[[3]],
      traits = traits.response,
      n.traits = dim(traits.response)[[2]],
      class_vec = class_vec,
      n.classes = length(unique(class_vec)),
      order_vec = order_vec,
      n.orders = length(unique(order_vec)))
    # )
    
    
    # Initial values 
    zst <- apply(y, c(1,3,4), function(x) max(x, na.rm = T))  # Observed occurrence as inits for z
    inits <- function(){ list(
      z = zst
    )}
    
    # Call JAGS, check convergence, summarize posteriors
    # n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (start.time <- Sys.time())
    out <- jags(bdata, inits, params, 
                "./code/jags/model_DCM.R", 
                n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (end.time <- Sys.time())
    elapsed.time <- difftime(end.time, start.time, units='mins')
    cat(paste(paste('Posterior computed in ', elapsed.time, sep=''), ' minutes', sep=''))
    
    #out$summary ; jags.View(out) ; print(out, 3) # not shown
    #par(mfrow = c(3,3)); traceplot(out, c("rich"))
    
    # Add the calculated summary statistic to the appropriate row
    (CPO$sumCPO[12] <- calc_cpo(out))
    
    for(i in 1:nc){
      lik <- as.matrix(out$samples[i])
      not.lik <- lik[,-c(grep("lik", colnames(lik)))]  # drop the samples of the 'lik' parameter
      out$samples[[i]] <- not.lik
    }
    
    msum <- out$summary
    head(msum)
    
    beep(sound = "coin")
    
    # Export model outputs 
    #saveRDS(out, "./data/model_outputs/DCM_univariate/univariateDCM_veg5km.rds")
    # write.csv(msum, "./data/model_summaries/DCM_univariate/univariateDCM_veg5km_summary.csv")
    
    
    
### Export Model CPO Results  
    # write.csv(CPO, "./data/model_summaries/DCM_univariate/univariateDCM_modelselectionCPO.csv")
    
    
#### Examine Univariate Results before Running Global Models ####
    CPO
    # top models
    # urbanization: 1 km
    # heterogeneity: 1 km 
    # veg: 1 km
    # quad urbanization: 1 km
    
    # Based solely on the univariate models, those would be the optimal scales for each covariate
    # However, we want to limit collinearity and some of these are highly correlated
    
    # across all seasons and scales, these variables are not highly correlated
    ggcorrplot(cor((data.site.long %>% select(c(
      imp, impquad, sdhi, ndvi
    ))), 
    use = "complete.obs", method = "pearson"), 
    hc.order = FALSE, 
    type = "lower",
    lab = TRUE,
    lab_size = 2.5,
    outline.color = "white") 
    
    # but the "optimized" scales create some potential collinearity issues
    ggcorrplot(cor((data.site %>% select(c(
      # imp100m, 
      imp1km, 
      # imp5km,
      # impquad100m, 
      impquad1km, 
      # impquad5km,
      sdhi100m, sdhi1km, sdhi5km,
      #ndvi_wd100m, ndvi_ww100m, ndvi_cw100m,
      ndvi_wd1km, ndvi_ww1km, ndvi_cw1km,
      #ndvi_wd5km, ndvi_ww5km, ndvi_cw5km
    ))), 
    use = "complete.obs", method = "pearson"), 
    hc.order = FALSE, 
    type = "lower",
    lab = TRUE,
    lab_size = 2.5,
    outline.color = "white")
    
    # At the 100m and 1km scales, patch diversity is more highly correlated  (-0.59 and -0.70) with quadratic urbanization
    # while 5 km patch diversity is less correlated (-0.42)
    # let's inspect this visually
    plot(data.site$impquad1km, data.site$sdhi100m)  # lots of 0s, due to the scale being so small
    plot(data.site$impquad1km, data.site$sdhi1km)   # a relatively strong correlation
    plot(data.site$impquad1km, data.site$sdhi5km)   # less of a correlation
    
    plot(data.site$imp1km, data.site$sdhi100m)  # lots of 0s, due to the scale being so small
    plot(data.site$imp1km, data.site$sdhi1km)   # a clearly peaked relationship
    plot(data.site$imp1km, data.site$sdhi5km)   # less of a peak and more of a negative relationship
    
    
    # NDVI is correlated with urbanization, but only during the cool-wet season
    
    
    