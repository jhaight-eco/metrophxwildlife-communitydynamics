rm(list=ls())
gc()


library(jagsUI)
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
    str(K)        # sampling effort = number of days each site was surveyed
  
    # Species covariates
    str(data.spp)
    dimnames(ysum)[[3]]  # the data is sorted by 'name_short'
    data.spp <- data.spp %>% arrange(name_short) # sort the species info to match
    
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
    
    # we will store these CPO summary statistic values in a small dataframe
    CPO <- data.frame(
      "model" = c("global (no quad urbanization)", "global (with quad urbanization)"),
      "sumCPO" = rep(NA, times = 2)
    )
    CPO


#### Setup for MCMC ####
    # Across all models, we will use these same settings and monitor the same parameters
    # MCMC settings
    na <- 1000	   # pre-burnin
    nb <- 37000	 # burn-in
    ni <- 40000  # iterations (including burn-in with jagsUI)
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
      # "eff.class.psi",
      # "eff.class.phi",
      # "eff.class.gamma",
      "eff.order.psi",
      "eff.order.phi",
      "eff.order.gamma",
      # species-level intercepts and slopes
      "beta.sp.psi1",
      "beta.sp.phi", 
      "beta.sp.gamma",
      "rho.sp.tot",
      # "rho.sp",
      # Community Composition Metrics
      "rich",
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
      # "psi"
    )
    


### Run Models ###
    #### Global Model, minus the Quadratic Urbanization term ####
        scale.model <- c("none",
                         "1km", 
                         # "1km",
                         "5km",
                         "1km")
        covariates.occ <- c(
          "intercept",
          "imp_std",
          # "impquad",
          "sdhi_std",
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
        
        for(l in 2:length(covariates.occ)){
          for(s in 1:n.season){    # set the occupancy covariate values from each season
            a.cov.occ[,s,l] <- data.site.long %>% # for each covariate past the first (the intercept)...
              filter(scale == scale.model[l]) %>%     # select the scale specified
              filter(season == unique(data.site.long$season)[s]) %>%  # select the season
              select(covariates.occ[l]) %>% as.matrix()          # select the covariates, then add to the array of covariate values
          }
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
          z = zst,
          beta.comm = rnorm(bdata$ncov.occ, -0.1, 0.1),
          tau.beta.comm = rep(dgamma(1, 1), bdata$ncov.occ),
          rho.comm = rnorm(bdata$ncov.det, -0.1, 0.1),
          tau.rho.comm = rep(dgamma(1, 1), bdata$ncov.det)#,
          # eff.traits.psi = matrix(
          #   rnorm(bdata$ncov.occ, -0.1, 0.1),
          #   ncol = bdata$n.traits, nrow = bdata$ncov.occ
          # ),
          # eff.traits.phi = matrix(
          #   rnorm(bdata$ncov.occ, -0.1, 0.1),
          #   ncol = bdata$n.traits, nrow = bdata$ncov.occ
          # ),
          # eff.traits.gamma = matrix(
          #   rnorm(bdata$ncov.occ, -0.1, 0.1),
          #   ncol = bdata$n.traits, nrow = bdata$ncov.occ
          # )#,
          # beta.species = matrix(
          #   rnorm(bdata$n.species * bdata$ncov.occ, 0, 0.1),
          #   ncol = bdata$n.species,
          #   nrow = bdata$ncov.occ
          # ),
          # rho.species = matrix(
          #   rnorm(bdata$n.species * bdata$ncov.det, -2.05, 0.1),
          #   ncol = bdata$n.species,
          #   nrow = bdata$ncov.det
          # )
        )}
    
        
    # Call JAGS, check convergence, summarize posteriors
        
    # Approximate run time: 20-30 min with 3k samples  
    (start.time <- Sys.time())
    out <- jags(bdata, inits, params, 
                "./code/2_0_model_DCM.R", 
                n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (end.time <- Sys.time())
    elapsed.time <- difftime(end.time, start.time, units='mins')
    cat(paste(paste('Posterior computed in ', elapsed.time, sep=''), ' minutes', sep=''))
    
    #out$summary ; jags.View(out) ; print(out, 3) # not shown
    #par(mfrow = c(3,3)); traceplot(out, c("rich"))
    
    msum <- out$summary
    head(msum)
    
    (CPO$sumCPO[1] <- calc_cpo(out))
    #beep(sound = "coin")
    
    
    # drop 'lik' from the samples, to reduce model object size before exporting 
    # now that the summary stat has been calculated, those aren't necessary and take up a lot of memory
    for(i in 1:nc){
      lik <- as.matrix(out$samples[i])
      not.lik <- lik[,-c(grep("lik", colnames(lik)))]  # drop the samples of the 'lik' parameter
      out$samples[[i]] <- not.lik
    }
    
    
    # Export model outputs 
    # saveRDS(out, "./data/model_outputs/globalDCMnoquadurb_sample3k.rds")
    # write.csv(msum, "./data/model_summaries/globalDCMnoquadurb_sample3k_summary.csv")
    
    

#### Global Model (including the quadratic urbanization term) ####
    scale.model <- c("none",
                     "1km", 
                     "1km",
                     "5km",
                     "1km")
    covariates.occ <- c(
      "intercept",
      "imp_std",
      "impquad",
      "sdhi_std",
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
    
    for(l in 2:length(covariates.occ)){
      for(s in 1:n.season){    # set the occupancy covariate values from each season
        a.cov.occ[,s,l] <- data.site.long %>% # for each covariate past the first (the intercept)...
          filter(scale == scale.model[l]) %>%     # select the scale specified
          filter(season == unique(data.site.long$season)[s]) %>%  # select the season
          select(covariates.occ[l]) %>% as.matrix()          # select the covariates, then add to the array of covariate values
      }
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
      z = zst,
      beta.comm = rnorm(bdata$ncov.occ, -0.1, 0.1),
      tau.beta.comm = rep(dgamma(1, 1), bdata$ncov.occ),
      rho.comm = rnorm(bdata$ncov.det, -0.1, 0.1),
      tau.rho.comm = rep(dgamma(1, 1), bdata$ncov.det)#,
      # eff.traits.psi = matrix(
      #   rnorm(bdata$ncov.occ, -0.1, 0.1),
      #   ncol = bdata$n.traits, nrow = bdata$ncov.occ
      # ),
      # eff.traits.phi = matrix(
      #   rnorm(bdata$ncov.occ, -0.1, 0.1),
      #   ncol = bdata$n.traits, nrow = bdata$ncov.occ
      # ),
      # eff.traits.gamma = matrix(
      #   rnorm(bdata$ncov.occ, -0.1, 0.1),
      #   ncol = bdata$n.traits, nrow = bdata$ncov.occ
      # )#,
      # beta.species = matrix(
      #   rnorm(bdata$n.species * bdata$ncov.occ, 0, 0.1),
      #   ncol = bdata$n.species,
      #   nrow = bdata$ncov.occ
      # ),
      # rho.species = matrix(
      #   rnorm(bdata$n.species * bdata$ncov.det, -2.05, 0.1),
      #   ncol = bdata$n.species,
      #   nrow = bdata$ncov.det
      # )
    )}
    
    # Call JAGS, check convergence, summarize posteriors
    # Approximate run time: 20-30 min with 3k samples
    (start.time <- Sys.time())
    out <- jags(bdata, inits, params, 
                "./code/2_0_model_DCM.R", 
                n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T) 
    (end.time <- Sys.time())
    elapsed.time <- difftime(end.time, start.time, units='mins')
    cat(paste(paste('Posterior computed in ', elapsed.time, sep=''), ' minutes', sep=''))
    
    #out$summary ; jags.View(out) ; print(out, 3) # not shown
    #par(mfrow = c(3,3)); traceplot(out, c("rich"))
    
    msum <- out$summary
    head(msum)
    
    (CPO$sumCPO[2] <- calc_cpo(out))
    #beep(sound = "coin")
    
    
    # drop 'lik' from the samples, to reduce model object size before exporting 
    # now that the summary stat has been calculated, those aren't necessary and take up a lot of memory
    for(i in 1:nc){
      lik <- as.matrix(out$samples[i])
      not.lik <- lik[,-c(grep("lik", colnames(lik)))]  # drop the samples of the 'lik' parameter
      out$samples[[i]] <- not.lik
    }
    
    
    # Export model outputs 
    # saveRDS(out, "./data/model_outputs/globalDCM_sample3k.rds")
    # write.csv(msum, "./data/model_summaries/globalDCM_sample3k_summary.csv")
    
    # saveRDS(out, "./data/model_outputs/globalDCM_sample30k.rds")
    # write.csv(msum, "./data/model_summaries/globalDCM_sample30k_summary.csv")



#### Export the model selection CPO summary table ####
  write.csv(CPO, "./data/model_summaries/global_modelselectionCPO.csv", row.names = FALSE)
