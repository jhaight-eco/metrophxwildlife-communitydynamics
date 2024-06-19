

#### Setup ####
  rm(list=ls()) # clear the environment
  gc()
  set.seed(54321)
  
  # Set working directory
  setwd("~/GitHub/metrophxwildlife-communitydynamics")
  
  # Load necessary packages
  library(jagsUI)   # for fitting models in JAGS
  library(ggplot2)
  library(abind)
  library(dplyr)
  library(beepr)



#### Import Model Input Dataset ####
  load("./data/model_inputs/ModelInputData_DCM.RData")
  str(data.site)
  

#### Import Model Output File ####
  # the full original version of the model is too large to upload to the GitHub repository
  # this file can be recreated using the script 'PhoenixCommunityDynamics_2_3_FittingGlobalDCM.R'
  out <- readRDS("./data/model_outputs/DCM_global/globalDCM_sample30k.rds")    


#### Summarize Species Richness ####
  str(
    tmp <- out$sims.list
    )
  rm(out)
  gc()
  
  pm_rich <- apply(tmp$rich, c(2,3), function(x)   mean(x, na.rm=TRUE))    # posterior mean
  sd_rich <- apply(tmp$rich, c(2,3), function(x)   sd(x, na.rm=TRUE))    # posterior mean
  pmed_rich <- apply(tmp$rich, c(2,3), function(x)   median(x, na.rm=TRUE))    # posterior median
  psd_rich <- apply(tmp$rich, c(2,3), function(x)   sd(x, na.rm=TRUE))    # posterior standard deviation
  cri_rich <- apply(tmp$rich, c(2,3), function(x)   quantile(x, prob = c(0.025, 0.975, 0.05, 0.95))) # posterior quantiles
  
  # add the predicted richness columns to the site dataset
  data.site$rich_wd_mean <- pm_rich[,1]
  data.site$rich_wd_med <- pmed_rich[,1]
  data.site$rich_wd_sd <- psd_rich[,1]
  data.site$rich_wd_2_5 <- c(t(cri_rich[1,,1]))
  data.site$rich_wd_97_5 <- c(t(cri_rich[2,,1]))
  
  data.site$rich_ww_mean <- pm_rich[,2]
  data.site$rich_ww_med <- pmed_rich[,2]
  data.site$rich_ww_sd <- psd_rich[,2]
  data.site$rich_ww_2_5 <- c(t(cri_rich[1,,2]))
  data.site$rich_ww_97_5 <- c(t(cri_rich[2,,2]))
  
  data.site$rich_cw_mean <- pm_rich[,3]
  data.site$rich_cw_med <- pmed_rich[,3]
  data.site$rich_cw_sd <- psd_rich[,3]
  data.site$rich_cw_2_5 <- c(t(cri_rich[1,,3]))
  data.site$rich_cw_97_5 <- c(t(cri_rich[2,,3]))
  

#### Plot Model Results ####
  # observed vs. estimated species richness
  # all points should be above a 1-1 line
  plot(data.site$rich_obs_wd, data.site$rich_wd_mean, xlim = c(0,16), ylim = c(0,16))
  plot(data.site$rich_obs_ww, data.site$rich_ww_mean, xlim = c(0,16), ylim = c(0,16))
  plot(data.site$rich_obs_cw, data.site$rich_cw_mean, xlim = c(0,16), ylim = c(0,16))
  
  # a quick look at estimated richness vs. urbanization in the three seasons
  plot(data.site$imp1km, data.site$rich_wd_mean, ylim = c(0,16))
  plot(data.site$imp1km, data.site$rich_ww_mean, ylim = c(0,16))
  plot(data.site$imp1km, data.site$rich_cw_mean, ylim = c(0,16))
    
#### Set Up Species Richness Model ####
  # make the seasonal design matrices for the analysis, stacking all the seasons into one array
  # vectorizing these simplifies the JAGS code
  dm.alpha.all <- abind(
    cbind(
      1,
      data.site$imp1km_std,
      data.site$imp1km_std^2,
      data.site$sdhi5km_std#,
      #ndvi.1km.std[,1]
    ),
    cbind(
      1,
      data.site$imp1km_std,
      data.site$imp1km_std^2,
      data.site$sdhi5km_std#,
      #ndvi.1km.std[,2]
    ),
    cbind(
      1,
      data.site$imp1km_std,
      data.site$imp1km_std^2,
      data.site$sdhi5km_std#,
      #ndvi.1km.std[,3]
    ),
    along = 3
  )
  
  # and give them useful names
  dimnames(dm.alpha.all)[[2]] <- c(
    "intercept",
    "urbanization",
    "urbanization quad",
    "heterogeneity"#,
    #"vegetation greenness"
  )
  
  # create the data lists for the analysis
  bdata <- list(
    alpha_z = pm_rich,
    alpha_sd_known = sd_rich,
    dm_alpha = dm.alpha.all,
    nsite = n.site,
    nseason = n.season,
    npar_alpha = ncol(dm.alpha.all)
  )
  
#### Fit Species Richness Model #####
  # Initial values function for JAGS
  inits_rich <- function() { list(
    beta = matrix(rep(runif(bdata$npar_alpha, 0, 0.5), n.season), ncol = n.season),
    my_re = matrix(rep(rnorm(bdata$nsite), n.season), ncol = n.season),
    re_tau = rep(dgamma(1,1,1), n.season)
  )
  }
  
  # MCMC settings
  # for jagsUI, number of samples = ni - nb
  nc <- 3       # number of chains
  nt <- 3       # number to samples thin by
  na <- 10000 # number of adaptations
  nb <- 370000   # number of burn-ins
  ni <- 460000   # number of iterations
  # for illustrative purposes, the next three parameters are 1/10th of their value in the manuscript
  #na <- 1000 # number of adaptations
  #nb <- 37000   # number of burn-ins
  #ni <- 40000   # number of iterations
  
  # fit model in JAGS, using 'jagsUI'
  (start.time <- Sys.time())
  m.sr <- jags(
    model.file = "./code/jags/model_alphadiversity_multiseason.R",
    data = bdata,
    n.chains = nc,
    parameters.to.save = c("beta", "re_sd"),
    n.adapt = na,
    n.burnin = nb,
    n.iter = ni,
    n.thin = nt,
    modules = "glm",   # log-normal model
    inits = inits_rich,
    parallel = TRUE
  )
  (end.time <- Sys.time())
  elapsed.time <- difftime(end.time, start.time, units='mins')
  cat(paste(paste('Posterior computed in ', elapsed.time, sep=''), ' minutes', sep=''))
  beep(sound = "coin")
  
  # summarize the model
  m.sr.sum <- data.frame(
    "parameter" = c(rep(colnames(dm.alpha.all), n.season), rep("residual error", 3), "deviance"),
    "season" = c(rep(c("spring","summer","winter"), each = ncol(dm.alpha.all)), c("spring","summer","winter"), NA),
    m.sr$summary
    )
  
  
  # view the summary
  m.sr.sum
  
  saveRDS(m.sr, "./data/model_outputs/metaanalysislogglm_speciesrichness_30k.rds")
  write.csv(m.sr.sum,"./data/model_summaries/communitymetaanalysisresults_richness.csv", row.names = FALSE)
  
  
  
#### Export Summarized Richness Estimates ####
  # include all the relevant site variables
  data.site$ndvi_wd100m_std <- ndvi.100m.std$ndvi_wd_std
  data.site$ndvi_ww100m_std <- ndvi.100m.std$ndvi_ww_std
  data.site$ndvi_cw100m_std <- ndvi.100m.std$ndvi_cw_std
  
  data.site$ndvi_wd1km_std <- ndvi.1km.std$ndvi_wd_std
  data.site$ndvi_ww1km_std <- ndvi.1km.std$ndvi_ww_std
  data.site$ndvi_cw1km_std <- ndvi.1km.std$ndvi_cw_std
  
  data.site$ndvi_wd5km_std <- ndvi.5km.std$ndvi_wd_std
  data.site$ndvi_ww5km_std <- ndvi.5km.std$ndvi_ww_std
  data.site$ndvi_cw5km_std <- ndvi.5km.std$ndvi_cw_std
  
  data.export <- data.site %>%
    dplyr::select(-c(
      X, urb100m_15, urb1km_15, veg100m_15, veg1km_15
    ))
  
  
  write.csv(data.export, "./data/model_outputs/data_PhoenixCommDynamics.csv", row.names = FALSE)
  
  