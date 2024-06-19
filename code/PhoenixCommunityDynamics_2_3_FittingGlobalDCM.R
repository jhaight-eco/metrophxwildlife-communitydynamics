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

# standardized versions of seasonally-varying covariates
# these were the only ones that were separate objects, because they were standardized by season
str(K.std)
str(ndvi.1km.std)
str(ndvihet.1km.std)

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

# we will store these CPO summary statistic values in a small dataframe
CPO <- data.frame(
  "model" = c("global (no quad urbanization)", "global (with quad urbanization"),
  "sumCPO" = rep(NA, times = 2)
)
CPO


#### Setup for MCMC ####
# Across all models, we will use these same settings and monitor the same parameters
# MCMC settings
na <- 100	   # pre-burnin
nb <- 370	 # burn-in
ni <- 400  # iterations (including burn-in with jagsUI)
nt <- 3
nc <- 3
str(ysum)   # check that this is integer data (will not run otherwise)



# Parameters monitored (could also add "z")
params <- c(
  # community-level intercepts
    "mu.lpsi1",
    "mu.lphi",
    "mu.lgamma",
    "mu.lp",
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
    "lp.year",
    "lp",
    "lpsi1",
    "lphi",
    "lgamma",
  # Species-level occurrences
    #"n.occ",
    #"z",
  # Likelihood
    "lik", 
  # Species-level Real Parameters
    #"p",
    #"gamma",
    #"psi1", 
    #"phi",
    "psi"
)

# Bundle data for BUGS
str(bdata <- list(
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
  #det1 = K$effort_std,
  det1 = as.matrix(K.std)
)
)

# Initial values 
zst <- apply(y, c(1,3,4), function(x) max(x, na.rm = T))  # Observed occurrence as inits for z
inits <- function(){ list(
  z = zst
)}


### Run Models ###
#### Global Model, minus the Quadratic Urbanization term ####
# Call JAGS, check convergence, summarize posteriors
# Approximate run time: 20-30 min on 
(start.time <- Sys.time())
out <- jags(bdata, inits, params, 
            "./code/jags/model_DCM_global_noquadurb.R", 
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
saveRDS(out, "./data/model_outputs/DCM_global/globalDCMnoquadurb_sample3k.rds")
write.csv(msum, "./data/model_outputs/DCM_global/globalDCMnoquadurb_sample3k_summary.csv")



#### Global Model (including the quadratic urbanization term) ####
# Call JAGS, check convergence, summarize posteriors
# Approximate run time: 20-30 min on 
(start.time <- Sys.time())
out <- jags(bdata, inits, params, 
            "./code/jags/model_DCM_global.R", 
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
saveRDS(out, "./data/model_outputs/DCM_global/globalDCM_sample3k.rds")
write.csv(msum, "./data/model_outputs/DCM_global/globalDCM_sample3k_summary.csv")





write.csv(CPO, "./data/model_outputs/DCM_global/global_modelselectionCPO.csv", row.names = FALSE)
