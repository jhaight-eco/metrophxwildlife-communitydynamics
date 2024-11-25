# Input data
    # y[] =  3D array with values for the number of survey occasions each species was detected at each sites, in each season
    # K[] = # of occasions in which each site was surveyed per season, 
    # nsite = number of sites surveyed
    # nsurvey = maximum number of survey occasions per site
    # nseason = number of sampling seasons
    # nspec = number of species modeled
    # psi.cov = matrix of standardized occupancy covariates for the linear predictor, with the first column representing the intercept (1)
    # ncov.occ = # of occupancy covariates + 1 (for the intercept)
    # rho.cov = matrix of standardized detection covariates for the linear predictor, with the first column representing the intercept (1)
    # ncov.det = # of detection covariates + 1 (for the intercept)
    # traits = matrix of standardized species trait covariates
    # n.traits = number of species trait covariates
    # order_vec = a vector of integers corresponding to each of the taxonomic orders represented in the community
    # n.orders = # of taxonomic orders represented in the community

# Modeled parameters
    # beta.comm.psi1 = community-level hyperparameters for initial occupancy, intercepts and slopes
    # beta.comm.phi = community-level hyperparameters for persistence, intercepts and slopes 
    # beta.comm.gamma = community-level hyperparameters for colonization, intercepts and slopes
    # beta.sp.psi1 = species-level parameters for initial occupancy, intercepts and slopes
    # beta.sp.phi = species-level parameters for persistence, intercepts and slopes 
    # beta.sp.gamma = species-level parameters for colonization, intercepts and slopes
    # rho.comm = community-level hyperparameters for detection probability, intercept and slope(s)
    # rho.sp.tot = species-level, season-averaged parameters for detection probability, intercept and slope(s)
    # rho.sp = species-level, season-specific parameters for detection probability, intercept and slope(s)
    # eff.traits.psi/eff.traits.phi/eff.traits.gamma = random effects of species traits on habitat use (included for hypothesis testing)
    # eff.order.psi/eff.order.phi/eff.order.gamma = random effect of taxonomic order on habitat use (added to partially control for phylogenetic effects)
    # lik = likelihood of each species observation outcome, used to calculate CPO

model {
# Specify priors
  # Priors for Occupancy (beta, random intercepts and slopes)
  # beta.comm.psi1[1] ~ dunif(0, 1)  
  for(l in 1:ncov.occ){
    beta.comm.psi1[l] ~ dnorm(0, 0.1)  # Community-mean intercept and slopes
    # beta.comm.psi1[l] ~ dt(0, 2.5, 1)
    tau.beta.comm.psi1[l] ~ dgamma(1, 1)    # Dispersion of community-mean parameters
    
    # Persistence: how likely a species is to remain at a site
    beta.comm.phi[l] ~ dnorm(0, 0.1)  # Community-mean intercept and slopes
    # beta.comm.phi[l] ~ dt(0, 2.5, 1)
    tau.beta.comm.phi[l] ~ dgamma(1, 1)    # Dispersion of community-mean parameters
    
    # Colonization: how likely a species is to use a site it wasn't already using
    beta.comm.gamma[l] ~ dnorm(0, 0.1)  # Community-mean intercept and slopes
    # beta.comm.gamma[l] ~ dt(0, 2.5, 1)
    tau.beta.comm.gamma[l] ~ dgamma(1, 1)    # Dispersion of community-mean parameters
    
    # Priors for Trait Effects on intercepts and slopes (for occupancy, persistence, and colonization)
      # Including these allows species-level intercepts and slopes (environmental responses) 
      # to vary from the community-average according to traits:
      for(m in 1:n.traits){
        eff.traits.psi[l,m] ~ dnorm(0, 0.1) # priors for species covariate effects
        eff.traits.phi[l,m] ~ dnorm(0, 0.1) # priors for species covariate effects
        eff.traits.gamma[l,m] ~ dnorm(0, 0.1) # priors for species covariate effects
      }
    
    # Priors for the taxonomic random effects
      # Structured similarly to the trait effects, including these allows for the model to account for some 
      # some random variation/error among species to be attributed to taxonomic/phylogenetic biases
      # in species responses to 
      tau.order.psi[l] ~ dgamma(0.1, 0.1)
      tau.order.phi[l] ~ dgamma(0.1, 0.1)
      tau.order.gamma[l] ~ dgamma(0.1, 0.1)
      for(k in 1:n.orders){
        eff.order.psi[l,k] ~ dnorm(0, tau.order.psi[l])
        eff.order.phi[l,k] ~ dnorm(0, tau.order.phi[l])
        eff.order.gamma[l,k] ~ dnorm(0, tau.order.gamma[l])
      }
    
    # structure variation of species-specific intercept and slope parameters from community-mean
    for(j in 1:nspec){
      # If species vary unpredictably from the community mean:
      # beta.species[l,j] ~ dnorm(beta.comm[l], tau.beta.comm[l])  # Species-specific
      # If species parameters are allowed to partially vary according to traits:
      beta.sp.psi1[l,j] ~ dnorm(beta.comm.psi1[l] 
                                + inprod(eff.traits.psi[l,1:n.traits], traits[j,1:n.traits]) 
                                # + eff.class.psi[l,class_vec[j]]
                                + eff.order.psi[l,order_vec[j]]
                                , tau.beta.comm.psi1[l])
      beta.sp.phi[l,j] ~ dnorm(beta.comm.phi[l] 
                               + inprod(eff.traits.phi[l,1:n.traits], traits[j,1:n.traits]) 
                               # + eff.class.phi[l,class_vec[j]]
                               + eff.order.phi[l,order_vec[j]]
                               , tau.beta.comm.phi[l])
      beta.sp.gamma[l,j] ~ dnorm(beta.comm.gamma[l] 
                                 + inprod(eff.traits.gamma[l,1:n.traits], traits[j,1:n.traits]) 
                                 # + eff.class.gamma[l,class_vec[j]]
                                 + eff.order.gamma[l,order_vec[j]]
                                 , tau.beta.comm.gamma[l])
    }
  }
  
  
  # Priors for Detection)
      for(l in 1:ncov.det){
        # rho.comm[l] ~ dt(0, 2.5, 1)
        rho.comm[l] ~ dnorm(0, 0.1)
        tau.rho.comm[l] ~ dgamma(1, 1)
        # sd.rho.comm[l] <- 1/sqrt(pow(tau.rho.comm[l])
        for(j in 1:nspec){
          rho.sp.tot[l,j] ~ dnorm(rho.comm[l], tau.rho.comm[l])
          tau.rho.sp.tot[l,j] <- tau.rho.comm[l]
          
          for(t in 1:nseason){
            rho.sp[l,j,t] ~ dnorm(rho.sp.tot[l,j], tau.rho.sp.tot[l,j])
          }
        }
      }
    
        
  # Ecological submodel: Define state conditional on parameters
  for(k in 1:nspec){       # loop over sites
    for(i in 1:nsite){ # loop over species
      
      # Occupancy/persistence/colonization varying across sites
      logit(psi1[i,k]) <- inprod(
        beta.sp.psi1[,k], 
        psi.cov[i,1,]
        ) 
      
      # State change probabilities
      for(t in 1:(nseason-1)){  # loop over seasons, but there should be one fewer state transitions than seasons
        # Persistence
        logit(phi[i,t,k]) <- inprod(beta.sp.phi[,k], psi.cov[i,t,])
          
        # and Colonization
        logit(gamma[i,t,k]) <- inprod(beta.sp.gamma[,k], psi.cov[i,t,])
      }
      
      
      # Initial conditions of system
      z[i,1,k] ~ dbern(psi1[i,k])   # presence/absence at start of study
      #z[i,1,k] ~ dbern(plogis(lpsi[k])) ??????? 
      # State transitions
      for(t in 2:nseason){
        z[i,t,k] ~ dbern(z[i, t-1, k] * phi[i,t-1,k] + (1-z[i, t-1, k]) * gamma[i,t-1,k])
      }
    }
  }
  
      
        
  # Observation model
  for(k in 1:nspec){       # Loop over species
    
    for(t in 1:nseason){   # Loop over seasons
      # rho0[k,t] ~ dnorm(rho0.tot[k], tau.rho0.tot[k])              # Detection
      # rho.sp[,k,t] ~ dnorm(rho.tot[k], tau.rho.tot[k]) 
    
      for(i in 1:nsite){         # Loop over sites
          # logit(p[i,t,k]) <- rho0[k,t] + rho1[k]*det1[i,t]       # Detection, with one covariate (total days sampled)
          logit(p[i,t,k]) <- inprod(rho.sp[,k,t],rho.cov[i,t,])
          
          y[i,t,k] ~ dbin((z[i,t,k]*p[i,t,k]), K[i,t])
        }
      }
  }
  
  # Derived parameters
    for(k in 1:nspec){    # loop over species
      
      # Species occupancy in each season
        for(i in 1:nsite){  
          #psi[1,k] <- psi1[,k]         # population occupancy
          psi[i,1,k] <- psi1[i,k]        # site occupancy
          
          for(t in 2:nseason){          # Loop over seasons
            psi[i,t,k] <- psi[i, t-1, k]*phi[i,t-1,k] + (1-psi[i, t-1, k])*gamma[i,t-1,k]
          }
        }
      
      # Number of occupied sites in each season
        n.occ[1,k] <- sum(z[,1,k])  
        for(t in 2:nseason){          # Loop over years
          n.occ[t,k] <- sum(z[,t,k])
          }
      
      
    }
    
    # Species richness and evenness in each season
      for(i in 1:nsite){
        for(t in 1:nseason){
          # Alpha diversity as local species richness (Hill Number 0)
          rich[i,t] <- sum(z[i,t,]) # Number of species occurring at each site, based on the occurrence matrix
          #rich[i,t] <- sum(psi[i,t,]) # Number of species occurring at each site, based on the occupancy probabilities
          
          # Species Diversity and Evenness (OPTIONAL)
            # for this, we will calculate Hill #1, then log-transform it to get Shannon Diversity
            # see Haight et al. 2023, Turrini and Knop 2015, or Boron et al 2019 to see how that works
            
            # calculating relative "abundances"
            # sum.psi[i,t] <- sum(psi[i,t,])  # sum of occupancy probabilities across species
            # sum.psi.checked[i,t] <- ifelse(sum.psi[i,t] == 0, 1E6, sum.psi[i,t]) # avoids dividing by 0 when calculating relative psi
            # 
            # for(k in 1:nspec){
            #   # relative.psi = relative occupancy: occupancy of each species divided by the across-species sum of probabilities
            #   relative.psi[i,t,k] <- psi[i,t,k]/sum.psi[i,t]
            #   log.relative.psi[i,t,k] <- ifelse(relative.psi[i,t,k]== 0,
            #                              log(relative.psi[i,t,k]+1E-6),
            #                              log(relative.psi[i,t,k]))
            # }
            # 
            # # Hill #1 (exponentiated Shannon Diversity)
            # hill1[i,t] <- exp(-sum(relative.psi[i,t,]*log.relative.psi[i,t,]))
            # #H[i,t] <- -sum(relative.psi[i,t,]*log.relative.psi[i,t,])   # Shannon Diversity Index
            # H[i,t] <- log(hill1[i,t])   # Shannon Diversity Index
            # 
            # # Evennness
            # #Pielou's Evenness = H/log(rich) or log(hill1)/log(rich)
            # #log.rich[i,t] <- ifelse(rich[i,t]== 0,
            # #                         log(rich[i,t]+1E-6),
            # #                         log(rich[i,t]))
            # log.sum.psi[i,t] <- ifelse(sum.psi[i,t]== 0,
            #                         log(sum.psi[i,t]+1E-6),
            #                         log(sum.psi[i,t]))
            # E[i,t] <- H[i,t]/log.sum.psi[i,t]    
        }
    
      }
      
  # Code for likelihood of each species observation outcome, used to calculate CPO
      for(i in 1:nsite){
        for(k in 1:nspec){ 
          for(t in 1:nseason){
            # A trick for calculating the binomial coefficient.
            BinCo[i,t,k] <- exp(logfact(K[i,t]) - (logfact(y[i,t,k]) + 
                                                     logfact(K[i,t] - y[i,t,k])))
            # the likelihood of observing what we did, given the model and the probabilities it estimated
            lik[i,t,k] <- ifelse(equals(y[i,t,k],0),   # if the species was not detected, then...
                                 # the likelihood of each non-observation = cumulative probability of 
                                 # (a) being in the study area, present at that site, but having been missed K times
                                 # (b) being in the study area, but not present at that site
                                 psi[i,t,k]*((1-p[i,t,k])^K[i,t]) + (1-psi[i,t,k]),
                                 # but if the species was detected, then likelihood = probability of that particular outcome: 
                                 # prob. of being in the study area, being present at the site, being detected a total of y times, 
                                 # but also not being detected on other survey occasions
                                 BinCo[i,t,k]*psi[i,t,k]*(p[i,t,k]^y[i,t,k]) * (1-p[i,t,k])^(K[i,t]-y[i,t,k])) 
          }
        }
      }
    
      
    
     
}