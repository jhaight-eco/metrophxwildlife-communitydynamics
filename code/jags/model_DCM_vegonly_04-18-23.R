model {
  
  # Specify priors (declare species levels effects as random, not dunif)
  for(k in 1:nspec){             # Loop over species
    # Intercepts
      # Without species traits effects
      #lpsi1[k] ~ dnorm(mu.lpsi1, tau.lpsi1)     # Initial occupancy
      #lphi[k] ~ dnorm(mu.lphi, tau.lphi)        # Persistence
      #lgamma[k] ~ dnorm(mu.lgamma, tau.lgamma)  # Colonization
      #lp[k] ~ dnorm(mu.lp, tau.lp)              # Detection
    
      # season-average detection probability
      lp.year[k] ~ dnorm(mu.lp, tau.lp)
      tau.lp.year[k] <- tau.lp  # the dispersion of detection for each species starts from 30
      
      
      # With species traits effects: 
      lpsi1[k] ~ dnorm(             # How common a species is varies based partly on their traits
        mu.lpsi1 + dint.psi.trait1 * trait1[k] + dint.psi.trait2 * trait2[k],     # 
        tau.lpsi1)
      
      lphi[k] ~ dnorm(             # Persistence: how likely a species is to persist at a site
        mu.lphi + dint.phi.trait1 * trait1[k]  + dint.phi.trait2 * trait2[k],
        tau.lphi)
      # Colonization: how likely a species is to colonize a site
      lgamma[k] ~ dnorm(
        mu.lgamma + dint.gamma.trait1*trait1[k]  + dint.gamma.trait2 * trait2[k], 
        tau.lgamma)
    
      # Slopes
      # Without species traits effects
      #beta1.psi[k] ~ dnorm(mu.beta1.psi, tau.mu.beta1.psi)
      
      # With species traits effects
      # # Species-specific slope for the first site-level covariate varies from the community-average based partly on the species trait(s)
      #beta1.psi[k] ~ dnorm(mu.beta1.psi + deff1.psi.trait1 * trait1[k] + deff1.psi.trait2 * trait2[k], tau.mu.beta1.psi)
      #beta1.phi[k] ~ dnorm(mu.beta1.phi + deff1.phi.trait1 * trait1[k] + deff1.phi.trait2 * trait2[k], tau.mu.beta1.phi)
      #beta1.gamma[k] ~ dnorm(mu.beta1.gamma + deff1.gamma.trait1 * trait1[k] + deff1.gamma.trait2 * trait2[k], tau.mu.beta1.gamma)
      
      #beta2.psi[k] ~ dnorm(mu.beta2.psi + deff2.psi.trait1 * trait1[k] + deff2.psi.trait2 * trait2[k], tau.mu.beta2.psi)
      #beta2.phi[k] ~ dnorm(mu.beta2.phi + deff2.phi.trait1 * trait1[k] + deff2.phi.trait2 * trait2[k], tau.mu.beta2.phi)
      #beta2.gamma[k] ~ dnorm(mu.beta2.gamma + deff2.gamma.trait1 * trait1[k] + deff2.gamma.trait2 * trait2[k], tau.mu.beta2.gamma)
      
      #beta3.psi[k] ~ dnorm(mu.beta3.psi + deff3.psi.trait1 * trait1[k] + deff3.psi.trait2 * trait2[k], tau.mu.beta3.psi)
      #beta3.phi[k] ~ dnorm(mu.beta3.phi + deff3.phi.trait1 * trait1[k] + deff3.phi.trait2 * trait2[k], tau.mu.beta3.phi)
      #beta3.gamma[k] ~ dnorm(mu.beta3.gamma + deff3.gamma.trait1 * trait1[k] + deff3.gamma.trait2 * trait2[k], tau.mu.beta3.gamma)
      
      beta4.psi[k] ~ dnorm(mu.beta4.psi + deff4.psi.trait1 * trait1[k] + deff4.psi.trait2 * trait2[k], tau.mu.beta4.psi)
      beta4.phi[k] ~ dnorm(mu.beta4.phi + deff4.phi.trait1 * trait1[k] + deff4.phi.trait2 * trait2[k], tau.mu.beta4.phi)
      beta4.gamma[k] ~ dnorm(mu.beta4.gamma + deff4.gamma.trait1 * trait1[k] + deff4.gamma.trait2 * trait2[k], tau.mu.beta4.gamma)
      
      alpha1[k] ~ dnorm(mu.alpha1, tau.mu.alpha1)
      
  }
  
  # Specify priors for species covariate effects 
      # Effects of two species traits on intercepts of occupancy, persistence, and colonization
      dint.psi.trait1 ~ dnorm(0,0.1)
      dint.psi.trait2 ~ dnorm(0,0.1)
      dint.phi.trait1 ~ dnorm(0,0.1)
      dint.phi.trait2 ~ dnorm(0,0.1)
      dint.gamma.trait1 ~ dnorm(0,0.1)
      dint.gamma.trait2 ~ dnorm(0,0.1)
            
      # Trait interactions with the effect of covariate 1 (urbanization) on occupancy, persistence, and colonization
      #deff1.psi.trait1 ~ dnorm(0,0.1)
      #deff1.psi.trait2 ~ dnorm(0,0.1)
      #deff1.phi.trait1 ~ dnorm(0,0.1)
      #deff1.phi.trait2 ~ dnorm(0,0.1) 
      #deff1.gamma.trait1 ~ dnorm(0,0.1)
      #deff1.gamma.trait2 ~ dnorm(0,0.1)
      
      # Trait interactions with the effect of covariate 2 (quadratic urbanization) on occupancy, persistence, and colonization
      #deff2.psi.trait1 ~ dnorm(0,0.1)
      #deff2.psi.trait2 ~ dnorm(0,0.1)
      #deff2.phi.trait1 ~ dnorm(0,0.1)
      #deff2.phi.trait2 ~ dnorm(0,0.1) 
      #deff2.gamma.trait1 ~ dnorm(0,0.1)
      #deff2.gamma.trait2 ~ dnorm(0,0.1)
      
      # Trait interactions with the effect of covariate 2 (quadratic urbanization) on occupancy, persistence, and colonization
      #deff3.psi.trait1 ~ dnorm(0,0.1)
      #deff3.psi.trait2 ~ dnorm(0,0.1)
      #deff3.phi.trait1 ~ dnorm(0,0.1)
      #deff3.phi.trait2 ~ dnorm(0,0.1) 
      #deff3.gamma.trait1 ~ dnorm(0,0.1)
      #deff3.gamma.trait2 ~ dnorm(0,0.1)
      
      # Trait interactions with the effect of covariate 3 (vegetation) on occupancy, persistence, and colonization
      deff4.psi.trait1 ~ dnorm(0,0.1)
      deff4.psi.trait2 ~ dnorm(0,0.1)
      deff4.phi.trait1 ~ dnorm(0,0.1)
      deff4.phi.trait2 ~ dnorm(0,0.1) 
      deff4.gamma.trait1 ~ dnorm(0,0.1)
      deff4.gamma.trait2 ~ dnorm(0,0.1)
  
  # Specify hyperpriors (priors for community-average hyperparameters)
      # Intercepts
        mu.lpsi1 <- logit(mean.psi1)         # Initial occupancy
        mean.psi1 ~ dunif(0,1)
        tau.lpsi1 <- pow(sd.lpsi1, -2) 
        sd.lpsi1 ~ dunif(0, 10) 
        mu.lphi <- logit(mean.phi) # Persistence 
        mean.phi ~ dunif(0, 1) 
        tau.lphi <- pow(sd.lphi, -2) 
        sd.lphi ~ dunif(0, 10) 
        mu.lgamma <- logit(mean.gamma) # Colonization 
        mean.gamma ~ dunif(0, 1) 
        tau.lgamma <- pow(sd.lgamma, -2) 
        sd.lgamma ~ dunif(0, 10) 
        mu.lp <- logit(mean.p) # Detection varies among season
        mean.p ~ dunif(0, 1) 
        tau.lp <- pow(sd.lp, -2) 
        sd.lp ~ dunif(0, 10)
    
      # Slopes
        # Site covariate 1
        #mu.beta1.psi ~ dnorm(0, 0.1)
        #tau.mu.beta1.psi <- pow(sd.mu.beta1.psi, -2) 
        #sd.mu.beta1.psi ~ dunif(0, 10) 
        #mu.beta1.phi ~ dnorm(0, 0.1)
        #tau.mu.beta1.phi <- pow(sd.mu.beta1.phi, -2) 
        #sd.mu.beta1.phi ~ dunif(0, 10) 
        #mu.beta1.gamma ~ dnorm(0, 0.1)
        #tau.mu.beta1.gamma <- pow(sd.mu.beta1.gamma, -2) 
        #sd.mu.beta1.gamma ~ dunif(0, 10) 
      
        # Site covariate 2
        #mu.beta2.psi ~ dnorm(0, 0.1)
        #tau.mu.beta2.psi <- pow(sd.mu.beta2.psi, -2) 
        #sd.mu.beta2.psi ~ dunif(0, 10) 
        #mu.beta2.phi ~ dnorm(0, 0.1)
        #tau.mu.beta2.phi <- pow(sd.mu.beta2.phi, -2) 
        #sd.mu.beta2.phi ~ dunif(0, 10) 
        #mu.beta2.gamma ~ dnorm(0, 0.1)
        #tau.mu.beta2.gamma <- pow(sd.mu.beta2.gamma, -2) 
        #sd.mu.beta2.gamma ~ dunif(0, 10) 
        
        # Site covariate 3
        #mu.beta3.psi ~ dnorm(0, 0.1)
        #tau.mu.beta3.psi <- pow(sd.mu.beta3.psi, -2) 
        #sd.mu.beta3.psi ~ dunif(0, 10) 
        #mu.beta3.phi ~ dnorm(0, 0.1)
        #tau.mu.beta3.phi <- pow(sd.mu.beta3.phi, -2) 
        #sd.mu.beta3.phi ~ dunif(0, 10) 
        #mu.beta3.gamma ~ dnorm(0, 0.1)
        #tau.mu.beta3.gamma <- pow(sd.mu.beta3.gamma, -2) 
        #sd.mu.beta3.gamma ~ dunif(0, 10) 
        
        # Site covariate 4 (varies by season)
        mu.beta4.psi ~ dnorm(0, 0.1)
        tau.mu.beta4.psi <- pow(sd.mu.beta4.psi, -2) 
        sd.mu.beta4.psi ~ dunif(0, 10) 
        mu.beta4.phi ~ dnorm(0, 0.1)
        tau.mu.beta4.phi <- pow(sd.mu.beta4.phi, -2) 
        sd.mu.beta4.phi ~ dunif(0, 10) 
        mu.beta4.gamma ~ dnorm(0, 0.1)
        tau.mu.beta4.gamma <- pow(sd.mu.beta4.gamma, -2) 
        sd.mu.beta4.gamma ~ dunif(0, 10) 
      
        # Detection covariate
        mu.alpha1 ~ dnorm(0, 0.1)
        tau.mu.alpha1 <- pow(sd.mu.alpha1, -2) 
        sd.mu.alpha1 ~ dunif(0, 10) 
  
  # Ecological submodel: Define state conditional on parameters
  for(k in 1:nspec){       # loop over sites
    
    # placing these linear predictors here assumes they are consistent across sites. Otherwise, move down below
    #logit(phi[k]) <- lphi[k]     # Persistence
    #logit(gamma[k]) <- lgamma[k] # Colonization
    
    for(i in 1:nsite){ # loop over species
      
      # With occupancy/persistence/colonization varies across sites, they need to be down here
      logit(psi1[i,k]) <- lpsi1[k] # Initial occupancy
      #+ beta1.psi[k]*cov.occ1[i]
      #+ beta2.psi[k]*cov.occ2[i]
      #+ beta3.psi[k]*cov.occ3[i]
      + beta4.psi[k]*cov.occ4[i,1]
      
      # Estimates of phi and gamma in the first season should probably be ignored, 
      # but they're necessary for having a time-varying covariate on persistence and colonization
      for(t in 1:(nseason-1)){  # loop over seasons, but there should be one fewer state transitions than seasons
        logit(phi[i,t,k]) <- lphi[k]     # Persistence
        #+ beta1.phi[k]*cov.occ1[i] 
        #+ beta2.phi[k]*cov.occ2[i]
        #+ beta3.phi[k]*cov.occ3[i]
        + beta4.phi[k]*cov.occ4[i,t]
        
        logit(gamma[i,t,k]) <- lgamma[k] # Colonization
        #+ beta1.gamma[k]*cov.occ1[i] 
        #+ beta2.gamma[k]*cov.occ2[i]
        #+ beta3.gamma[k]*cov.occ3[i]
        + beta4.gamma[k]*cov.occ4[i,t]
      }
        
      
      # Initial conditions of system
      z[i,1,k] ~ dbern(psi1[i,k])   # presence/absence at start of study
      #z[i,1,k] ~ dbern(plogis(lpsi[k])) ??????? 
      # State transitions
      for(t in 2:nseason){
        z[i,t,k] ~ dbern(z[i, t-1, k] * phi[i,t-1,k] + (1-z[i, t-1, k]) * gamma[i,t-1,k])
        # SHOULD THIS BE BASED ON THE PERSISTENCE/COLONIZATION PROBABILITIES IN THE BEFORE SEASON OR THE AFTER SEASON?
        # MAYBE PROBABILITY OF PERSISTING BASED ON BEFORE CONDITIONS AND COLONIZING BASED ON AFTER CONDITIONS?
      }
    }
  }
  
      
        
  # Observation model
  for(k in 1:nspec){       # Loop over species
    
    for(t in 1:nseason){   # Loop over seasons
      
        lp[k,t] ~ dnorm(lp.year[k], tau.lp.year[k])              # Detection
    
      for(i in 1:nsite){         # Loop over sites
        #logit(p[i,k]) <- lp[k] + alpha1[k]*det1[i]       # Detection, with one covariate (total days sampled)
        
          #mu.p[i,k] <- p[i,k]*z[i,t,k] 
          #mu.p[i,j,r] <- p[i,j,r]*Z[i,j,r] # from the UWIN code   
          logit(p[i,t,k]) <- lp[k,t] + alpha1[k]*det1[i,t]       # Detection, with one covariate (total days sampled)
          
          y[i,t,k] ~ dbin((z[i,t,k]*p[i,t,k]), K[i,t])
          #y[i,t,k] ~ dbin(mu.p[i,t,k], K[i,t]) 
          
          #for(j in 1:nsurvey){   # Loop over surveys
          #  y[i,j,t,k] ~ dbern(z[i,t,k]*p[i,k])      
           # y[i,j,t,k] ~ dbern(mu.p[i,t, k]) 
          #}
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
          # Alpha species diversity (local species richness = Hill Number 0)
          rich[i,t] <- sum(z[i,t,]) # Number of species occurring at each site, based on the occurrence matrix
          #rich[i,t] <- sum(psi[i,t,]) # Number of species occurring at each site, based on the occupancy probabilities
          
          # Species Diversity and Evenness
          # for this, we will essentially have to calculate Hill #1, log() that to get Shannon Diversity
          # see the UWIN code and Turrini and Knop 2015 or Boron et al 2019 to see how that works
          
          # calculating relative "abundances"
          sum.psi[i,t] <- sum(psi[i,t,])  # sum of occupancy probabilities across species
          sum.psi.checked[i,t] <- ifelse(sum.psi[i,t] == 0, 1E6, sum.psi[i,t]) # avoids dividing by 0 when calculating relative psi
          
          for(k in 1:nspec){
            # relative.psi = relative occupancy: occupancy of each species divided by the across-species sum of probabilities
            relative.psi[i,t,k] <- psi[i,t,k]/sum.psi[i,t]
            log.relative.psi[i,t,k] <- ifelse(relative.psi[i,t,k]== 0,
                                       log(relative.psi[i,t,k]+1E-6),
                                       log(relative.psi[i,t,k]))
          }
          
          # Hill #1 (exponentiated Shannon Diversity)
          hill1[i,t] <- exp(-sum(relative.psi[i,t,]*log.relative.psi[i,t,]))
          #H[i,t] <- -sum(relative.psi[i,t,]*log.relative.psi[i,t,])   # Shannon Diversity Index
          H[i,t] <- log(hill1[i,t])   # Shannon Diversity Index
          
          # Evennness
          #Pielou's Evenness = H/log(rich) or log(hill1)/log(rich)
          #log.rich[i,t] <- ifelse(rich[i,t]== 0,
          #                         log(rich[i,t]+1E-6),
          #                         log(rich[i,t]))
          log.sum.psi[i,t] <- ifelse(sum.psi[i,t]== 0,
                                  log(sum.psi[i,t]+1E-6),
                                  log(sum.psi[i,t]))
          E[i,t] <- H[i,t]/log.sum.psi[i,t]    
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