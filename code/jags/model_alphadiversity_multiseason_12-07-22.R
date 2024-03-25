model{
  for(i in 1:nsite){
    for(s in 1:nseason){
      # Priors for site-level random effects
      my_re[i,s] ~ dnorm(0, re_tau[s])
      
      log(mu_alpha[i,s]) <- inprod(  # log-normal model. Using inprod() allows for vector multiplication (more streamlined equation)
        beta[,s],                    # intercept and covariate slopes
        dm_alpha[i,,s]               # design matrix of the covariate data
      ) + my_re[i,s]                 # site-level random effect, i.e. the standard residual error term in log-linear models
      
      # site-level alpha diversity estimates are random parameters with 
      # standard deviation of alpha diversity estimates carried over directly from the occupancy model
      alpha_z[i,s] ~ dnorm(
        mu_alpha[i,s],
        pow(alpha_sd_known[i,s], -2)
      )
    
      
    }
  }
  
  # loop across seasons, but not over sites
  for(s in 1:nseason){
    # Priors for intercept and covariate slopes
    for(h in 1:npar_alpha){
      beta[h,s] ~ dnorm(0, 0.01)
    }
    
    # dispersion of random effect parameter, aka the residual standard error
    re_tau[s] ~ dgamma(1,1)
    re_sd[s] <- 1 / sqrt(re_tau[s])
  }
  
}