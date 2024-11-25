#### 
# Haight, Jeffrey D.
# #.#

#### Setup ####
rm(list=ls()) # clear the environment
gc()
set.seed(4321)

# Set working directory
setwd("~/GitHub/metrophxwildlife-communitydynamics")

  # Load necessary packages
  #library(dplyr)
  library(tidyverse)
  
  # modeling packages
  library(jagsUI)
  library(runjags)
  
  # visualization packages
  library(flextable)
  library(ggplot2)
  library(gghighlight)
  library(viridis)
  library(RColorBrewer)
  library(lterpalettefinder)   # devtools::install_github("lter/lterpalettefinder")
  library(png)

  # import some color palettes
  load("./data/model_inputs/arizonapalettes.RData")
  # Warm-Dry: color = 2, fill = 3
  # Warm-Wet: color = 4, fill = 5
  # Cool-Wet: color = 6, fill = 1

##### Import Data #####
  # The main dataset containing all the site and regional covariates 
  # and summarized community composition parameters
  data <- read.csv("./data/model_outputs/data_PhoenixCommDynamics.csv")
  
  # import the meta-analysis models and their parameter summaries
  m.sr <- readRDS("./data/model_outputs/metaanalysislogglm_speciesrichness_30k.rds")
  
  m.sr.sum <- read.csv("./data/model_summaries/communitymetaanalysisresults_richness.csv")
  
  # get the mcmc steps all in one matrix (for each season)
  str(m.sr$sims.list)
  m.sr.mcmc.s1 <- cbind(
    m.sr$sims.list$beta[,,1],
    m.sr$sims.list$re_sd[,1],
    m.sr$sims.list$deviance)
  
  m.sr.mcmc.s2 <- cbind(
    m.sr$sims.list$beta[,,2],
    m.sr$sims.list$re_sd[,2],
    m.sr$sims.list$deviance)
  
  m.sr.mcmc.s3 <- cbind(
    m.sr$sims.list$beta[,,3],
    m.sr$sims.list$re_sd[,3],
    m.sr$sims.list$deviance)
  
  # rename the columns
  colnames(m.sr.mcmc.s1) <- c("int", 
                              "urb",
                              "urbquad",
                              "pdiv", 
                              # "veg", 
                              "re", "dev")
  colnames(m.sr.mcmc.s2) <- c("int", 
                              "urb",
                              "urbquad",
                              "pdiv", 
                              # "veg", 
                              "re", "dev")
  colnames(m.sr.mcmc.s3) <- c("int", 
                              "urb",
                              "urbquad",
                              "pdiv", 
                              # "veg", 
                              "re", "dev")
  
  
  colnames(m.sr.sum)[1] <- c("Predictor")
  
  m.sr.sum

  
  # check where the peak/valley of the urbanization curve would theoretically be in each season
  (mean(m.sr.mcmc.s1[,2])/mean(m.sr.mcmc.s1[,3])*-1)*sd(data$imp1km)+mean(data$imp1km)  # effectively straight
  (mean(m.sr.mcmc.s2[,2])/mean(m.sr.mcmc.s2[,3])*-1)*sd(data$imp1km)+mean(data$imp1km)  # very, very slightly convex
  (mean(m.sr.mcmc.s3[,2])/mean(m.sr.mcmc.s3[,3])*-1)*sd(data$imp1km)+mean(data$imp1km)  # also very, very slightly convex
  # as we can see, there is no real peak,  despite the inclusion of the quadratic term
  
#### Predict Richness at Select Urbanization Levels
  urb.sel <- seq(0, 80, 10)
  urb.sel.std <- (urb.sel - mean(data$imp1km))/sd(data$imp1km)  # standardized values from 0 to 80% urbanization
    
  
  dm.urb.sel <- cbind(
    1, # intercept
    urb.sel.std, #urban
    urb.sel.std^2, # urbanization quadratic term
    0, # patch diversity
    0  # veg productivity
  )
  
  preds.sr.urb1 <- exp(m.sr.mcmc.s1[,1:5] %*% t(dm.urb.sel))
  preds.sr.urb1.cri <- apply(
    preds.sr.urb1,
    2,
    quantile, 
    probs = c(0.025,0.5,0.975)
  ) %>% t() %>% data.frame()
  preds.sr.urb1.mean <- apply(preds.sr.urb1, 2,  function(x) mean(x, na.rm = TRUE))
  
  preds.sr.urb2 <- exp(m.sr.mcmc.s2[,1:5] %*% t(dm.urb.sel))
  preds.sr.urb2.cri <- apply(
    preds.sr.urb2,
    2,
    quantile, 
    probs = c(0.025,0.5,0.975)
  ) %>% t() %>% data.frame()
  preds.sr.urb2.mean <- apply(preds.sr.urb2, 2,  function(x) mean(x, na.rm = TRUE))
  
  preds.sr.urb3 <- exp(m.sr.mcmc.s3[,1:5] %*% t(dm.urb.sel))
  preds.sr.urb3.cri <- apply(
    preds.sr.urb3,
    2,
    quantile, 
    probs = c(0.025,0.5,0.975)
  ) %>% t() %>% data.frame()
  preds.sr.urb3.mean <- apply(preds.sr.urb3, 2,  function(x) mean(x, na.rm = TRUE))
  
  
  plot(urb.sel, preds.sr.urb1.mean, ylim = c(0,12))
  plot(urb.sel, preds.sr.urb2.mean, ylim = c(0,12))
  plot(urb.sel, preds.sr.urb3.mean, ylim = c(0,12))
  
  # Where does the steepest slope lie?
    # Warm-Dry
    (preds.sr.urb1.mean[2] - preds.sr.urb1.mean[1])/10  # 0 to 10
    (preds.sr.urb1.mean[3] - preds.sr.urb1.mean[2])/10  # 10 to 20
    (preds.sr.urb1.mean[4] - preds.sr.urb1.mean[3])/10  # 20 to 30
    (preds.sr.urb1.mean[5] - preds.sr.urb1.mean[4])/10  # 30 to 40
    (preds.sr.urb1.mean[6] - preds.sr.urb1.mean[5])/10  # 40 to 50
    (preds.sr.urb1.mean[7] - preds.sr.urb1.mean[6])/10  # 50 to 60
    (preds.sr.urb1.mean[8] - preds.sr.urb1.mean[7])/10  # 60 to 70
    (preds.sr.urb1.mean[9] - preds.sr.urb1.mean[8])/10  # 70 to 80
    # Warm-Wet
    (preds.sr.urb2.mean[2] - preds.sr.urb2.mean[1])/10  # 0 to 10
    (preds.sr.urb2.mean[3] - preds.sr.urb2.mean[2])/10  # 10 to 20
    (preds.sr.urb2.mean[4] - preds.sr.urb2.mean[3])/10  # 20 to 30
    (preds.sr.urb2.mean[5] - preds.sr.urb2.mean[4])/10  # 30 to 40
    (preds.sr.urb2.mean[6] - preds.sr.urb2.mean[5])/10  # 40 to 50
    (preds.sr.urb2.mean[7] - preds.sr.urb2.mean[6])/10  # 50 to 60
    (preds.sr.urb2.mean[8] - preds.sr.urb2.mean[7])/10  # 60 to 70
    (preds.sr.urb2.mean[9] - preds.sr.urb2.mean[8])/10  # 70 to 80
    # Cool-Wet
    (preds.sr.urb3.mean[2] - preds.sr.urb3.mean[1])/10  # 0 to 10
    (preds.sr.urb3.mean[3] - preds.sr.urb3.mean[2])/10  # 10 to 20
    (preds.sr.urb3.mean[4] - preds.sr.urb3.mean[3])/10  # 20 to 30
    (preds.sr.urb3.mean[5] - preds.sr.urb3.mean[4])/10  # 30 to 40
    (preds.sr.urb3.mean[6] - preds.sr.urb3.mean[5])/10  # 40 to 50
    (preds.sr.urb3.mean[7] - preds.sr.urb3.mean[6])/10  # 50 to 60
    (preds.sr.urb3.mean[8] - preds.sr.urb3.mean[7])/10  # 60 to 70
    (preds.sr.urb3.mean[9] - preds.sr.urb3.mean[8])/10  # 70 to 80
    
  
  
  # richness values at zero impervious surface cover 
  preds.sr.urb1.mean[1]; preds.sr.urb1.cri[1,]
  preds.sr.urb2.mean[1]; preds.sr.urb2.cri[1,]
  preds.sr.urb3.mean[1]; preds.sr.urb3.cri[1,]
  
  # richness at 30% impervious surface cover
  preds.sr.urb1.mean[4]; preds.sr.urb1.cri[4,]
  preds.sr.urb2.mean[4]; preds.sr.urb2.cri[4,]
  preds.sr.urb3.mean[4]; preds.sr.urb3.cri[4,]
  
#### Main Figure: Richness Predicted Across Gradients 
    npred <- 200    # how many values to predict across
    
    cov1 <- data$imp1km
    cov1.std <- data$imp1km_std
    
    cov3 <- data$sdhi5km
    cov3.std <- data$sdhi5km_std
    
    #cov4 <- cbind(data$ndvi_wd1km, data$ndvi_ww1km, data$ndvi_cw1km)
    #cov4.std <- cbind(data$ndvi_wd1km_std, data$ndvi_ww1km_std, data$ndvi_cw1km_std)
  
# 1) Predict across Urbanization Gradient
  # Values to predict off of
  # Hold all other covariates constant at their mean (~0 since they are standardized)
  dm.urb <- cbind(
    # Season 1
    1,# intercept
    seq(min(cov1.std),max(cov1.std), length.out = npred), #urban
    seq(min(cov1.std),max(cov1.std), length.out = npred)^2, # urbanization quadratic term
    0#, # patch diversity
    #0  # veg productivity. This was standardized by season, so we can use the same value here across all seasons
  )
  
  preds.sr.urb1 <- exp(m.sr.mcmc.s1[,1:ncol(dm.urb)] %*% t(dm.urb))   # exponentiate, because that's the reverse of the log-link
  preds.sr.urb1 <- apply(
    preds.sr.urb1,
    2,
    quantile, 
    probs = c(0.025,0.5,0.975)
  ) %>% t() %>% data.frame()
  colnames(preds.sr.urb1)[1:3] <- c("lower95", "median", "upper95")
  # this is not a season-varying covariate, but the 
  
  # add a column for the covariate
  preds.sr.urb1$impervious <- seq(min(cov1), max(cov1), length.out = npred)
  preds.sr.urb1$impervious_std <- seq(min(cov1.std), max(cov1.std), length.out = npred)
  str(preds.sr.urb1)
  
  # repeat for each season, using the different betas in each season
    # Season 2: Warm-Wet
    preds.sr.urb2 <- exp(m.sr.mcmc.s2[,1:ncol(dm.urb)] %*% t(dm.urb))   # note that it is now m.sr.mcmc.s2
    preds.sr.urb2 <- apply(
      preds.sr.urb2,
      2,
      quantile, 
      probs = c(0.025,0.5,0.975)
    ) %>% t() %>% data.frame()
    colnames(preds.sr.urb2)[1:3] <- c("lower95", "median", "upper95")
    
    preds.sr.urb2$impervious <- seq(min(cov1), max(cov1), length.out = npred)
    preds.sr.urb2$impervious_std <- seq(min(cov1.std), max(cov1.std), length.out = npred)
    str(preds.sr.urb2)
    
    # Season 3: Cool-Wet
    preds.sr.urb3 <- exp(m.sr.mcmc.s3[,1:ncol(dm.urb)] %*% t(dm.urb))   # m.sr.mcmc.s3
    preds.sr.urb3 <- apply(
      preds.sr.urb3,
      2,
      quantile, 
      probs = c(0.025,0.5,0.975)
    ) %>% t() %>% data.frame()
    colnames(preds.sr.urb3)[1:3] <- c("lower95", "median", "upper95")
    
    preds.sr.urb3$impervious <- seq(min(cov1), max(cov1), length.out = npred)
    preds.sr.urb3$impervious_std <- seq(min(cov1.std), max(cov1.std), length.out = npred)
    str(preds.sr.urb3)
    
  

  
# 2) Predict across Patch Diversity Gradient
    dm.pdiv <- cbind(
      # Season 1
      1,# intercept
      0, #urban
      0, # urbanization quadratic term
      seq(min(cov3.std),max(cov3.std), length.out = npred)#, # patch diversity
      #0  # veg productivity. This was standardized by season, so we can use the same value here across all seasons
    )
    
    preds.sr.pdiv1 <- exp(m.sr.mcmc.s1[,1:ncol(dm.pdiv)] %*% t(dm.pdiv))   # exponentiate, because that's the reverse of the log-link
    preds.sr.pdiv1 <- apply(
      preds.sr.pdiv1,
      2,
      quantile, 
      probs = c(0.025,0.5,0.975)
    ) %>% t() %>% data.frame()
    colnames(preds.sr.pdiv1)[1:3] <- c("lower95", "median", "upper95")
    # this is not a season-varying covariate, but the 
    
    # add a column for the covariate
    preds.sr.pdiv1$sdhi <- seq(min(cov3), max(cov3), length.out = npred)
    preds.sr.pdiv1$sdhi_std <- seq(min(cov3.std), max(cov3.std), length.out = npred)
    str(preds.sr.pdiv1)
    
    # repeat for each season, using the different betas in each season
    # Season 2: Warm-Wet
    preds.sr.pdiv2 <- exp(m.sr.mcmc.s2[,1:ncol(dm.pdiv)] %*% t(dm.pdiv))   # note that it is now m.sr.mcmc.s2
    preds.sr.pdiv2 <- apply(
      preds.sr.pdiv2,
      2,
      quantile, 
      probs = c(0.025,0.5,0.975)
    ) %>% t() %>% data.frame()
    colnames(preds.sr.pdiv2)[1:3] <- c("lower95", "median", "upper95")
    
    preds.sr.pdiv2$sdhi <- seq(min(cov3), max(cov3), length.out = npred)
    preds.sr.pdiv2$sdhi_std <- seq(min(cov3.std), max(cov3.std), length.out = npred)
    str(preds.sr.pdiv2)
    
    # Season 3: Cool-Wet
    preds.sr.pdiv3 <- exp(m.sr.mcmc.s3[,1:ncol(dm.pdiv)] %*% t(dm.pdiv))   # m.sr.mcmc.s3
    preds.sr.pdiv3 <- apply(
      preds.sr.pdiv3,
      2,
      quantile, 
      probs = c(0.025,0.5,0.975)
    ) %>% t() %>% data.frame()
    colnames(preds.sr.pdiv3)[1:3] <- c("lower95", "median", "upper95")
    
    preds.sr.pdiv3$sdhi <- seq(min(cov3), max(cov3), length.out = npred)
    preds.sr.pdiv3$sdhi_std <- seq(min(cov3.std), max(cov3.std), length.out = npred)
    str(preds.sr.pdiv3)

# 3) Predict across Vegetation Productivity Gradients
  #dm.veg1 <- cbind(
    # Season 1
  #  1,# intercept
  #  mean(cov1.std), #urban
  #  mean(cov1.std)^2, # urban quadratic term
  #  mean(cov3.std), # patch diversity
  #  seq(min(cov4.std[,1]),max(cov4.std[,1]), length.out = npred)  # veg productivity
  #)
  #dm.veg2 <- cbind(
    # Season 2
  #  1,# intercept
  #  mean(cov1.std), #urban
  #  mean(cov1.std)^2, # urban quadratic term
  #  mean(cov3.std), # patch diversity
  #  seq(min(cov4.std[,2]),max(cov4.std[,2]), length.out = npred)  # veg productivity
  #)
  #dm.veg3 <- cbind(
    # Season 3
  #  1,# intercept
  #  mean(cov1.std), #urban
  #  mean(cov1.std)^2, # urban quadratic term
  #  mean(cov3.std), # patch diversity
  #  seq(min(cov4.std[,3]),max(cov4.std[,3]), length.out = npred)  # veg productivity
  #)
  
  # predict for the first season's NDVI values
  #    preds.sr.veg1 <- exp(m.sr.mcmc.s1[,1:5] %*% t(dm.veg1))   # exponentiate, because that's the reverse of the log-link
  #    preds.sr.veg1 <- apply(
  #      preds.sr.veg1,
  #      2,
  #      quantile, 
  #      probs = c(0.025,0.5,0.975)
  #    ) %>% t() %>% data.frame()
  #    colnames(preds.sr.veg1)[1:3] <- c("lower95", "median", "upper95")
      
      # add a column for the covariate
  #    preds.sr.veg1$ndvi <- seq(min(cov4[,1]), max(cov4[,1]), length.out = npred)
  #    preds.sr.veg1$ndvi_std <- seq(min(cov4.std[,1]), max(cov4.std[,1]), length.out = npred)
  #    str(preds.sr.veg1)
    
  # repeat for the other seasons
      # Season 2
   #   preds.sr.veg2 <- exp(m.sr.mcmc.s2[,1:5] %*% t(dm.veg2))   # exponentiate, because that's the reverse of the log-link
  #    preds.sr.veg2 <- apply(
  #      preds.sr.veg2,
  #      2,
  #      quantile, 
  #      probs = c(0.025,0.5,0.975)
  #    ) %>% t() %>% data.frame()
  #    colnames(preds.sr.veg2)[1:3] <- c("lower95", "median", "upper95")
  #    preds.sr.veg2$ndvi <- seq(min(cov4[,2]), max(cov4[,2]), length.out = npred)
   ##   preds.sr.veg2$ndvi_std <- seq(min(cov4.std[,2]), max(cov4.std[,2]), length.out = npred)
      #str(preds.sr.veg2)
      
      # Season 3
  #    preds.sr.veg3 <- exp(m.sr.mcmc.s3[,1:5] %*% t(dm.veg3))   # exponentiate, because that's the reverse of the log-link
  #    preds.sr.veg3 <- apply(
  #      preds.sr.veg3,
  #      2,
  #      quantile, 
  #      probs = c(0.025,0.5,0.975)
  #    ) %>% t() %>% data.frame()
  #    colnames(preds.sr.veg3)[1:3] <- c("lower95", "median", "upper95")
  #    preds.sr.veg3$ndvi <- seq(min(cov4[,3]), max(cov4[,3]), length.out = npred)
  #    preds.sr.veg3$ndvi_std <- seq(min(cov4.std[,3]), max(cov4.std[,3]), length.out = npred)
      #str(preds.sr.veg3)


# reshape the data for easier plotting
data.long <- rbind(
  data %>% mutate(season = "1_wd") %>% select(site, season, imp1km, sdhi5km, "ndvi" = ndvi_wd1km, "rich_mean" = rich_wd_mean, "rich_2_5"= rich_wd_2_5, "rich_97_5"= rich_wd_97_5),
  data %>% mutate(season = "2_ww") %>% select(site, season, imp1km, sdhi5km, "ndvi" = ndvi_ww1km, "rich_mean"= rich_ww_mean, "rich_2_5"= rich_ww_2_5, "rich_97_5"= rich_ww_97_5),
  data %>% mutate(season = "3_cw") %>% select(site, season, imp1km, sdhi5km, "ndvi" = ndvi_cw1km, "rich_mean"= rich_cw_mean, "rich_2_5"= rich_cw_2_5, "rich_97_5"= rich_cw_97_5)
) #%>% str()


preds.sr.urb <- rbind(
  preds.sr.urb1 %>% mutate(season = "1_wd"), 
  preds.sr.urb2 %>% mutate(season = "2_ww"), 
  preds.sr.urb3 %>% mutate(season = "3_cw"))  #  %>% str()
      
##### Plot Actual Richness Estimates Across Urbanization #####
      palette_ggdemo(wildflower)
      
      
      # poppy orange and brittlebush yellow
      wildflower[c(16,18)] 
      
      plot.wd <- ggplot() +
        geom_pointrange(data = data, aes(x = imp1km, y = rich_wd_mean, ymin = rich_wd_2_5, ymax = rich_wd_97_5),
                   col = "#ff8701", fatten = 1)+
        stat_smooth(data = data, aes(x = imp1km, y = rich_wd_mean), fill = "#f2ff0d", col = "#ff8701",
                    method = "lm", formula = y ~ x + I(x^2), se = FALSE, size = 0.5) +
        #geom_smooth(data = preds.sr.urb1, aes(x = impervious, y = median, ymin = lower95, ymax = upper95,),
        #            fill = "#f2ff0d", color = "#ff8701", se = FALSE, linetype = "dashed")+
        
        theme_bw() + 
        coord_cartesian(xlim=c(0, max(data$imp1km)), ylim=c(0,13))+
        labs(x = "Urbanization \n(% Impervious Surface)", y = "Species Richness")  +
        theme(axis.text.x = element_text(face = "bold"), 
              axis.text.y = element_text(face = "bold"), 
              axis.ticks = element_blank(),
              axis.title.x = element_text(face = "bold"), 
              axis.title.y = element_text(face = "bold"),
        ) 
      plot.wd
      
      
      # prickly pear pinks
      wildflower[c(4,5)]
      
      
      plot.ww <- ggplot() +
        geom_pointrange(data = data, aes(x = imp1km, y = rich_ww_mean, ymin = rich_ww_2_5, ymax = rich_ww_97_5),
                        col = "#d10073", fatten = 1)+
        stat_smooth(data = data, aes(x = imp1km, y = rich_ww_mean), fill = "#ff23c1", col = "#d10073",
                    method = "lm", formula = y ~ x + I(x^2), se = FALSE, size = 0.5) +
        #geom_smooth(data = preds.sr.urb1, aes(x = impervious, y = median, ymin = lower95, ymax = upper95,),
        #            fill = "#ff23c1", color = "#d10073", se = FALSE, linetype = "dashed")+
        
        
        theme_bw() + 
        coord_cartesian(xlim=c(0, max(data$imp1km)), ylim=c(0,13))+
        labs(x = "Urbanization \n(% Impervious Surface)", y = "Species Richness")  +
        theme(axis.text.x = element_text(face = "bold"), 
              axis.text.y = element_text(face = "bold"), 
              axis.ticks = element_blank(),
              axis.title.x = element_text(face = "bold"), 
              axis.title.y = element_text(face = "bold"),
        ) 
      
      plot.ww
      
      wildflower[c(20,21)] 
      
      plot.cw <- ggplot() +
        geom_pointrange(data = data, aes(x = imp1km, y = rich_cw_mean, ymin = rich_cw_2_5, ymax = rich_cw_97_5),
                        col = "#2f3ca3", fatten = 1)+
        stat_smooth(data = data, aes(x = imp1km, y = rich_cw_mean), fill = "#7291fb", col = "#2f3ca3",
                    method = "lm", formula = y ~ x + I(x^2), se = FALSE, size = 0.5) +
        #geom_smooth(data = preds.sr.urb1, aes(x = impervious, y = median, ymin = lower95, ymax = upper95,),
        #            fill = "#7291fb", color = "#2f3ca3", se = FALSE, linetype = "dashed")+   # predict
        
        theme_bw() + 
        coord_cartesian(xlim=c(0, max(data$imp1km)), ylim=c(0,18))+
        labs(x = "Urbanization \n(% Impervious Surface)", y = "Species Richness")  +
        theme(axis.text.x = element_text(face = "bold"), 
              axis.text.y = element_text(face = "bold"), 
              axis.ticks = element_blank(),
              axis.title.x = element_text(face = "bold"), 
              axis.title.y = element_text(face = "bold"),
        ) 
      
      plot.cw
      
      ggsave("./figures/richness_actual_vs_urb_1warmdry.png",
             plot.wd,
             width = 3,
             height = 2.5,
             units = "in",
             dpi = 300)
      
      ggsave("./figures/richness_actual_vs_urb_2warmwet.png",
             plot.ww,
             width = 3,
             height = 2.5,
             units = "in",
             dpi = 300)
      
      ggsave("./figures/richness_actual_vs_urb_3coolwet.png",
             plot.cw,
             width = 3,
             height = 2.5,
             units = "in",
             dpi = 300)
      
      
      
      plot.all <- ggplot() +
        geom_pointrange(data = data, aes(x = imp1km, y = rich_wd_mean, ymin = rich_wd_2_5, ymax = rich_wd_97_5),
                        col = "#ff8701", fatten = 1, position = "jitter")+
        stat_smooth(data = data, aes(x = imp1km, y = rich_wd_mean), fill = "#f2ff0d", col = "#ff8701",
                    method = "lm", formula = y ~ x + I(x^2), se = FALSE, size = 0.5) +
        geom_pointrange(data = data, aes(x = imp1km, y = rich_ww_mean, ymin = rich_ww_2_5, ymax = rich_ww_97_5),
                        col = "#d10073", fatten = 1)+
        stat_smooth(data = data, aes(x = imp1km, y = rich_ww_mean), fill = "#ff23c1", col = "#d10073",
                    method = "lm", formula = y ~ x + I(x^2), se = FALSE, size = 0.5) +
        geom_pointrange(data = data, aes(x = imp1km, y = rich_cw_mean, ymin = rich_cw_2_5, ymax = rich_cw_97_5),
                        col = "#2f3ca3", fatten = 1)+
        stat_smooth(data = data, aes(x = imp1km, y = rich_cw_mean), fill = "#7291fb", col = "#2f3ca3",
                    method = "lm", formula = y ~ x + I(x^2), se = FALSE, size = 0.5) +
        #geom_smooth(data = preds.sr.urb1, aes(x = impervious, y = median, ymin = lower95, ymax = upper95,),
        #            fill = "#7291fb", color = "#2f3ca3", se = FALSE, linetype = "dashed")+   # predict
        
        theme_bw() + 
        coord_cartesian(xlim=c(0, max(data$imp1km)), ylim=c(0,16))+
        labs(x = "Urbanization \n(% Impervious Surface)", y = "Species Richness")  +
        theme(axis.text.x = element_text(face = "bold"), 
              axis.text.y = element_text(face = "bold"), 
              axis.ticks = element_blank(),
              axis.title.x = element_text(face = "bold"), 
              axis.title.y = element_text(face = "bold"),
        ) 
      plot.all
      
      
      plot.all <- ggplot(data = data.long, aes(x = imp1km, y = rich_mean, ymin = rich_2_5, ymax = rich_97_5, group = season, color = season)) +
        geom_pointrange(fatten = 1, position = position_dodge2(width = 3))+    #, position = "jitter"
        # geom_errorbar(fatten = 1, position = position_dodge(width = 3))+
        # geom_point(position = position_dodge2(width = 3))+
        stat_smooth(data = data.long, aes(x = imp1km, y = rich_mean),
                    method = "lm", formula = y ~ x + I(x^2), se = FALSE, size = 0.5) +
        scale_color_manual(values = c( "#ff8701", "#d10073", "#2f3ca3")) +
        scale_fill_manual(values = c( "#f2ff0d", "#ff23c1", "#7291fb")) +
        theme_bw() + 
        coord_cartesian(xlim=c(0, max(data$imp1km)), ylim=c(0,13))+
        labs(x = "Urbanization \n(% Impervious Surface)", y = "Species Richness")  +
        theme(axis.text.x = element_text(face = "bold"), 
              axis.text.y = element_text(face = "bold"), 
              axis.ticks = element_blank(),
              axis.title.x = element_text(face = "bold"), 
              axis.title.y = element_text(face = "bold"),
        ) 
      plot.all
      
      
      ggsave("./figures/richness_actual_vs_urb_allseasons.png",
             plot.all,
             width = 6,
             height = 5,
             units = "in",
             dpi = 300)
      
      
##### Plot Predicted Richness across Hypothetical Gradients  #####
  ##### Urbanization vs. Species Richness Predictions by Season #####
      # Predicted species richness
      palette_ggdemo(wildflower)
      
      plot.wd <- ggplot() +
        theme_bw(base_size = 11) + 
        #geom_errorbar(data = data, aes(x = imp1km, y = rich_wd_mean, ymin = rich_wd_2_5, ymax = rich_wd_97_5), 
        #              lwd = 0.3, alpha = 0.3) +
        #geom_point(data = data, aes(x = imp1km, y = rich_wd_mean, ymin = rich_wd_2_5, ymax = rich_wd_97_5))+
        geom_ribbon(data = preds.sr.urb1, aes(x = impervious, y = median, ymin = lower95, ymax = upper95), alpha = 0.5, 
                    fill = "#f2ff0d") +
        geom_smooth(data = preds.sr.urb1, aes(x = impervious, y = median, ymin = lower95, ymax = upper95,),
                    fill = "#f2ff0d", color = "#ff8701", se = FALSE, size = 0.5, linetype = "dashed")+
        #stat_smooth(data = data, aes(x = imp1km, y = rich_wd_mean), fill = "#f2ff0d", col = "#ff8701",
        #            method = "lm", formula = y ~ x + I(x^2), se = FALSE, size = 0.5) +
        coord_cartesian(xlim=c(0, max(data$imp1km)), ylim=c(0,13))+
        labs(x = "Urbanization \n(% Impervious Surface)", y = "Species Richness")  +
        theme(axis.text.x = element_text(face = "bold"), 
              axis.text.y = element_text(face = "bold"),
              axis.ticks = element_blank(),
              axis.title.x = element_text(face = "bold"), 
              axis.title.y = element_text(face = "bold"),
        ) 
      
      plot.wd
      
      # 32, 54, 16, fill/color
      plot.ww <- ggplot() +
        theme_bw() + 
        geom_ribbon(data = preds.sr.urb2, aes(x = impervious, y = median, ymin = lower95, ymax = upper95), alpha = 0.5,
                    fill = "#ff23c1") + #"#ff23c1"
        geom_smooth(data = preds.sr.urb2, aes(x = impervious, y = median, ymin = lower95, ymax = upper95),
                    fill = "#ff23c1", col = "#d10073", se = FALSE, size = 0.5, linetype = "dashed")+
        #stat_smooth(data = data, aes(x = imp1km, y = rich_wd_mean), fill = "#ff23c1", col = "#d10073",
        #            method = "lm", formula = y ~ x + I(x^2), se = FALSE, size = 0.5) +
        coord_cartesian(xlim=c(0, max(data$imp1km)), ylim=c(0,13))+
        labs(x = "Urbanization \n(% Impervious Surface)", y = "Species Richness")  +
        theme(axis.text.x = element_text(face = "bold"), 
              axis.text.y = element_text(face = "bold"), 
              axis.ticks = element_blank(),
              axis.title.x = element_text(face = "bold"), 
              axis.title.y = element_text(face = "bold"),
        ) 
      
      plot.ww
      
      plot.cw <- ggplot() +
        geom_ribbon(data = preds.sr.urb3, aes(x = impervious, y = median, ymin = lower95, ymax = upper95), alpha = 0.5,
                    fill = "#7291fb") +
        geom_smooth(data = preds.sr.urb3, aes(x = impervious, y = median, ymin = lower95, ymax = upper95),
                    fill = "#7291fb", color = "#2f3ca3", se = FALSE, size = 0.5, linetype = "dashed")+
        #stat_smooth(data = data, aes(x = imp1km, y = rich_wd_mean), fill = "#7291fb", col = "#2f3ca3",
        #            method = "lm", formula = y ~ x + I(x^2), se = FALSE, size = 0.5) +
        theme_bw() + 
        coord_cartesian(xlim=c(0, max(data$imp1km)), ylim=c(0,13))+
        labs(x = "Urbanization \n(% Impervious Surface)", y = "Species Richness")  +
        theme(axis.text.x = element_text(face = "bold"), 
              axis.text.y = element_text(face = "bold"), 
              axis.ticks = element_blank(),
              axis.title.x = element_text(face = "bold"), 
              axis.title.y = element_text(face = "bold"),
        ) 
      
      plot.cw
  
  ggsave("./figures/richness_pred_vs_urb_1warmdry.png",
         plot.wd,
         width = 3,
         height = 2.5,
         units = "in",
         dpi = 300)
  
  ggsave("./figures/richness_pred_vs_urb_2warmwet.png",
         plot.ww,
         width = 3,
         height = 2.5,
         units = "in",
         dpi = 300)
  
  ggsave("./figures/richness_pred_vs_urb_3coolwet.png",
         plot.cw,
         width = 3,
         height = 2.5,
         units = "in",
         dpi = 300)
  
  plot.all <- ggplot() +
    theme_bw(base_size = 11) + 
    geom_ribbon(data = preds.sr.urb1, aes(x = impervious, y = median, ymin = lower95, ymax = upper95), alpha = 0.3, 
                fill = "#f2ff0d") +
    geom_smooth(data = preds.sr.urb1, aes(x = impervious, y = median, ymin = lower95, ymax = upper95,),
                fill = "#f2ff0d", color = "#ff8701", se = FALSE, linetype = "solid")+
    
    geom_ribbon(data = preds.sr.urb2, aes(x = impervious, y = median, ymin = lower95, ymax = upper95), alpha = 0.3,
                fill = "#ff23c1") +
    geom_smooth(data = preds.sr.urb2, aes(x = impervious, y = median, ymin = lower95, ymax = upper95),
                fill = "#ff23c1", col = "#d10073", se = FALSE, linetype = "dotted")+
    
    geom_ribbon(data = preds.sr.urb3, aes(x = impervious, y = median, ymin = lower95, ymax = upper95), alpha = 0.3,
                fill = "#7291fb") +
    geom_smooth(data = preds.sr.urb3, aes(x = impervious, y = median, ymin = lower95, ymax = upper95),
                fill = "#7291fb", color = "#2f3ca3",  se = FALSE, linetype = "dashed")+
    coord_cartesian(xlim=c(0, max(data$imp1km)), ylim=c(0,13))+
    labs(x = "Urbanization \n(% Impervious Surface)", y = "Species Richness")  +
    theme(axis.text.x = element_text(face = "bold"), 
          axis.text.y = element_text(face = "bold"), 
          axis.title.x = element_text(face = "bold", size = 12), 
          axis.title.y = element_text(face = "bold", size = 12),
    ) 
  plot.all
  

##### Plot Actual AND Predicted Richness across Hypothetical Gradients  #####
  ##### Urbanization vs. Species Richness Predictions by Season #####
  # Predicted species richness
  # palette_ggdemo(wildflower)
  
  plot.wd <- ggplot() +
    theme_bw(base_size = 11) + 
    #geom_errorbar(data = data, aes(x = imp1km, y = rich_wd_mean, ymin = rich_wd_2_5, ymax = rich_wd_97_5), 
    #              lwd = 0.3, alpha = 0.3) +
    geom_point(data = data, aes(x = imp1km, y = rich_wd_mean, ymin = rich_wd_2_5, ymax = rich_wd_97_5), color = "#ff8701")+
    geom_ribbon(data = preds.sr.urb1, aes(x = impervious, y = median, ymin = lower95, ymax = upper95), alpha = 0.5, 
                fill = "#f2ff0d") +
    geom_smooth(data = preds.sr.urb1, aes(x = impervious, y = median, ymin = lower95, ymax = upper95,),
                fill = "#f2ff0d", color = "#ff8701", se = FALSE, size = 0.5, linetype = "dashed")+
    #stat_smooth(data = data, aes(x = imp1km, y = rich_wd_mean), fill = "#f2ff0d", col = "#ff8701",
    #            method = "lm", formula = y ~ x + I(x^2), se = FALSE, size = 0.5) +
    coord_cartesian(xlim=c(0, max(data$imp1km)), ylim=c(0,13))+
    labs(x = "Urbanization \n(% Impervious Surface)", y = "Species Richness")  +
    theme(axis.text.x = element_text(face = "bold"), 
          axis.text.y = element_text(face = "bold"),
          axis.ticks = element_blank(),
          axis.title.x = element_text(face = "bold"), 
          axis.title.y = element_text(face = "bold"),
    ) 
  
  plot.wd
  
  # 32, 54, 16, fill/color
  plot.ww <- ggplot() +
    theme_bw() + 
    geom_ribbon(data = preds.sr.urb2, aes(x = impervious, y = median, ymin = lower95, ymax = upper95), alpha = 0.5,
                fill = "#ff23c1") + #"#ff23c1"
    geom_smooth(data = preds.sr.urb2, aes(x = impervious, y = median, ymin = lower95, ymax = upper95),
                fill = "#ff23c1", col = "#d10073", se = FALSE, size = 0.5, linetype = "dashed")+
    #stat_smooth(data = data, aes(x = imp1km, y = rich_wd_mean), fill = "#ff23c1", col = "#d10073",
    #            method = "lm", formula = y ~ x + I(x^2), se = FALSE, size = 0.5) +
    coord_cartesian(xlim=c(0, max(data$imp1km)), ylim=c(0,12))+
    labs(x = "Urbanization \n(% Impervious Surface)", y = "Species Richness")  +
    theme(axis.text.x = element_text(face = "bold"), 
          axis.text.y = element_text(face = "bold"), 
          axis.ticks = element_blank(),
          axis.title.x = element_text(face = "bold"), 
          axis.title.y = element_text(face = "bold"),
    ) 
  
  plot.ww
  
  plot.cw <- ggplot() +
    geom_ribbon(data = preds.sr.urb3, aes(x = impervious, y = median, ymin = lower95, ymax = upper95), alpha = 0.5,
                fill = "#7291fb") +
    geom_smooth(data = preds.sr.urb3, aes(x = impervious, y = median, ymin = lower95, ymax = upper95),
                fill = "#7291fb", color = "#2f3ca3", se = FALSE, size = 0.5, linetype = "dashed")+
    #stat_smooth(data = data, aes(x = imp1km, y = rich_wd_mean), fill = "#7291fb", col = "#2f3ca3",
    #            method = "lm", formula = y ~ x + I(x^2), se = FALSE, size = 0.5) +
    theme_bw() + 
    coord_cartesian(xlim=c(0, max(data$imp1km)), ylim=c(0,12))+
    labs(x = "Urbanization \n(% Impervious Surface)", y = "Species Richness")  +
    theme(axis.text.x = element_text(face = "bold"), 
          axis.text.y = element_text(face = "bold"), 
          axis.ticks = element_blank(),
          axis.title.x = element_text(face = "bold"), 
          axis.title.y = element_text(face = "bold"),
    ) 
  
  plot.cw
  
  ggsave("./figures/richness_pred_vs_urb_1warmdry.png",
         plot.wd,
         width = 3,
         height = 2.5,
         units = "in",
         dpi = 300)
  
  ggsave("./figures/richness_pred_vs_urb_2warmwet.png",
         plot.ww,
         width = 3,
         height = 2.5,
         units = "in",
         dpi = 300)
  
  ggsave("./figures/richness_pred_vs_urb_3coolwet.png",
         plot.cw,
         width = 3,
         height = 2.5,
         units = "in",
         dpi = 300)
 

  
  plot.all <- ggplot(data = preds.sr.urb, aes(group = season, color = season, 
                                              linetype = season, fill = season)) +
    theme_bw(base_size = 11) + 
    geom_ribbon(aes(x = impervious, y = median, ymin = lower95, ymax = upper95), alpha = 0.3, color = "transparent") +
    geom_smooth(aes(x = impervious, y = median, ymin = lower95, ymax = upper95,),
                se = FALSE, method = "gam") +
    scale_color_manual(values = c( "#ff8701", "#d10073", "#2f3ca3")) +
    scale_fill_manual(values = c( "#f2ff0d", "#ff23c1", "#7291fb"))+
    geom_point(data = data.long, aes(x = imp1km, y = rich_mean)) +
    coord_cartesian(xlim=c(0, max(data$imp1km)), ylim=c(0,13))+
    labs(x = "Urbanization \n(% Impervious Surface)", y = "Species Richness")  +
    theme(axis.text.x = element_text(face = "bold"), 
          axis.text.y = element_text(face = "bold"), 
          legend.position = "none",
          axis.title.x = element_text(face = "bold", size = 12), 
          axis.title.y = element_text(face = "bold", size = 12),
    ) 
  plot.all
  
  
  ggsave("./figures/richness__vs_urb_allseasons_actualandpredicted.png",
         plot.all,
         width = 4,
         height = 4,
         units = "in",
         dpi = 300)
  
  