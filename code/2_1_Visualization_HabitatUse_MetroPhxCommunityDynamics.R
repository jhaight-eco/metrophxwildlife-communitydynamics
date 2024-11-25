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
  
  # visualization packages
  library(ggplot2)
  library(viridis)
  library(RColorBrewer)
  library(png)
  library(scales)


##### Import Data #####
  # all the data input into the models
  load("./data/model_inputs/ModelInputData_DCM.RData")
  # The main dataset containing all the site and regional covariates 
  # and summarized community composition parameters
  data <- read.csv("./data/model_outputs/data_PhoenixCommDynamics.csv")

  
  # Community-average dynamic use parameters (occupancy, persistence, and colonization)
  # predicted across each environmental covariate with all others held constant
  data.pred.comm <- read.csv("./data/model_outputs/predictions_commparams.csv")
  
  # Species-level parameter predictions based on the actual trait values of the 20 species
  pred.spp.urb <- read.csv("./data/model_outputs/predictions_spp_urban.csv")    
  pred.spp.pdiv <- read.csv("./data/model_outputs/predictions_spp_patchdiversity.csv")    
  pred.spp.veg <- read.csv("./data/model_outputs/predictions_spp_vegetation.csv")
  
  # Species-level parameter predictions based on three idealized trait levels (with the other trait held constant)
  data.pred.bm <- read.csv("./data/model_outputs/predictions_ideal_3spp_size.csv")
  data.pred.di <- read.csv("./data/model_outputs/predictions_ideal_3spp_diet.csv")
  
  
  # To plot the trait variation without the 'moderate' trait species, remove those from the predictions dataset
  data.pred.bm <- data.pred.bm %>% filter(species != "sp2")
  data.pred.di <- data.pred.di %>% filter(species != "sp2")
  

  
#### Plot Community-average Predictions ####
  ##### Across a gradient of urbanization #####
  plot.psi.urb <- ggplot(data = data.pred.comm, aes(x = impervious, y = psi_urb_med)) +
    theme_classic(base_size = 20)+
    geom_ribbon(aes(ymin = psi_urb_low95, ymax = psi_urb_high95), fill = "gray40", alpha = 0.3) +
    geom_ribbon(aes(ymin = psi_urb_low85, ymax = psi_urb_high85), fill = "gray40", alpha = 0.5) +
    geom_smooth(se = FALSE, color = "black")+
    #geom_point(color = "gray4")+
    scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75)) +
    coord_cartesian(ylim = c(0,1)) +
    labs(x = "Urbanization \n(% Impervious Surface)", y = "Occupancy") +
    theme(axis.text.x = element_text(face = "bold"), 
          axis.text.y = element_text(face = "bold"), 
          axis.ticks = element_blank(),
          axis.title.x = element_text(face = "bold", size = 18), 
          axis.title.y = element_text(face = "bold", size = 18)
    ) 
  plot.psi.urb  # substantial relationship
  
  plot.phi.urb <- ggplot(data = data.pred.comm, aes(x = impervious, y = phi_urb_med)) +
    theme_classic(base_size = 20)+
    geom_ribbon(aes(ymin = phi_urb_low95, ymax = phi_urb_high95), fill = "gray40", alpha = 0.3) +
    geom_ribbon(aes(ymin = phi_urb_low85, ymax = phi_urb_high85), fill = "gray40", alpha = 0.5) +
    geom_smooth(se = FALSE, color = "black")+
    #geom_point(color = "gray4")+
    scale_y_continuous(labels = label_number(accuracy = 0.01)) +
    coord_cartesian(ylim = c(0,1)) +
    labs(x = "Urbanization \n(% Impervious Surface)", y = "Persistence") +
    theme(axis.text.x = element_text(face = "bold"), 
          axis.text.y = element_text(face = "bold"), 
          axis.ticks = element_blank(),
          axis.title.x = element_text(face = "bold", size = 18), 
          axis.title.y = element_text(face = "bold", size = 18)
    ) 
  plot.phi.urb
  
  plot.gamma.urb <- ggplot(data = data.pred.comm, aes(x = impervious, y = gamma_urb)) +
    theme_classic(base_size = 20)+
    geom_ribbon(aes(ymin = gamma_urb_low95, ymax = gamma_urb_high95), fill = "gray40", alpha = 0.3) +
    geom_ribbon(aes(ymin = gamma_urb_low85, ymax = gamma_urb_high85), fill = "gray40", alpha = 0.5) +
    geom_smooth(se = FALSE, color = "black")+
    #geom_point(color = "gray4")+
    scale_y_continuous(labels = label_number(accuracy = 0.01)) +
    coord_cartesian(ylim = c(0,1)) +
    labs(x = "Urbanization \n(% Impervious Surface)", y = "Colonization") +
    theme(axis.text.x = element_text(face = "bold"), 
          axis.text.y = element_text(face = "bold"), 
          axis.ticks = element_blank(),
          axis.title.x = element_text(face = "bold", size = 18), 
          axis.title.y = element_text(face = "bold", size = 18)
    ) 
  plot.gamma.urb
  
  ggsave("./figures/pred_1psi_urb_comm.png",
         plot.psi.urb,
         width = 4,
         height = 4,
         units = "in",
         dpi = 300)
  
  ggsave("./figures/pred_2phi_urb_comm.png",
         plot.phi.urb,
         width = 4,
         height = 4,
         units = "in",
         dpi = 300)
  
  ggsave("./figures/pred_3gamma_urb_comm.png",
         plot.gamma.urb,
         width = 4,
         height = 4,
         units = "in",
         dpi = 300)
  
  
  ##### Across a gradient of patch diversity #####
  plot.psi.pdiv <- ggplot(data = data.pred.comm, aes(x = het, y = psi_pdiv)) +
    theme_classic(base_size = 20)+
    geom_ribbon(aes(ymin = psi_pdiv_low95, ymax = psi_pdiv_high95), fill =  "gray40", alpha = 0.3) +
    geom_ribbon(aes(ymin = psi_pdiv_low85, ymax = psi_pdiv_high85), fill =  "gray40", alpha = 0.5) +
    geom_smooth(se = FALSE, color =  "black")+
    #geom_point(color = "gray4")+
    scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75)) +
    coord_cartesian(ylim = c(0,1)) +
    labs(x = "Patch Diversity \n(Shannon Diversity Index)", y = "Occupancy") +
    theme(axis.text.x = element_text(face = "bold"), 
          axis.text.y = element_text(face = "bold"), 
          axis.ticks = element_blank(),
          axis.title.x = element_text(face = "bold", size = 18), 
          axis.title.y = element_text(face = "bold", size = 18)
    ) 
  plot.psi.pdiv
  
  plot.phi.pdiv <- ggplot(data = data.pred.comm, aes(x = het, y = phi_pdiv)) +
    theme_classic(base_size = 20)+
    geom_ribbon(aes(ymin = phi_pdiv_low95, ymax = phi_pdiv_high95), fill =  "gray40", alpha = 0.3) +
    geom_ribbon(aes(ymin = phi_pdiv_low85, ymax = phi_pdiv_high85), fill =  "gray40", alpha = 0.5) +
    geom_smooth(se = FALSE, color =  "black")+
    #geom_point(color = "gray4")+
    scale_y_continuous(labels = label_number(accuracy = 0.01)) +
    coord_cartesian(ylim = c(0,1)) +
    labs(x = "Patch Diversity \n(Shannon Diversity Index)", y = "Persistence") +
    theme(axis.text.x = element_text(face = "bold"), 
          axis.text.y = element_text(face = "bold"), 
          axis.ticks = element_blank(),
          axis.title.x = element_text(face = "bold", size = 18), 
          axis.title.y = element_text(face = "bold", size = 18)
    ) 
  plot.phi.pdiv
  
  plot.gamma.pdiv <- ggplot(data = data.pred.comm, aes(x = het, y = gamma_pdiv)) +
    theme_classic(base_size = 20)+
    geom_ribbon(aes(ymin = gamma_pdiv_low95, ymax = gamma_pdiv_high95), fill =  "gray40", alpha = 0.3) +
    geom_ribbon(aes(ymin = gamma_pdiv_low85, ymax = gamma_pdiv_high85), fill =  "gray40", alpha = 0.5) +
    geom_smooth(se = FALSE, color =  "black")+
    #geom_point(color = "gray4")+
    scale_y_continuous(labels = label_number(accuracy = 0.01)) +
    coord_cartesian(ylim = c(0,1)) +
    labs(x = "Patch Diversity \n(Shannon Diversity Index)", y = "Colonization") +
    theme(axis.text.x = element_text(face = "bold"), 
          axis.text.y = element_text(face = "bold"), 
          axis.ticks = element_blank(),
          axis.title.x = element_text(face = "bold", size = 18), 
          axis.title.y = element_text(face = "bold", size = 18)
    ) 
  plot.gamma.pdiv
  
  ggsave("./figures/pred_1psi_pdiv_comm.png",
         plot.psi.pdiv,
         width = 4,
         height = 4,
         units = "in",
         dpi = 300)
  
  ggsave("./figures/pred_2phi_pdiv_comm.png",
         plot.phi.pdiv,
         width = 4,
         height = 4,
         units = "in",
         dpi = 300)
  
  ggsave("./figures/pred_3gamma_pdiv_comm.png",
         plot.gamma.pdiv,
         width = 4,
         height = 4,
         units = "in",
         dpi = 300)
  
  ##### Across a gradient of vegetation greenness #####
  plot.psi.veg <- ggplot(data = data.pred.comm, aes(x = veg, y = psi_veg)) +
    theme_classic(base_size = 20)+
    geom_ribbon(aes(ymin = psi_veg_low95, ymax = psi_veg_high95), fill =  "gray40", alpha = 0.3) +
    geom_ribbon(aes(ymin = psi_veg_low85, ymax = psi_veg_high85), fill =  "gray40", alpha = 0.5) +
    geom_smooth(se = FALSE, color =  "black")+
    #geom_point(color =  "gray4")+
    scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75)) +
    coord_cartesian(ylim = c(0,1)) +
    labs(x = "Vegetation Greenness \n(NDVI)", y = "Occupancy") +
    theme(axis.text.x = element_text(face = "bold"), 
          axis.text.y = element_text(face = "bold"), 
          axis.ticks = element_blank(),
          axis.title.x = element_text(face = "bold", size = 18), 
          axis.title.y = element_text(face = "bold", size = 18)
    ) 
  plot.psi.veg 
  
  plot.phi.veg <- ggplot(data = data.pred.comm, aes(x = veg, y = phi_veg)) +
    theme_classic(base_size = 20)+
    geom_ribbon(aes(ymin = phi_veg_low95, ymax = phi_veg_high95), fill =  "gray40", alpha = 0.3) +
    geom_ribbon(aes(ymin = phi_veg_low85, ymax = phi_veg_high85), fill =  "gray40", alpha = 0.5) +
    geom_smooth(se = FALSE, color =  "black")+
    #geom_point(color =  "gray4")+
    scale_y_continuous(labels = label_number(accuracy = 0.01)) +
    coord_cartesian(ylim = c(0,1)) +
    labs(x = "Vegetation Greenness \n(NDVI)", y = "Persistence") +
    theme(axis.text.x = element_text(face = "bold"), 
          axis.text.y = element_text(face = "bold"), 
          axis.ticks = element_blank(),
          axis.title.x = element_text(face = "bold", size = 18), 
          axis.title.y = element_text(face = "bold", size = 18)
    ) 
  plot.phi.veg 
  
  plot.gamma.veg <- ggplot(data = data.pred.comm, aes(x = veg, y = gamma_veg)) +
    theme_classic(base_size = 20)+
    geom_ribbon(aes(ymin = gamma_veg_low95, ymax = gamma_veg_high95), fill =  "gray40", alpha = 0.3) +
    geom_ribbon(aes(ymin = gamma_veg_low85, ymax = gamma_veg_high85), fill =  "gray40", alpha = 0.5) +
    geom_smooth(se = FALSE, color =  "black")+
    #geom_point(color =  "gray4")+
    scale_y_continuous(labels = label_number(accuracy = 0.01)) +
    coord_cartesian(ylim = c(0,1)) +
    labs(x = "Vegetation Greenness \n(NDVI)", y = "Colonization") +
    theme(axis.text.x = element_text(face = "bold"), 
          axis.text.y = element_text(face = "bold"), 
          axis.ticks = element_blank(),
          axis.title.x = element_text(face = "bold", size = 18), 
          axis.title.y = element_text(face = "bold", size = 18)
    ) 
  plot.gamma.veg 
   
  
  ggsave("./figures/pred_1psi_veg_comm.png",
         plot.psi.veg,
         width = 4,
         height = 4,
         units = "in",
         dpi = 300)
  
  ggsave("./figures/pred_2phi_veg_comm.png",
         plot.phi.veg,
         width = 4,
         height = 4,
         units = "in",
         dpi = 300)
  
  ggsave("./figures/pred_3gamma_veg_comm.png",
         plot.gamma.veg,
         width = 4,
         height = 4,
         units = "in",
         dpi = 300)

#### Plot Species-Level Predictions for Species with Contrasting Body Size ####   
  ##### Size interaction with urbanization relationship #####
  # Use/Occupancy
  plot.psi.urb <- ggplot(data = data.pred.bm, aes(x = impervious, y = psi_urb_med, 
                                                  group = logmass, fill = logmass, col = logmass)) +
    geom_smooth(method = "gam", lwd = 2, col = "black", se = FALSE) +
    geom_smooth(method = "gam", lwd = 1, se = FALSE) +
    theme_classic(base_size = 20) +
    geom_ribbon(aes(ymin = psi_urb_low95, ymax = psi_urb_high95), alpha = 0.4, color = NA) +
    scale_color_viridis(option = "plasma", direction = -1, begin = 0.25, end = 0.8,
                        trans = "log10", breaks = log10(c(0.5, 1, 5, 20)*1000), labels = c("0.5", "1", "5", "20")) +
    scale_fill_viridis(option = "plasma", direction = -1, begin = 0.25, end = 0.8,
                       trans = "log10", breaks = log10(c(0.5, 1, 5, 20)*1000), labels = c("0.5", "1", "5", "20")) +
    coord_cartesian(ylim = c(0,1)) +
    labs(x = "Urbanization \n(% Impervious Surface)", y = "Occupancy", color = "Body \nMass \n(kg)") +
    #gghighlight(species %in% c( "cottontail_desert", "mountain_lion", "bobcat", "jackrabbit_bt")) + 
    theme(axis.text.x = element_text(face = "bold"), 
          axis.text.y = element_text(face = "bold"), 
          legend.position = "none",
          legend.title = element_text(face = "bold"),
          legend.text = element_text(face = "bold"),
          axis.ticks = element_blank(),
          axis.title.x = element_text(face = "bold", size = 18), 
          axis.title.y = element_text(face = "bold", size = 18))
  plot.psi.urb
  
  ggsave("./figures/traitvarideal_urbsize_vs_1use.png",
         plot.psi.urb,
         width = 4,
         height = 4,
         units = "in",
         dpi = 300)
  
  # Persistence
  plot.phi.urb <- ggplot(data = data.pred.bm, aes(x = impervious, y = phi_urb_med, 
                                                  group = logmass, fill = logmass, col = logmass)) +
    geom_smooth(method = "gam", lwd = 2, col = "black", se = FALSE) +
    geom_smooth(method = "gam", lwd = 1, se = FALSE) +
    theme_classic(base_size = 20) +
    geom_ribbon(aes(ymin = phi_urb_low95, ymax = phi_urb_high95), alpha = 0.4, color = NA) +
    scale_color_viridis(option = "plasma", direction = -1, begin = 0.25, end = 0.8,
                        trans = "log10", breaks = log10(c(0.5, 1, 5, 20)*1000), labels = c("0.5", "1", "5", "20")) +
    scale_fill_viridis(option = "plasma", direction = -1, begin = 0.25, end = 0.8,
                       trans = "log10", breaks = log10(c(0.5, 1, 5, 20)*1000), labels = c("0.5", "1", "5", "20")) +
    coord_cartesian(ylim = c(0,1)) +
    labs(x = "Urbanization \n(% Impervious Surface)", y = "Persistence", color = "Body \nMass \n(kg)") +
    theme(axis.text.x = element_text(face = "bold"), 
          axis.text.y = element_text(face = "bold"), 
          legend.position = "none",
          legend.title = element_text(face = "bold"),
          legend.text = element_text(face = "bold"),
          axis.ticks = element_blank(),
          axis.title.x = element_text(face = "bold", size = 18), 
          axis.title.y = element_text(face = "bold", size = 18))
  plot.phi.urb
  
  ggsave("./figures/traitvarideal_urbsize_vs_2persistence.png",
         plot.phi.urb,
         width = 4,
         height = 4,
         units = "in",
         dpi = 300)
  
  # Colonization
  plot.gamma.urb <- ggplot(data = data.pred.bm, aes(x = impervious, y = gamma_urb_med, 
                                                    group = logmass, fill = logmass, col = logmass)) +
    geom_smooth(method = "gam", lwd = 2, col = "black", se = FALSE) +
    geom_smooth(method = "gam", lwd = 1, se = FALSE) +
    theme_classic(base_size = 20) +
    geom_ribbon(aes(ymin = gamma_urb_low95, ymax = gamma_urb_high95), alpha = 0.4, color = NA) +
    scale_color_viridis(option = "plasma", direction = -1, begin = 0.25, end = 0.8,
                        trans = "log10", breaks = log10(c(0.5, 1, 5, 20)*1000), labels = c("0.5", "1", "5", "20")) +
    scale_fill_viridis(option = "plasma", direction = -1, begin = 0.25, end = 0.8,
                       trans = "log10", breaks = log10(c(0.5, 1, 5, 20)*1000), labels = c("0.5", "1", "5", "20")) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Urbanization \n(% Impervious Surface)", y = "Colonization", color = "Body \nMass \n(kg)") +
    theme(axis.text.x = element_text(face = "bold"), 
          axis.text.y = element_text(face = "bold"), 
          # legend.position = "none",
          legend.title = element_text(face = "bold"),
          legend.text = element_text(face = "bold"),
          axis.ticks = element_blank(),
          axis.title.x = element_text(face = "bold", size = 18), 
          axis.title.y = element_text(face = "bold", size = 18))
  plot.gamma.urb
  
  ggsave("./figures/traitvarideal_urbsize_vs_3colonization.png",
         plot.gamma.urb,
         width = 4,
         height = 4,
         units = "in",
         dpi = 300)
  
  ##### Size interaction with patch diversity relationship #####
  
  setwd("~/GitHub/metrophxwildlife-communitydynamics")     
  
  # Use/Occupancy
  plot.psi.pdiv <- ggplot(data = data.pred.bm, aes(x = sdhi, y = psi_pdiv_med, 
                                                   group = logmass, fill = logmass, col = logmass)) +
    geom_smooth(method = "gam", lwd = 2, col = "black", se = FALSE) +
    geom_smooth(method = "gam", lwd = 1, se = FALSE) +
    theme_classic(base_size = 20) +
    geom_ribbon(aes(ymin = psi_pdiv_low95, ymax = psi_pdiv_high95), alpha = 0.4, color = NA) +
    scale_color_viridis(option = "plasma", direction = -1, begin = 0.25, end = 0.8,
                        trans = "log10", breaks = log10(c(0.5, 1, 5, 20)*1000), labels = c("0.5", "1", "5", "20")) +
    scale_fill_viridis(option = "plasma", direction = -1, begin = 0.25, end = 0.8,
                       trans = "log10", breaks = log10(c(0.5, 1, 5, 20)*1000), labels = c("0.5", "1", "5", "20")) +
    coord_cartesian(ylim = c(0,1)) +
    labs(x = "Patch Diversity \n(Shannon Diversity Index)", y = "Occupancy", color = "Body \nMass \n(kg)") +
    #gghighlight(species %in% c( "cottontail_desert", "mountain_lion", "bobcat", "jackrabbit_bt")) + 
    theme(axis.text.x = element_text(face = "bold"), 
          axis.text.y = element_text(face = "bold"), 
          legend.position = "none",
          legend.title = element_text(face = "bold"),
          legend.text = element_text(face = "bold"),
          axis.ticks = element_blank(),
          axis.title.x = element_text(face = "bold", size = 18), 
          axis.title.y = element_text(face = "bold", size = 18))
  plot.psi.pdiv
  
  ggsave("./figures/traitvarideal_pdivsize_vs_1use.png",
         plot.psi.pdiv,
         width = 4,
         height = 4,
         units = "in",
         dpi = 300)
  
  # Persistence
  plot.phi.pdiv <- ggplot(data = data.pred.bm, aes(x = sdhi, y = phi_pdiv_med, 
                                                   group = logmass, fill = logmass, col = logmass)) +
    geom_smooth(method = "gam", lwd = 2, col = "black", se = FALSE) +
    geom_smooth(method = "gam", lwd = 1, se = FALSE) +
    theme_classic(base_size = 20) +
    geom_ribbon(aes(ymin = phi_pdiv_low95, ymax = phi_pdiv_high95), alpha = 0.4, color = NA) +
    scale_color_viridis(option = "plasma", direction = -1, begin = 0.25, end = 0.8,
                        trans = "log10", breaks = log10(c(0.5, 1, 5, 20)*1000), labels = c("0.5", "1", "5", "20")) +
    scale_fill_viridis(option = "plasma", direction = -1, begin = 0.25, end = 0.8,
                       trans = "log10", breaks = log10(c(0.5, 1, 5, 20)*1000), labels = c("0.5", "1", "5", "20")) +
    coord_cartesian(ylim = c(0,1)) +
    labs(x = "Patch Diversity \n(Shannon Diversity Index)", y = "Persistence", color = "Body \nMass \n(kg)") +
    theme(axis.text.x = element_text(face = "bold"), 
          axis.text.y = element_text(face = "bold"), 
          legend.position = "none",
          legend.title = element_text(face = "bold"),
          legend.text = element_text(face = "bold"),
          axis.ticks = element_blank(),
          axis.title.x = element_text(face = "bold", size = 18), 
          axis.title.y = element_text(face = "bold", size = 18))
  plot.phi.pdiv     # substantial relationship in model without quadratic urbanization term
  
  ggsave("./figures/traitvarideal_pdivsize_vs_2persistence.png",
         plot.phi.pdiv,
         width = 4,
         height = 4,
         units = "in",
         dpi = 300)
  
  # Colonization
  plot.gamma.pdiv <- ggplot(data = data.pred.bm, aes(x = sdhi, y = gamma_pdiv_med, 
                                                     group = logmass, fill = logmass, col = logmass)) +
    geom_smooth(method = "gam", lwd = 2, col = "black", se = FALSE) +
    geom_smooth(method = "gam", lwd = 1, se = FALSE) +
    theme_classic(base_size = 20) +
    geom_ribbon(aes(ymin = gamma_pdiv_low95, ymax = gamma_pdiv_high95), alpha = 0.4, color = NA) +
    scale_color_viridis(option = "plasma", direction = -1, begin = 0.25, end = 0.8,
                        trans = "log10", breaks = log10(c(0.5, 1, 5, 20)*1000), labels = c("0.5", "1", "5", "20")) +
    scale_fill_viridis(option = "plasma", direction = -1, begin = 0.25, end = 0.8,
                       trans = "log10", breaks = log10(c(0.5, 1, 5, 20)*1000), labels = c("0.5", "1", "5", "20")) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Patch Diversity \n(Shannon Diversity Index)", y = "Colonization", color = "Body \nMass \n(kg)") +
    theme(axis.text.x = element_text(face = "bold"), 
          axis.text.y = element_text(face = "bold"), 
          legend.position = "none",
          legend.title = element_text(face = "bold"),
          legend.text = element_text(face = "bold"),
          axis.ticks = element_blank(),
          axis.title.x = element_text(face = "bold", size = 18), 
          axis.title.y = element_text(face = "bold", size = 18))
  plot.gamma.pdiv
  
  ggsave("./figures/traitvarideal_hetsize_vs_3colonization.png",
         plot.gamma.pdiv,
         width = 4,
         height = 4,
         units = "in",
         dpi = 300)
  
  
  ##### Size interaction with vegetation greenness relationship #####
  # Use/Occupancy
  plot.psi.veg <- ggplot(data = data.pred.bm, aes(x = veg, y = psi_veg_med, 
                                                  group = logmass, fill = logmass, col = logmass)) +
    geom_smooth(method = "gam", lwd = 2, col = "black", se = FALSE) +
    geom_smooth(method = "gam", lwd = 1, se = FALSE) +
    theme_classic(base_size = 20) +
    geom_ribbon(aes(ymin = psi_veg_low95, ymax = psi_veg_high95), alpha = 0.4, color = NA) +
    scale_color_viridis(option = "plasma", direction = -1, begin = 0.25, end = 0.8,
                        trans = "log10", breaks = log10(c(0.5, 1, 5, 20)*1000), labels = c("0.5", "1", "5", "20")) +
    scale_fill_viridis(option = "plasma", direction = -1, begin = 0.25, end = 0.8,
                       trans = "log10", breaks = log10(c(0.5, 1, 5, 20)*1000), labels = c("0.5", "1", "5", "20")) +
    coord_cartesian(ylim = c(0,1)) +
    labs(x = "Vegetation Greenness \n(NDVI)", y = "Occupancy", color = "Body \nMass \n(kg)") +
    #gghighlight(species %in% c( "cottontail_desert", "mountain_lion", "bobcat", "jackrabbit_bt")) + 
    theme(axis.text.x = element_text(face = "bold"), 
          axis.text.y = element_text(face = "bold"), 
          legend.position = "none",
          legend.title = element_text(face = "bold"),
          legend.text = element_text(face = "bold"),
          axis.ticks = element_blank(),
          axis.title.x = element_text(face = "bold", size = 18), 
          axis.title.y = element_text(face = "bold", size = 18))
  plot.psi.veg
  
  ggsave("./figures/traitvarideal_vegsize_vs_1use.png",
         plot.psi.veg,
         width = 4,
         height = 4,
         units = "in",
         dpi = 300)
  
  # Persistence
  plot.phi.veg <- ggplot(data = data.pred.bm, aes(x = veg, y = phi_veg_med, 
                                                  group = logmass, fill = logmass, col = logmass)) +
    geom_smooth(method = "gam", lwd = 2, col = "black", se = FALSE) +
    geom_smooth(method = "gam", lwd = 1, se = FALSE) +
    theme_classic(base_size = 20) +
    geom_ribbon(aes(ymin = phi_veg_low95, ymax = phi_veg_high95), alpha = 0.4, color = NA) +
    scale_color_viridis(option = "plasma", direction = -1, begin = 0.25, end = 0.8,
                        trans = "log10", breaks = log10(c(0.5, 1, 5, 20)*1000), labels = c("0.5", "1", "5", "20")) +
    scale_fill_viridis(option = "plasma", direction = -1, begin = 0.25, end = 0.8,
                       trans = "log10", breaks = log10(c(0.5, 1, 5, 20)*1000), labels = c("0.5", "1", "5", "20")) +
    coord_cartesian(ylim = c(0,1)) +
    labs(x = "Vegetation Greenness \n(NDVI)", y = "Persistence", color = "Body \nMass \n(kg)") +
    theme(axis.text.x = element_text(face = "bold"), 
          axis.text.y = element_text(face = "bold"), 
          legend.position = "none",
          legend.title = element_text(face = "bold"),
          legend.text = element_text(face = "bold"),
          axis.ticks = element_blank(),
          axis.title.x = element_text(face = "bold", size = 18), 
          axis.title.y = element_text(face = "bold", size = 18))
  plot.phi.veg
  
  ggsave("./figures/traitvarideal_vegsize_vs_2persistence.png",
         plot.phi.veg,
         width = 4,
         height = 4,
         units = "in",
         dpi = 300)
  
  # Colonization
  plot.gamma.veg <- ggplot(data = data.pred.bm, aes(x = veg, y = gamma_veg_med, 
                                                    group = logmass, fill = logmass, col = logmass)) +
    geom_smooth(method = "gam", lwd = 2, col = "black", se = FALSE) +
    geom_smooth(method = "gam", lwd = 1, se = FALSE) +
    theme_classic(base_size = 20) +
    geom_ribbon(aes(ymin = gamma_veg_low95, ymax = gamma_veg_high95), alpha = 0.4, color = NA) +
    scale_color_viridis(option = "plasma", direction = -1, begin = 0.25, end = 0.8,
                        trans = "log10", breaks = log10(c(0.5, 1, 5, 20)*1000), labels = c("0.5", "1", "5", "20")) +
    scale_fill_viridis(option = "plasma", direction = -1, begin = 0.25, end = 0.8,
                       trans = "log10", breaks = log10(c(0.5, 1, 5, 20)*1000), labels = c("0.5", "1", "5", "20")) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Vegetation Greenness \n(NDVI)", y = "Colonization", color = "Body \nMass \n(kg)") +
    theme(axis.text.x = element_text(face = "bold"), 
          axis.text.y = element_text(face = "bold"), 
          legend.position = "none",
          legend.title = element_text(face = "bold"),
          legend.text = element_text(face = "bold"),
          axis.ticks = element_blank(),
          axis.title.x = element_text(face = "bold", size = 18), 
          axis.title.y = element_text(face = "bold", size = 18))
  plot.gamma.veg
  
  ggsave("./figures/traitvarideal_vegsize_vs_3colonization.png",
         plot.gamma.veg,
         width = 4,
         height = 4,
         units = "in",
         dpi = 300)
  
    
  
#### Plot Species-Level Predictions for Species with Contrasting Diet Diversity ####  
  ##### Diet interaction with urbanization relationships #####
  # Use/Occupancy
  plot.psi.urb <- ggplot(data = data.pred.di, aes(x = impervious, y = psi_urb_med, 
                                                  group = dietdiv, fill = dietdiv, col = dietdiv)) +
    geom_smooth(method = "gam", lwd = 2, col = "black", se = FALSE) +
    geom_smooth(method = "gam", lwd = 1, se = FALSE) +
    theme_classic(base_size = 20) +
    geom_ribbon(aes(ymin = psi_urb_low95, ymax = psi_urb_high95), alpha = 0.4, color = NA) +
    #scale_color_distiller(palette = "PuBuGn", direction = -1)+
    #scale_fill_distiller(palette = "PuBuGn", direction = -1)+
    # scale_color_manual(values = pal_viridis(direction = -1, option = "viridis"))
    scale_color_viridis(option = "viridis", direction = -1, end = 0.5, begin = 0.95) +
    scale_fill_viridis(option = "viridis", direction = -1, end = 0.5, begin = 0.95) +
    scale_y_continuous(labels = label_number(accuracy = 0.01)) +
    coord_cartesian(ylim = c(0,1)) +
    labs(x = "Urbanization \n(% Impervious Surface)", y = "Occupancy", color = "Diet \nDiversity \nIndex") +
    #gghighlight(species %in% c( "cottontail_desert", "mountain_lion", "bobcat", "jackrabbit_bt")) + 
    theme(axis.text.x = element_text(face = "bold"), 
          axis.text.y = element_text(face = "bold"), 
          # legend.position = "none",
          legend.title = element_text(face = "bold"),
          legend.text = element_text(face = "bold"),
          axis.ticks = element_blank(),
          axis.title.x = element_text(face = "bold", size = 18), 
          axis.title.y = element_text(face = "bold", size = 18))
  plot.psi.urb
  
  ggsave("./figures/traitvarideal_urbdiet_vs_1use.png",
         plot.psi.urb,
         width = 4,
         height = 4,
         units = "in",
         dpi = 300)
  
  # Persistence
  plot.phi.urb <- ggplot(data = data.pred.di, aes(x = impervious, y = phi_urb_med, 
                                                  group = dietdiv, fill = dietdiv, col = dietdiv)) +
    geom_smooth(method = "gam", lwd = 2, col = "black", se = FALSE) +
    geom_smooth(method = "gam", lwd = 1, se = FALSE) +
    theme_classic(base_size = 20) +
    geom_ribbon(aes(ymin = phi_urb_low95, ymax = phi_urb_high95), alpha = 0.4, color = NA) +
    scale_color_viridis(option = "viridis", direction = -1, end = 0.5, begin = 0.95) +
    scale_fill_viridis(option = "viridis", direction = -1, end = 0.5, begin = 0.95) +
    coord_cartesian(ylim = c(0,1)) +
    labs(x = "Urbanization \n(% Impervious Surface)", y = "Persistence", color = "Diet \nDiversity \nIndex") +
    theme(axis.text.x = element_text(face = "bold"), 
          axis.text.y = element_text(face = "bold"), 
          legend.position = "none",
          legend.title = element_text(face = "bold"),
          legend.text = element_text(face = "bold"),
          axis.ticks = element_blank(),
          axis.title.x = element_text(face = "bold", size = 18), 
          axis.title.y = element_text(face = "bold", size = 18))
  plot.phi.urb
  
  ggsave("./figures/traitvarideal_urbdiet_vs_2persistence.png",
         plot.phi.urb,
         width = 4,
         height = 4,
         units = "in",
         dpi = 300)
  
  # Colonization
  plot.gamma.urb <- ggplot(data = data.pred.di, aes(x = impervious, y = gamma_urb_med, 
                                                    group = dietdiv, fill = dietdiv, 
                                                    col = dietdiv#, linetype = as.factor(dietdiv)
                                                    )) +
    geom_smooth(method = "gam", lwd = 2, col = "black", se = FALSE) +
    geom_smooth(method = "gam", lwd = 1, se = FALSE) +
    theme_classic(base_size = 20) +
    geom_ribbon(aes(ymin = gamma_urb_low95, ymax = gamma_urb_high95), alpha = 0.4, color = NA) +
    # scale_linetype_manual(values = c("solid", "dashed")) +
    scale_color_viridis(option = "viridis", direction = -1, end = 0.5, begin = 0.95) +
    scale_fill_viridis(option = "viridis", direction = -1, end = 0.5, begin = 0.95) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Urbanization \n(% Impervious Surface)", y = "Colonization", color = "Diet \nDiversity \nIndex") +
    theme(axis.text.x = element_text(face = "bold"), 
          axis.text.y = element_text(face = "bold"), 
          legend.position = "none",
          legend.title = element_text(face = "bold"),
          legend.text = element_text(face = "bold"),
          axis.ticks = element_blank(),
          axis.title.x = element_text(face = "bold", size = 18), 
          axis.title.y = element_text(face = "bold", size = 18))
  plot.gamma.urb     # substantial relationship in model without quadratic urbanization term
  
  ggsave("./figures/traitvarideal_urbdiet_vs_3colonization.png",
         plot.gamma.urb,
         width = 4,
         height = 4,
         units = "in",
         dpi = 300)
  
  
  ##### Diet interaction with patch diversity relationships #####
  # Use/Occupancy
  plot.psi.pdiv <- ggplot(data = data.pred.di, aes(x = sdhi, y = psi_pdiv_med, 
                                                   group = dietdiv, fill = dietdiv, col = dietdiv)) +
    geom_smooth(method = "gam", lwd = 2, col = "black", se = FALSE) +
    geom_smooth(method = "gam", lwd = 1, se = FALSE) +
    theme_classic(base_size = 20) +
    geom_ribbon(aes(ymin = psi_pdiv_low95, ymax = psi_pdiv_high95), alpha = 0.4, color = NA) +
    scale_color_viridis(option = "viridis", direction = -1, end = 0.5, begin = 0.95) +
    scale_fill_viridis(option = "viridis", direction = -1, end = 0.5, begin = 0.95) +
    scale_y_continuous(labels = label_number(accuracy = 0.01)) +
    coord_cartesian(ylim = c(0,1)) +
    labs(x = "Patch Diversity \n(Shannon Diversity Index)", y = "Occupancy", color = "Diet \nDiversity \nIndex") +
    #gghighlight(species %in% c( "cottontail_desert", "mountain_lion", "bobcat", "jackrabbit_bt")) + 
    theme(axis.text.x = element_text(face = "bold"), 
          axis.text.y = element_text(face = "bold"), 
          legend.position = "none",
          legend.title = element_text(face = "bold"),
          legend.text = element_text(face = "bold"),
          axis.ticks = element_blank(),
          axis.title.x = element_text(face = "bold", size = 18), 
          axis.title.y = element_text(face = "bold", size = 18))
  plot.psi.pdiv
  
  ggsave("./figures/traitvarideal_pdivdiet_vs_1use.png",
         plot.psi.pdiv,
         width = 4,
         height = 4,
         units = "in",
         dpi = 300)
  
  # Persistence
  plot.phi.pdiv <- ggplot(data = data.pred.di, aes(x = sdhi, y = phi_pdiv_med, 
                                                   group = dietdiv, fill = dietdiv, col = dietdiv)) +
    geom_smooth(method = "gam", lwd = 2, col = "black", se = FALSE) +
    geom_smooth(method = "gam", lwd = 1, se = FALSE) +
    theme_classic(base_size = 20) +
    geom_ribbon(aes(ymin = phi_pdiv_low95, ymax = phi_pdiv_high95), alpha = 0.4, color = NA) +
    scale_color_viridis(option = "viridis", direction = -1, end = 0.5, begin = 0.95) +
    scale_fill_viridis(option = "viridis", direction = -1, end = 0.5, begin = 0.95) +
    coord_cartesian(ylim = c(0,1)) +
    labs(x = "Patch Diversity \n(Shannon Diversity Index)", y = "Persistence", color = "Diet \nDiversity \nIndex") +
    theme(axis.text.x = element_text(face = "bold"), 
          axis.text.y = element_text(face = "bold"), 
          legend.position = "none",
          legend.title = element_text(face = "bold"),
          legend.text = element_text(face = "bold"),
          axis.ticks = element_blank(),
          axis.title.x = element_text(face = "bold", size = 18), 
          axis.title.y = element_text(face = "bold", size = 18))
  plot.phi.pdiv
  
  ggsave("./figures/traitvarideal_pdivdiet_vs_2persistence.png",
         plot.phi.pdiv,
         width = 4,
         height = 4,
         units = "in",
         dpi = 300)
  
  # Colonization
  plot.gamma.pdiv <- ggplot(data = data.pred.di, aes(x = sdhi, y = gamma_pdiv_med, 
                                                     group = dietdiv, fill = dietdiv, col = dietdiv)) +
    geom_smooth(method = "gam", lwd = 2, col = "black", se = FALSE) +
    geom_smooth(method = "gam", lwd = 1, se = FALSE) +
    theme_classic(base_size = 20) +
    geom_ribbon(aes(ymin = gamma_pdiv_low95, ymax = gamma_pdiv_high95), alpha = 0.4, color = NA) +
    scale_color_viridis(option = "viridis", direction = -1, end = 0.5, begin = 0.95) +
    scale_fill_viridis(option = "viridis", direction = -1, end = 0.5, begin = 0.95) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Patch Diversity \n(Shannon Diversity Index)", y = "Colonization", color = "Diet \nDiversity \nIndex") +
    theme(axis.text.x = element_text(face = "bold"), 
          axis.text.y = element_text(face = "bold"), 
          legend.position = "none",
          legend.title = element_text(face = "bold"),
          legend.text = element_text(face = "bold"),
          axis.ticks = element_blank(),
          axis.title.x = element_text(face = "bold", size = 18), 
          axis.title.y = element_text(face = "bold", size = 18))
  plot.gamma.pdiv
  
  ggsave("./figures/traitvarideal_pdivdiet_vs_3colonization.png",
         plot.gamma.pdiv,
         width = 4,
         height = 4,
         units = "in",
         dpi = 300)
  
  ##### Diet interaction with vegetation greenness relationship #####
  # Use/Occupancy
  plot.psi.veg <- ggplot(data = data.pred.di, aes(x = veg, y = psi_veg_med, 
                                                  group = dietdiv, fill = dietdiv, col = dietdiv)) +
    geom_smooth(method = "gam", lwd = 2, col = "black", se = FALSE) +
    geom_smooth(method = "gam", lwd = 1, se = FALSE) +
    theme_classic(base_size = 20) +
    geom_ribbon(aes(ymin = psi_veg_low95, ymax = psi_veg_high95), alpha = 0.4, color = NA) +
    scale_color_viridis(option = "viridis", direction = -1, end = 0.5, begin = 0.95) +
    scale_fill_viridis(option = "viridis", direction = -1, end = 0.5, begin = 0.95) +
    scale_y_continuous(labels = label_number(accuracy = 0.01)) +
    coord_cartesian(ylim = c(0,1)) +
    labs(x = "Vegetation Greenness \n(NDVI)", y = "Occupancy", color = "Diet \nDiversity \nIndex") +
    #gghighlight(species %in% c( "cottontail_desert", "mountain_lion", "bobcat", "jackrabbit_bt")) + 
    theme(axis.text.x = element_text(face = "bold"), 
          axis.text.y = element_text(face = "bold"), 
          legend.position = "none",
          legend.title = element_text(face = "bold"),
          legend.text = element_text(face = "bold"),
          axis.ticks = element_blank(),
          axis.title.x = element_text(face = "bold", size = 18), 
          axis.title.y = element_text(face = "bold", size = 18))
  plot.psi.veg
  
  ggsave("./figures/traitvarideal_vegdiet_vs_1use.png",
         plot.psi.veg,
         width = 4,
         height = 4,
         units = "in",
         dpi = 300)
  
  # Persistence
  plot.phi.veg <- ggplot(data = data.pred.di, aes(x = veg, y = phi_veg_med, group = dietdiv, fill = dietdiv, col = dietdiv)) +
     geom_smooth(method = "gam", lwd = 2, col = "black", se = FALSE) +
     geom_smooth(method = "gam", lwd = 1, se = FALSE) +
    theme_classic(base_size = 20) +
    geom_ribbon(aes(ymin = phi_veg_low95, ymax = phi_veg_high95), alpha = 0.4, color = NA) +
    scale_color_viridis(option = "viridis", direction = -1, end = 0.5, begin = 0.95) +
    scale_fill_viridis(option = "viridis", direction = -1, end = 0.5, begin = 0.95) +
    coord_cartesian(ylim = c(0,1)) +
    labs(x = "Vegetation Greenness \n(NDVI)", y = "Persistence", color = "Diet \nDiversity \nIndex") +
    theme(axis.text.x = element_text(face = "bold"), 
          axis.text.y = element_text(face = "bold"), 
          legend.position = "none",
          legend.title = element_text(face = "bold"),
          legend.text = element_text(face = "bold"),
          axis.ticks = element_blank(),
          axis.title.x = element_text(face = "bold", size = 18), 
          axis.title.y = element_text(face = "bold", size = 18))
  plot.phi.veg
  
  ggsave("./figures/traitvarideal_vegdiet_vs_2persistence.png",
         plot.phi.veg,
         width = 4,
         height = 4,
         units = "in",
         dpi = 300)
  
  # Colonization
  plot.gamma.veg <- ggplot(data = data.pred.di, aes(x = veg, y = gamma_veg_med, 
                                                    group = dietdiv, fill = dietdiv, col = dietdiv)) +
     geom_smooth(method = "gam", lwd = 2, col = "black", se = FALSE) +
     geom_smooth(method = "gam", lwd = 1, se = FALSE) +
    theme_classic(base_size = 20) +
    geom_ribbon(aes(ymin = gamma_veg_low95, ymax = gamma_veg_high95), alpha = 0.4, color = NA) +
    scale_color_viridis(option = "viridis", direction = -1, end = 0.5, begin = 0.95) +
    scale_fill_viridis(option = "viridis", direction = -1, end = 0.5, begin = 0.95) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Vegetation Greenness \n(NDVI)", y = "Colonization", color = "Diet \nDiversity \nIndex") +
    theme(axis.text.x = element_text(face = "bold"), 
          axis.text.y = element_text(face = "bold"), 
          legend.position = "none",
          legend.title = element_text(face = "bold"),
          legend.text = element_text(face = "bold"),
          axis.ticks = element_blank(),
          axis.title.x = element_text(face = "bold", size = 18), 
          axis.title.y = element_text(face = "bold", size = 18))
  plot.gamma.veg
  
  ggsave("./figures/traitvarideal_vegdiet_vs_3colonization.png",
         plot.gamma.veg,
         width = 4,
         height = 4,
         units = "in",
         dpi = 300)
  
  
 
  
#### SUPPLEMENTARY FIGURE: PREDICTIONS FOR ACTUAL SPECIES ####
##### plot spp predictions vs urbanization #####  
  spp.plot <- c("coyote")
  # spp.plot <- c("kangaroo_rat")
  spp.plot <- names.short
  for (s in 1:length(spp.plot)){
    spp.name <- data.spp$name[which(data.spp$name_short == spp.plot[s])]
    
    plot.pred <- pred.spp.urb %>% filter(name_short == spp.plot[s]) %>% 
      ggplot(aes(x = impervious, y = psi_med)) +
      theme_classic() +
      geom_point(aes(x = impervious, y = psi_med), size = 0.5)+
      geom_smooth(method = "gam", color = "black", aes(x = impervious, y = psi_med), se = FALSE)+
      geom_ribbon(aes(x = impervious, y = psi_med, ymin = psi_low95, ymax = psi_upp95),
                  alpha = 0.4) +
      scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75)) +
      coord_cartesian(ylim=c(0,1))+
      labs(x = "Urbanization \n(% Impervious Surface)", y = "Occupancy", title = spp.name) +
      theme(axis.text.x = element_text(face = "bold", size = 14),
            axis.text.y = element_text(face = "bold", size = 14),
            axis.title.x = element_text(face = "bold", size = 18),
            axis.title.y = element_text(face = "bold", size = 18),
            legend.position = "none",
            legend.margin = margin(t = 5, r = 74, b = 5, l = 5, unit = "pt") )
    # print(plot.pred)
    
    ggsave(paste(paste("~/GitHub/metrophxwildlife-communitydynamics/figures/spp_pred/pred_1psi_urb_", spp.plot[s], sep = ""), ".png", sep = ""),
           plot.pred,
           width = 4,
           height = 4,
           units = "in",
           dpi = 300)
    
  }
  
  
  
  # Persistence
  spp.plot <- names.short
  for (s in 1:length(spp.plot)){
    spp.name <- data.spp$name[which(data.spp$name_short == spp.plot[s])]
    
    plot.pred <- pred.spp.urb %>% filter(name_short == spp.plot[s]) %>% 
      ggplot(aes(x = impervious, y = phi_med)) +
      theme_classic() +
      # geom_point(aes(x = impervious, y = psi_med))+
      geom_smooth(method = "gam", color = "black", aes(x = impervious, y = phi_med), se = FALSE)+
      geom_ribbon(aes(x = impervious, y = phi_med, ymin = phi_low95, ymax = phi_upp95),
                  alpha = 0.4) +
      scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75)) +
      coord_cartesian(ylim=c(0,1))+
      labs(x = "Urbanization \n(% Impervious Surface)", y = "Persistence", title = spp.name) +
      theme(axis.text.x = element_text(face = "bold", size = 14),
            axis.text.y = element_text(face = "bold", size = 14),
            axis.title.x = element_text(face = "bold", size = 18),
            axis.title.y = element_text(face = "bold", size = 18),
            legend.position = "none",
            legend.margin = margin(t = 5, r = 74, b = 5, l = 5, unit = "pt") )
    # print(plot.pred)
    
    ggsave(paste(paste("~/GitHub/metrophxwildlife-communitydynamics/figures/spp_pred/pred_2phi_urb_", spp.plot[s], sep = ""), ".png", sep = ""),
           plot.pred,
           width = 4,
           height = 4,
           units = "in",
           dpi = 300)
    
  }
  
  
  # Colonization
  spp.plot <- names.short
  for (s in 1:length(spp.plot)){
    spp.name <- data.spp$name[which(data.spp$name_short == spp.plot[s])]
    
    plot.pred <- pred.spp.urb %>% filter(name_short == spp.plot[s]) %>% 
      ggplot(aes(x = impervious, y = gamma_med)) +
      theme_classic() +
      # geom_point(aes(x = impervious, y = psi_med))+
      geom_smooth(method = "gam", color = "black", aes(x = impervious, y = gamma_med), se = FALSE)+
      geom_ribbon(aes(x = impervious, y = gamma_med, ymin = gamma_low95, ymax = gamma_upp95),
                  alpha = 0.4) +
      scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75)) +
      coord_cartesian(ylim=c(0,1))+
      labs(x = "Urbanization \n(% Impervious Surface)", y = "Colonization", title = spp.name) +
      theme(axis.text.x = element_text(face = "bold", size = 14),
            axis.text.y = element_text(face = "bold", size = 14),
            axis.title.x = element_text(face = "bold", size = 18),
            axis.title.y = element_text(face = "bold", size = 18),
            legend.position = "none",
            legend.margin = margin(t = 5, r = 74, b = 5, l = 5, unit = "pt") )
    # print(plot.pred)
    
    ggsave(paste(paste("~/GitHub/metrophxwildlife-communitydynamics/figures/spp_pred/pred_3gamma_urb_", spp.plot[s], sep = ""), ".png", sep = ""),
           plot.pred,
           width = 4,
           height = 4,
           units = "in",
           dpi = 300)
    
  }
  
  
  
##### plot species predictions vs patch diversity #####  
  
  # Occupancy
  spp.plot <- c("coyote")
  # spp.plot <- c("kangaroo_rat")
  spp.plot <- names.short
  for (s in 1:length(spp.plot)){
    spp.name <- data.spp$name[which(data.spp$name_short == spp.plot[s])]
    
    plot.pred <- pred.spp.pdiv %>% filter(name_short == spp.plot[s]) %>% 
      ggplot(aes(x = sdhi, y = psi_med)) +
      theme_classic() +
      # geom_point(aes(x = sdhi, y = psi_med))+
      geom_smooth(method = "gam", color = "black", aes(x = sdhi, y = psi_med), se = FALSE)+
      geom_ribbon(aes(x = sdhi, y = psi_med, ymin = psi_low95, ymax = psi_upp95),
                  alpha = 0.4) +
      scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75)) +
      coord_cartesian(ylim=c(0,1))+
      labs(x = "Patch Diversity \n(Shannon Diversity Index)", y = "Occupancy", title = spp.name) +
      theme(axis.text.x = element_text(face = "bold", size = 14),
            axis.text.y = element_text(face = "bold", size = 14),
            axis.title.x = element_text(face = "bold", size = 18),
            axis.title.y = element_text(face = "bold", size = 18),
            legend.position = "none",
            legend.margin = margin(t = 5, r = 74, b = 5, l = 5, unit = "pt") )
    # print(plot.pred)
    
    ggsave(paste(paste("~/GitHub/metrophxwildlife-communitydynamics/figures/spp_pred/pred_1psi_pdiv_", spp.plot[s], sep = ""), ".png", sep = ""),
           plot.pred,
           width = 4,
           height = 4,
           units = "in",
           dpi = 300)
    
  }
  
  # Persistence 
  spp.plot <- names.short
  for (s in 1:length(spp.plot)){
    spp.name <- data.spp$name[which(data.spp$name_short == spp.plot[s])]
    
    plot.pred <- pred.spp.pdiv %>% filter(name_short == spp.plot[s]) %>% 
      ggplot(aes(x = sdhi, y = phi_med)) +
      theme_classic() +
      # geom_point(aes(x = sdhi, y = psi_med))+
      geom_smooth(method = "gam", color = "black", aes(x = sdhi, y = phi_med), se = FALSE)+
      geom_ribbon(aes(x = sdhi, y = phi_med, ymin = phi_low95, ymax = phi_upp95),
                  alpha = 0.4) +
      scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75)) +
      coord_cartesian(ylim=c(0,1))+
      labs(x = "Patch Diversity \n(Shannon Diversity Index)", y = "Persistence", title = spp.name) +
      theme(axis.text.x = element_text(face = "bold", size = 14),
            axis.text.y = element_text(face = "bold", size = 14),
            axis.title.x = element_text(face = "bold", size = 18),
            axis.title.y = element_text(face = "bold", size = 18),
            legend.position = "none",
            legend.margin = margin(t = 5, r = 74, b = 5, l = 5, unit = "pt") )
    # print(plot.pred)
    
    ggsave(paste(paste("~/GitHub/metrophxwildlife-communitydynamics/figures/spp_pred/pred_2phi_pdiv_", spp.plot[s], sep = ""), ".png", sep = ""),
           plot.pred,
           width = 4,
           height = 4,
           units = "in",
           dpi = 300)
    
  }
  
  # Colonization 
  spp.plot <- names.short
  for (s in 1:length(spp.plot)){
    spp.name <- data.spp$name[which(data.spp$name_short == spp.plot[s])]
    
    plot.pred <- pred.spp.pdiv %>% filter(name_short == spp.plot[s]) %>% 
      ggplot(aes(x = sdhi, y = gamma_med)) +
      theme_classic() +
      # geom_point(aes(x = sdhi, y = psi_med))+
      geom_smooth(method = "gam", color = "black", aes(x = sdhi, y = gamma_med), se = FALSE)+
      geom_ribbon(aes(x = sdhi, y = gamma_med, ymin = gamma_low95, ymax = gamma_upp95),
                  alpha = 0.4) +
      scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75)) +
      coord_cartesian(ylim=c(0,1))+
      labs(x = "Patch Diversity \n(Shannon Diversity Index)", y = "Colonization", title = spp.name) +
      theme(axis.text.x = element_text(face = "bold", size = 14),
            axis.text.y = element_text(face = "bold", size = 14),
            axis.title.x = element_text(face = "bold", size = 18),
            axis.title.y = element_text(face = "bold", size = 18),
            legend.position = "none",
            legend.margin = margin(t = 5, r = 74, b = 5, l = 5, unit = "pt") )
    # print(plot.pred)
    
    ggsave(paste(paste("~/GitHub/metrophxwildlife-communitydynamics/figures/spp_pred/pred_3gamma_pdiv_", spp.plot[s], sep = ""), ".png", sep = ""),
           plot.pred,
           width = 4,
           height = 4,
           units = "in",
           dpi = 300)
    
  }  
  
##### plot spp predictions vs vegetation #####  
  # Initial occupancy
  spp.plot <- c("coyote")
  # spp.plot <- c("kangaroo_rat")
  spp.plot <- names.short
  for (s in 1:length(spp.plot)){
    spp.name <- data.spp$name[which(data.spp$name_short == spp.plot[s])]
    
    plot.pred <- pred.spp.veg %>% filter(name_short == spp.plot[s]) %>% 
      ggplot(aes(x = veg, y = psi_med)) +
      theme_classic() +
      # geom_point(aes(x = veg, y = psi_med))+
      geom_smooth(method = "gam", color = "black", aes(x = veg, y = psi_med), se = FALSE)+
      geom_ribbon(aes(x = veg, y = psi_med, ymin = psi_low95, ymax = psi_upp95),
                  alpha = 0.4) +
      scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75)) +
      coord_cartesian(ylim=c(0,1))+
      labs(x = "Vegetation Greenness \n(NDVI)", y = "Occupancy", title = spp.name) +
      theme(axis.text.x = element_text(face = "bold", size = 14),
            axis.text.y = element_text(face = "bold", size = 14),
            axis.title.x = element_text(face = "bold", size = 18),
            axis.title.y = element_text(face = "bold", size = 18),
            legend.position = "none",
            legend.margin = margin(t = 5, r = 74, b = 5, l = 5, unit = "pt") )
    # print(plot.pred)
    
    ggsave(paste(paste("~/GitHub/metrophxwildlife-communitydynamics/figures/spp_pred/pred_1psi_veg_", spp.plot[s], sep = ""), ".png", sep = ""),
           plot.pred,
           width = 4,
           height = 4,
           units = "in",
           dpi = 300)
    
  }
  
  
  # Persistence
  spp.plot <- names.short
  for (s in 1:length(spp.plot)){
    spp.name <- data.spp$name[which(data.spp$name_short == spp.plot[s])]
    
    plot.pred <- pred.spp.veg %>% filter(name_short == spp.plot[s]) %>% 
      ggplot(aes(x = veg, y = phi_med)) +
      theme_classic() +
      # geom_point(aes(x = veg, y = psi_med))+
      geom_smooth(method = "gam", color = "black", aes(x = veg, y = phi_med), se = FALSE)+
      geom_ribbon(aes(x = veg, y = phi_med, ymin = phi_low95, ymax = phi_upp95),
                  alpha = 0.4) +
      scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75)) +
      coord_cartesian(ylim=c(0,1))+
      labs(x = "Vegetation Greenness \n(NDVI)", y = "Persistence", title = spp.name) +
      theme(axis.text.x = element_text(face = "bold", size = 14),
            axis.text.y = element_text(face = "bold", size = 14),
            axis.title.x = element_text(face = "bold", size = 18),
            axis.title.y = element_text(face = "bold", size = 18),
            legend.position = "none",
            legend.margin = margin(t = 5, r = 74, b = 5, l = 5, unit = "pt") )
    # print(plot.pred)
    
    ggsave(paste(paste("~/GitHub/metrophxwildlife-communitydynamics/figures/spp_pred/pred_2phi_veg_", spp.plot[s], sep = ""), ".png", sep = ""),
           plot.pred,
           width = 4,
           height = 4,
           units = "in",
           dpi = 300)
    
  }
  
  
  # Colonization
  spp.plot <- names.short
  for (s in 1:length(spp.plot)){
    spp.name <- data.spp$name[which(data.spp$name_short == spp.plot[s])]
    
    plot.pred <- pred.spp.veg %>% filter(name_short == spp.plot[s]) %>% 
      ggplot(aes(x = veg, y = gamma_med)) +
      theme_classic() +
      # geom_point(aes(x = veg, y = psi_med))+
      geom_smooth(method = "gam", color = "black", aes(x = veg, y = gamma_med), se = FALSE)+
      geom_ribbon(aes(x = veg, y = gamma_med, ymin = gamma_low95, ymax = gamma_upp95),
                  alpha = 0.4) +
      scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75)) +
      coord_cartesian(ylim=c(0,1))+
      labs(x = "Vegetation Greenness \n(NDVI)", y = "Colonization", title = spp.name) +
      theme(axis.text.x = element_text(face = "bold", size = 14),
            axis.text.y = element_text(face = "bold", size = 14),
            axis.title.x = element_text(face = "bold", size = 18),
            axis.title.y = element_text(face = "bold", size = 18),
            legend.position = "none",
            legend.margin = margin(t = 5, r = 74, b = 5, l = 5, unit = "pt") )
    # print(plot.pred)
    
    ggsave(paste(paste("~/GitHub/metrophxwildlife-communitydynamics/figures/spp_pred/pred_3gamma_veg_", spp.plot[s], sep = ""), ".png", sep = ""),
           plot.pred,
           width = 4,
           height = 4,
           units = "in",
           dpi = 300)
    
  }
  
  