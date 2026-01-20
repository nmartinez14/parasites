## Script for plotting all of the important explanatory variables.
rm(list=ls())
source("lab_paths.R")
setwd(local.path)
## Packages required
library(ggpubr)
library(gridExtra)
library(tidyverse)
library(brms)
library(tidybayes)
setwd("parasites/analysis")
load(file="saved/spec_weights.Rdata")
source("src/misc.R")
source("src/ggplotThemes.R")


site.screened <- spec.orig %>%
  group_by(Site, Year, SampleRound) %>%
  summarise(SiteScreened = sum(!is.na(Apidae)))

spec.orig <- left_join(spec.orig, site.screened)
spec.net <- left_join(spec.net, site.screened)
## ***********************************************************************
## scaling/unscaling labs
## ***********************************************************************
# original data, subsetted to unique values
spec.uni.orig <- spec.orig[spec.orig$Weights ==1,]
## scaled data subsetted to unique values
spec.uni <- spec.net[spec.net$Weights ==1,]

## use unscaled data to have nice axis labels, convert to scaled for
## the axes

## lat (logged)
labs.lat.x <- pretty(c(spec.uni.orig$Lat),
                      n=10)
axis.lat.x <-  standardize.axis(labs.lat.x, spec.uni.orig$Lat)

## flower div (not logged)
labs.flower.div <- (pretty(spec.uni.orig$MeanFloralDiversity, n=5))
axis.flower.div <-  standardize.axis(labs.flower.div,
                                     spec.uni.orig$MeanFloralDiversity)
## HB abund (logged + 1)
labs.HB.abund <- (pretty(c(spec.uni.orig$Net_HBAbundance), n=5))
axis.HB.abund <-  standardize.axis(labs.HB.abund, spec.uni.orig$Net_HBAbundance)
## bombus abund (logged + 1)
labs.bombus.abund <- (pretty(c(spec.uni.orig$Net_BombusAbundance), n=5))
axis.bombus.abund <-  standardize.axis(labs.bombus.abund, spec.uni.orig$Net_BombusAbundance)

## bee diversity (not logged)
labs.bee.div <- (pretty(c(spec.uni.orig$Net_BeeDiversity), n=5))
axis.bee.div <-  standardize.axis(labs.bee.div,
                                  spec.uni.orig$Net_BeeDiversity)

## use all the species data or just bombus? 
bombus.par <- spec.orig[spec.orig$WeightsPar==1 & spec.orig$Genus == "Bombus", ]
apis.par <- spec.orig[spec.orig$WeightsPar==1 & spec.orig$Genus == "Apis", ]
## rare.degree (logged)
labs.degree <- (pretty(bombus.par$rare.degree, n=10))
axis.degree <-  standardize.axis(labs.degree,
                                  bombus.par$rare.degree)

## ***************************************************************************
# Load model for bee diversity
load(file="saved/parasiteFit_Bombus_CrithidiaPresenceApicystisSpp_bee_div_lat.Rdata")
fit.bombus.bd <- fit.parasite.bombus

## Generate newdata draws

cond.effects <- conditional_effects(fit.bombus.bd)

## ***************************************************************************
## Crithidia ~ bee diversity

crithidia_beediv <-
  cond.effects[["CrithidiaPresence.CrithidiaPresence_Net_BeeDiversity"]]

p1.parasite <- ggplot(crithidia_beediv, aes(x = Net_BeeDiversity, 
                                            y= estimate__)) +
  geom_line(aes(x = Net_BeeDiversity, y= estimate__), color = "#3182bd", 
            size = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = "#3182bd",
              alpha=0.4)+
  labs(x = "Bee diversity", y = "Crithidia prevalence",
       fill = "Credible interval") +
  #theme_dark_black() +
  theme_ms() +
  #theme(legend.position = "bottom") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=16),
        text = element_text(size=16))+
  geom_jitter(data=spec.uni,
              aes(y= CrithidiaParasitismRate, x=Net_BeeDiversity,
                  color = SiteScreened),
              width=0.05) +
  scale_color_gradient(low = "grey80", high = "grey20") +
  labs(color = "Screened individuals")


################################################################################
## Apicystis ~ bee diversity Bombus
################################################################################

apicystis_beediv <-cond.effects[["ApicystisSpp.ApicystisSpp_Net_BeeDiversity"]]

p2.parasite <- ggplot(apicystis_beediv, aes(x = Net_BeeDiversity, 
                                            y = estimate__)) +
  geom_line(aes(x = Net_BeeDiversity, y= estimate__), 
            linewidth = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__),
               alpha=0.2)+
  labs(x = "Bee diversity", y = "Apicystis prevalence") +
  #theme_dark_black()+
  theme_ms() +
  #theme(legend.position = "bottom") +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16))+
  geom_jitter(data=spec.uni,
              aes(y= ApicystisParasitismRate, x=Net_BeeDiversity,
                  color = SiteScreened),
              width=0.05) +
  scale_color_gradient(low = "grey80", high = "grey20") +
  labs(color = "Screened individuals")


## ***************************************************************************
# Load model for floral diversity
load(file="saved/parasiteFit_bombus_CrithidiaPresenceApicystisSpp_floral_div_lat.Rdata")
fit.bombus.fd <- fit.parasite.bombus

## Generate newdata draws

cond.effects <- conditional_effects(fit.bombus.fd)

## ***************************************************************************
## Crithidia ~ floral diversity


crithidia_floraldiv <-cond.effects[["CrithidiaPresence.CrithidiaPresence_MeanFloralDiversity"]]


p3.parasite <- ggplot(crithidia_floraldiv, aes(x = MeanFloralDiversity, 
                                               y= estimate__)) +
  geom_line(aes(x = MeanFloralDiversity, y= estimate__),
                color = "#3182bd", size = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = "#3182bd",
              alpha=0.4) +
    labs(x = "Floral diversity", y = "Crithidia prevalence") +
    #theme_dark_black() +
  theme_ms() +
    #theme(legend.position = "bottom") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size=16),
          text = element_text(size=16)) +
  geom_jitter(data=spec.uni,
              aes(y= CrithidiaParasitismRate, x=MeanFloralDiversity,
                  color = SiteScreened),
              width=0.05) +
  scale_color_gradient(low = "grey80", high = "grey20") +
  labs(color = "Screened individuals")


## ***************************************************************************
## Apicystis ~ floral diversity

apicystis_floraldiv <- 
  cond.effects[["ApicystisSpp.ApicystisSpp_MeanFloralDiversity"]]


p4.parasite <- ggplot(apicystis_floraldiv, aes(x = MeanFloralDiversity, 
                                               y = estimate__)) +
  geom_line(aes(x = MeanFloralDiversity, y= estimate__), 
            size = 1.5, color = "#3182bd") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), 
              alpha=0.4, fill = "#3182bd")+
  labs(x = "Floral diversity", y = "Apicystis prevalence") +
  #theme_dark_black()+
  theme_ms() +
  #theme(legend.position = "bottom") +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  geom_jitter(data=spec.uni,
              aes(y= ApicystisParasitismRate, x=MeanFloralDiversity, 
                  color = SiteScreened),
              width=0.05) +
  scale_color_gradient(low = "grey80", high = "grey20") +
  labs(color = "Screened individuals")
    
    
    
parasite.dilution <- ggarrange(p3.parasite, p1.parasite, p4.parasite, p2.parasite,
                            labels = c("A", "B", "C","D"),
                            #font.label = list(color = "white"),
                            ncol = 2, nrow = 2,
                            common.legend = TRUE,
                            legend = "bottom")

ggsave(parasite.dilution, file="figures/parasite_diversity.pdf", 
       height=6, width=10)

## ***************************************************************************
# Load model for bombus abundance
load(file="saved/parasiteFit_Bombus_CrithidiaPresenceApicystisSpp_bombus_abund_lat.Rdata")
fit.bombus.ba <- fit.parasite.bombus

## Generate newdata draws

cond.effects <- conditional_effects(fit.bombus.ba)
## ***************************************************************************
## crithidia ~ bombus abundance

crithidia_beeabun <-
  cond.effects[["CrithidiaPresence.CrithidiaPresence_Net_BombusAbundance"]]

p5.parasite <- ggplot(crithidia_beeabun, aes(x = Net_BombusAbundance, 
                                             y = estimate__)) +
  geom_line(aes(x = Net_BombusAbundance, y= estimate__), 
            size = 1.5, color = "darkgoldenrod3") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha=0.4, 
              fill = "darkgoldenrod3") +
  scale_fill_manual(labels ="Bombus 0.95")+
  labs(x = "Bombus abundance (log)", y = "Crithidia prevalence") +
  #theme_dark_black()+
  theme_ms() +
  #theme(legend.position = "bottom") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
   geom_jitter(data=spec.uni,
               aes(y= CrithidiaParasitismRate, x=Net_BombusAbundance, 
                   color = SiteScreened),
               width=0.05) +
  scale_color_gradient(low = "grey80", high = "grey20") +
  labs(color = "Screened individuals")



################################################################################
## Apicysits ~ bombus abundance
################################################################################
apicystis_bombusabun <-
  cond.effects[["ApicystisSpp.ApicystisSpp_Net_BombusAbundance"]]

p6.parasite <- ggplot(apicystis_bombusabun, aes(x = Net_BombusAbundance, 
                                                y = estimate__)) +
  geom_line(aes(x = Net_BombusAbundance, y= estimate__), 
            size = 1.5, color = "#3182bd") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), 
              alpha=0.4, fill = "#3182bd")+
  labs(x = "Bombus abundance (log)", y = "Apicystis prevalence") +
  #theme_dark_black()+
  theme_ms() +
  #theme(legend.position = "bottom") +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16))+  
  geom_jitter(data=spec.uni,
              aes(y= ApicystisParasitismRate, x=Net_BombusAbundance, 
                  color = SiteScreened),
              width=0.05) +
  scale_color_gradient(low = "grey80", high = "grey20") +
  labs(color = "Screened individuals")


## ***************************************************************************
# Load model for apis abundance
load(file="saved/parasiteFit_Bombus_CrithidiaPresenceApicystisSpp_hb_abund_lat.Rdata")
fit.bombus.ha <- fit.parasite.bombus

## Generate newdata draws
cond.effects <- conditional_effects(fit.bombus.ha)

################################################################################
## Crithidia ~ apis abundance
################################################################################

crithidia_hbabun <-
  cond.effects[["CrithidiaPresence.CrithidiaPresence_Net_HBAbundance"]]

p7.parasite <- ggplot(crithidia_hbabun, aes(x = Net_HBAbundance, 
                                            y = estimate__)) +
  geom_line(aes(x = Net_HBAbundance, y= estimate__), size = 1.5, 
            color = "#3182bd") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = "#3182bd", 
              alpha=0.4)+
  labs(x = "Apis abundance (log)", y = "Crithidia prevalence") +
  #theme_dark_black()+
  theme_ms() +
  #theme(legend.position = "bottom") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) + 
  geom_jitter(data=spec.uni,
              aes(y= CrithidiaParasitismRate, x=Net_HBAbundance, 
                  color = SiteScreened),
              width=0.05) +
  scale_color_gradient(low = "grey80", high = "grey20") +
  labs(color = "Screened individuals")


################################################################################
## Apicysits ~ apis abundance
################################################################################

apicystis_hbabun <-
  cond.effects[["ApicystisSpp.ApicystisSpp_Net_HBAbundance"]]

p8.parasite <- ggplot(apicystis_hbabun, aes(x = Net_HBAbundance, 
                                            y = estimate__)) +
  geom_line(aes(x = Net_HBAbundance, y= estimate__), size = 1.5, 
            color = "#3182bd") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = "#3182bd", 
              alpha=0.4)+
  labs(x = "Apis abundance (log)", y = "Apicystis prevalence") +
  #theme_dark_black()+
  theme_ms() +
  #theme(legend.position = "bottom") +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
         text = element_text(size=16)) +
  geom_jitter(data=spec.uni,
              aes(y= ApicystisParasitismRate, x=Net_HBAbundance, 
                  color = SiteScreened),
              width=0.05) +
  scale_color_gradient(low = "grey80", high = "grey20") +
  labs(color = "Screened individuals")




parasite.amplification <- ggarrange(p5.parasite, p7.parasite, p6.parasite, p8.parasite,
                                    nrow= 2, ncol = 2,
                                    #font.label = list(color = "white"),
                                    labels = c("A", "B", "C", "D"),
                                    common.legend = TRUE,
                                    legend = "bottom")

ggsave(parasite.amplification, file="figures/parasite_amplification.pdf",
       height=6, width=10)




