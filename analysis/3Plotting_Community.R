## Script for plotting all of the community level explanatory variables.
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
## all bee abund (logged)
labs.bee.abund <- (pretty(c(spec.uni.orig$Net_BeeAbundance), n=5))
axis.bee.abund <-  standardize.axis(labs.bee.abund, spec.uni.orig$Net_BeeAbundance)
## bee diversity (not logged)
labs.bee.div <- (pretty(c(spec.uni.orig$Net_BeeDiversity), n=5))
axis.bee.div <-  standardize.axis(labs.bee.div,
                                  spec.uni.orig$Net_BeeDiversity)


## ***************************************************************************
# Load model for latitude
load(file="saved/parasiteFit_bombus_CrithidiaPresenceApicystisSpp_lat.Rdata")
fit.bombus.l <- fit.parasite

## Generate newdata draws
cond.effects <- conditional_effects(fit.bombus.l)

## Community level visuals
## ***********************************************************************
## bee community diversity and latitude
## ***********************************************************************

lat_beediv <-
  cond.effects[["NetBeeDiversity.NetBeeDiversity_Lat"]]

p1 <- ggplot(lat_beediv, aes(x = Lat, y = estimate__)) +
  geom_line(aes(x = Lat, y= estimate__), size = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha=0.4) +
  labs(x = "Latitude (log)", y = "Bee diversity",
       fill = "Credible interval")+
  theme_ms() +
  #theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.lat.x,
    labels =  labs.lat.x) +
  scale_y_continuous(
    breaks = axis.bee.div,
    labels =  labs.bee.div) +
  theme(axis.title.x = element_text(size=16),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=12),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  geom_point(data= spec.uni,
             aes(y= Net_BeeDiversity, x=Lat), cex=2)

ggsave(p1, file="figures/Lat_beediv.pdf",
       height=4, width=5)

################################################################################
## Plant community diversity and latitude
################################################################################

lat_floraldiv <-
  cond.effects[["MeanFloralDiversity.MeanFloralDiversity_Lat"]]

p2 <- ggplot(lat_floraldiv, aes(x = Lat, y = estimate__)) +
  geom_line(aes(x = Lat, y= estimate__), size = 1.5, 
            color = "darkgoldenrod3") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), 
              alpha=0.4, fill = "darkgoldenrod3") +
  labs(x = "Latitude (log)", y = "Mean floral diversity",
       fill = "Credible interval") +
  theme_ms() +
  #theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.lat.x,
    labels =  labs.lat.x) +
  scale_y_continuous(
    breaks = axis.flower.div,
    labels =  labs.flower.div) +
  theme(axis.title.x = element_text(size=16),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=12),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  ##theme_dark_black()+
  geom_point(data=spec.uni,
             aes(y= MeanFloralDiversity, x=Lat), cex=2)

ggsave(p2, file="figures/Lat_floraldiv.pdf",
       height=4, width=5)


################################################################################
## Bee abundance and lat
################################################################################

lat_bombusabund <-
  cond.effects[["NetBombusAbundance.NetBombusAbundance_Lat"]]


p3 <- ggplot(lat_bombusabund, aes(x = Lat, y = estimate__)) +
  geom_line(aes(x = Lat, y= estimate__), size = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha=0.4) +
  labs(x = "Latitude (log)", y = "Bombus abundance(log)",
       fill = "Credible interval") +
  theme_ms() +
  #theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.lat.x,
    labels =  labs.lat.x) +
  scale_y_continuous(
    breaks = axis.bombus.abund,
    labels =  labs.bombus.abund) +
  theme(axis.title.x = element_text(size=16),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=12),
        axis.title.y = element_text(size=14),
        text = element_text(size=16)) +
  ##theme_dark_black()+
  geom_point(data=spec.uni,
             aes(y= Net_BombusAbundance, x=Lat), cex=2)



## Honeybee abundance
lat_apisabund <-
  cond.effects[["NetHBAbundance.NetHBAbundance_Lat"]]

p4 <- ggplot(lat_apisabund, aes(x = Lat, y = estimate__)) +
  geom_line(aes(x = Lat, y= estimate__), size = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha=0.4) +
  labs(x = "Latitude (log)", y = "Apis abundance (log)",
       fill = "Credible interval") +
  theme_ms() +
  #theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.lat.x,
    labels =  labs.lat.x) +
  scale_y_continuous(
    breaks = axis.HB.abund,
    labels =  labs.HB.abund) +
  theme(axis.title.x = element_text(size=16),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=12),
        axis.title.y = element_text(size=15),
        text = element_text(size=16)) +
  ##theme_dark_black()+
  geom_point(data=spec.uni,
             aes(y= Net_HBAbundance, x=Lat), cex=2)



lat_community <- ggarrange(p2,p1,p3,p4, #plots that are going to be included in this multipanel figure
                       labels = c("A", "B", "C","D"), #labels given each panel 
                       ncol = 2, nrow = 2 #adjust plot space 
                       )
ggsave(lat_community, file="figures/lat_community.pdf",
       height=8, width=12)

## ***********************************************************************
## lat and crithidia
## ***********************************************************************

crithidia_lat <-
  cond.effects[["CrithidiaPresence.CrithidiaPresence_Lat"]]

p5 <- ggplot(crithidia_lat, aes(x = Lat, y = estimate__)) +
  geom_line(aes(x = Lat, y= estimate__), size = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__),
              alpha=0.4) +
  scale_fill_manual(labels ="Bombus 0.95") +
  labs(x = "Latitude (log)", y = "Crithidia prevalence") +
  theme_ms() +
  #theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.lat.x,
    labels =  labs.lat.x) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=12),
        axis.title.y = element_text(size=16),
        text = element_text(size=16),
        plot.title = element_text(color = "black")) +
  geom_point(data= spec.uni,aes(y= CrithidiaParasitismRate, x= Lat), cex=2)




## ***********************************************************************
## lat and apicystis
## ***********************************************************************

apicystis_lat <-
  cond.effects[["ApicystisSpp.ApicystisSpp_Lat"]]

p6 <- ggplot(apicystis_lat, aes(x = Lat, y= estimate__)) +
  geom_line(aes(x = Lat, y= estimate__), size = 1.5, color = "#3182bd") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = "#3182bd",
              alpha=0.4) +
  scale_fill_manual(labels ="Bombus 0.95") +
  labs(x = "Latitude (log)", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  #theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.lat.x,
    labels =  labs.lat.x) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=12),
        axis.title.y = element_text(size=16),
        text = element_text(size=16))+
  geom_jitter(data= spec.uni,
              aes(y= ApicystisParasitismRate, x= Lat), cex=2) 



################################################################################
## crithida ~ lat in Apis
## ***************************************************************************
# Load model for lat
load(file="saved/parasiteFit_apis_CrithidiaPresenceApicystisSpp_lat.Rdata")
fit.apis.l <- fit.parasite

## Generate newdata draws
cond.effects <- conditional_effects(fit.apis.l)

crithidia_lat_apis <-
  cond.effects[["CrithidiaPresence.CrithidiaPresence_Lat"]]

p7 <- ggplot(crithidia_lat_apis, aes(x = Lat, y = estimate__)) +
  geom_line(aes(x = Lat, y= estimate__), size= 1.5, 
            linetype = "dotdash") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha=0.4 
              ) +
  scale_fill_manual(labels ="Apis 0.95") +
  labs(x = "Latitude (log)", y = "Crithidia prevalence") +
  theme_ms() +
  #theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.lat.x,
    labels =  labs.lat.x) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=12),
        axis.title.y = element_text(size=16),
        text = element_text(size=16), 
        plot.title = element_text(color = "black")) +
  geom_jitter(data= spec.uni,
              aes(y= CrithidiaParasitismRate, x= Lat), cex=2) 



################################################################################
## Lat and apicystis in apis
################################################################################
apicystis_lat_apis <-
  cond.effects[["ApicystisSpp.ApicystisSpp_Lat"]]

p8 <- ggplot(apicystis_lat_apis, aes(x = Lat, y = estimate__)) +
  geom_line(aes(x = Lat, y= estimate__), size = 1.5, 
            linetype = "dotdash", color = "darkgoldenrod3") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), 
              alpha=0.4, fill = "darkgoldenrod3") +
  scale_fill_manual(labels ="Apis 0.95") +
  labs(x = "Latitude (log)", y = "Apicystis prevalence",
       fill = "Credible interval") +
  theme_ms() +
  #theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.lat.x,
    labels =  labs.lat.x) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=12),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  geom_jitter(data=spec.uni,
              aes(y= ApicystisParasitismRate, x= Lat), cex=2) 




lat_parasite <- ggarrange(p5, p7, p6, p8,
                         p2,p1,p3,p4,
                         labels = c("A", "B", "C","D", "E", "F", "G", "H"),
                         ncol = 4, nrow = 2, 
                         heights = 1, widths = c(3, 3, 3, 3))

ggsave(lat_parasite, file="figures/lat_parasites.pdf",
       height=6, width=14)
################################################################################
## Plant diversity and bee diversity
################################################################################
beediv_floraldiv <-
  cond.effects[["NetBeeDiversity.NetBeeDiversity_MeanFloralDiversity"]]

plantdiv_beediv <- ggplot(beediv_floraldiv, aes(x = MeanFloralDiversity, 
                                                y = estimate__)) +
  geom_line(aes(x = MeanFloralDiversity, y= estimate__), 
            size = 1.5, color = "#3182bd") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), 
              alpha=0.4, fill = "#3182bd") +
  labs(y = "Bee Species Diversity", x = "Floral Diversity",
       fill = "Credible interval") +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = axis.flower.div,
    labels =  labs.flower.div) +
  scale_y_continuous(
    breaks = axis.bee.div,
    labels =  labs.bee.div) +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  theme_ms() +
  ##theme_dark_black()+
  geom_point(data=spec.uni,
             aes(y= Net_BeeDiversity, x= MeanFloralDiversity), cex=2)

ggsave(plantdiv_beediv, file="figures/Beediv_floraldiv.pdf",
       height=4, width=5)

