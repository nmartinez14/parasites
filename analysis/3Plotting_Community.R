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

site.screened <- spec.orig %>%
  group_by(Site, Year, SampleRound) %>%
  summarise(SiteScreened = sum(!is.na(Apidae)))

spec.net <- left_join(spec.net, site.screened)
## data subsetted to unique values
spec.uni <- spec.net[spec.net$Weights ==1,]


## ***************************************************************************
# Load model for latitude
load(file="saved/parasiteFit_Bombus_CrithidiaPresenceApicystisSpp_lat_lat.Rdata")
fit.bombus.l <- fit.parasite.bombus

## Generate newdata draws
cond.effects <- conditional_effects(fit.bombus.l)

## Community level visuals
## ***********************************************************************
## bee community diversity and latitude
## ***********************************************************************

lat_beediv <-
  cond.effects[["scaleNetBeeDiversity.scaleNetBeeDiversity_Lat"]]

p1 <- ggplot(lat_beediv, aes(x = Lat, y = estimate__)) +
  geom_line(aes(x = Lat, y= estimate__), size = 1.5, 
            color = "darkgoldenrod3") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), 
              alpha=0.4, fill = "darkgoldenrod3") +
  labs(x = "Latitude (log)", y = "Bee diversity")+
  #theme_dark_black()+
  theme_ms() +
  #theme(legend.position = "bottom") +
  theme(axis.title.x = element_text(size=16),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=12),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  geom_jitter(data= spec.uni,
             aes(y= scale(Net_BeeDiversity), x=Lat), cex=2)+
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 6))+
  coord_cartesian(
    xlim = range(spec.uni$Lat, na.rm = TRUE)
  )



################################################################################
## Plant community diversity and latitude
################################################################################

lat_floraldiv <-
  cond.effects[["scaleMeanFloralDiversity.scaleMeanFloralDiversity_Lat"]]

p2 <- ggplot(lat_floraldiv, aes(x = Lat, y = estimate__)) +
  geom_line(aes(x = Lat, y= estimate__), size = 1.5, 
            color = "#3182bd") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), 
              alpha=0.4, fill = "#3182bd") +
  labs(x = "Latitude (log)", y = "Mean floral diversity") +
  #theme_dark_black() +
  theme_ms() +
  #theme(legend.position = "bottom") +
  theme(axis.title.x = element_text(size=16),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=12),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  geom_jitter(data=spec.uni,
             aes(y= scale(MeanFloralDiversity), x=Lat), cex=2)+
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 6))+
  coord_cartesian(
    xlim = range(spec.uni$Lat, na.rm = TRUE)
  )



################################################################################
## Bee abundance and lat
################################################################################

lat_bombusabund <-
  cond.effects[["scaleNetBombusAbundance.scaleNetBombusAbundance_Lat"]]


p3 <- ggplot(lat_bombusabund, aes(x = Lat, y = estimate__)) +
  geom_line(aes(x = Lat, y= estimate__), size = 1.5, color = "#3182bd") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha=0.4, fill = "#3182bd") +
  labs(x = "Latitude (log)", y = "Bombus abundance(log)") +
  #theme_dark_black() +
  theme_ms() +
  #theme(legend.position = "bottom") +
  theme(axis.title.x = element_text(size=16),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=12),
        axis.title.y = element_text(size=14),
        text = element_text(size=16)) +
  geom_jitter(data=spec.uni,
             aes(y= scale(Net_BombusAbundance), x=Lat), cex=2)+
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 6))+
  coord_cartesian(
    xlim = range(spec.uni$Lat, na.rm = TRUE)
  )



## Honeybee abundance
lat_apisabund <- cond.effects[["scaleNetHBAbundance.scaleNetHBAbundance_Lat"]]

p4 <- ggplot(lat_apisabund, aes(x = Lat, y = estimate__)) +
  geom_line(aes(x = Lat, y= estimate__), size = 1.5) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha=0.4) +
  labs(x = "Latitude (log)", y = "Apis abundance (log)") +
  #theme_dark_black() +
   theme_ms() +
  #theme(legend.position = "bottom") +
  theme(axis.title.x = element_text(size=16),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=12),
        axis.title.y = element_text(size=15),
        text = element_text(size=16)) +
  geom_jitter(data=spec.uni,
             aes(y= scale(Net_HBAbundance), x=Lat), cex=2)+
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 6))+
  coord_cartesian(
    xlim = range(spec.uni$Lat, na.rm = TRUE)
  )




# lat_community <- ggarrange(p2,p1,p3,p4, #plots that are going to be included in this multipanel figure
#                        labels = c("A", "B", "C","D"), #labels given each panel 
#                        ncol = 2, nrow = 2, #adjust plot space 
#                        font.label = list(color = "white"))

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
              alpha=0.2) +
  labs(x = "Latitude (log)", y = "Crithidia prevalence") +
  #theme_dark_black() +
  theme_ms() +
  #theme(legend.position = "bottom") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=12),
        axis.title.y = element_text(size=16),
        text = element_text(size=16),
        plot.title = element_text(color = "black")) +
  geom_jitter(data= spec.uni,aes(y= CrithidiaParasitismRate, x= Lat,
                                 color = SiteScreened),
             width=0.05) +
  scale_color_gradient(low = "grey80", high = "grey20") +
  labs(color = "Screened individuals")+
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 6))+
  coord_cartesian(
    xlim = range(crithidia_lat$Lat, na.rm = TRUE)
  )

## ***********************************************************************
## lat and apicystis
## ***********************************************************************

apicystis_lat <-
  cond.effects[["ApicystisSpp.ApicystisSpp_Lat"]]

p6 <- ggplot(apicystis_lat, aes(x = Lat, y= estimate__)) +
  geom_line(aes(x = Lat, y= estimate__), size = 1.5, color = "#3182bd") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = "#3182bd",
              alpha=0.4) +
  labs(x = "Latitude (log)", y = "Apicystis prevalence") +
  #theme_dark_black() +
  theme_ms() +
  #theme(legend.position = "bottom") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=12),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  geom_jitter(data= spec.uni,
              aes(y= ApicystisParasitismRate, x= Lat,
                  color = SiteScreened),
              width=0.05) +
  scale_color_gradient(low = "grey80", high = "grey20") +
  labs(color = "Screened individuals")+
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 6))+
  coord_cartesian(
    xlim = range(crithidia_lat$Lat, na.rm = TRUE)
  )



################################################################################
## crithida ~ lat in Apis
## ***************************************************************************
# Load model for lat
load(file="saved/parasiteFit_Apis_CrithidiaPresenceApicystisSpp_lat_lat.Rdata")
fit.apis.l <- fit.parasite.apis

## Generate newdata draws
cond.effects <- conditional_effects(fit.apis.l)

crithidia_lat_apis <-
  cond.effects[["CrithidiaPresence.CrithidiaPresence_Lat"]]

p7 <- ggplot(crithidia_lat_apis, aes(x = Lat, y = estimate__)) +
  geom_line(aes(x = Lat, y= estimate__), size= 1.5, 
            linetype = "dotdash") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha=0.2) +
  labs(x = "Latitude (log)", y = "Crithidia prevalence") +
  #theme_dark_black() +
  theme_ms() +
  #theme(legend.position = "bottom") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=12),
        axis.title.y = element_text(size=16),
        text = element_text(size=16), 
        plot.title = element_text(color = "black")) +
  geom_jitter(data= spec.uni,
              aes(y= CrithidiaParasitismRate, x= Lat,
                  color = SiteScreened),
              width=0.05) +
  scale_color_gradient(low = "grey80", high = "grey20") +
  labs(color = "Screened individuals") +
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 6))+
  coord_cartesian(
    xlim = range(crithidia_lat_apis$Lat, na.rm = TRUE)
  )



################################################################################
## Lat and apicystis in apis
################################################################################
apicystis_lat_apis <-
  cond.effects[["ApicystisSpp.ApicystisSpp_Lat"]]

p8 <- ggplot(apicystis_lat_apis, aes(x = Lat, y = estimate__)) +
  geom_line(aes(x = Lat, y= estimate__), size = 1.5, 
            linetype = "dotdash", color = "#3182bd") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), 
              alpha=0.4, fill = "#3182bd") +
  labs(x = "Latitude (log)", y = "Apicystis prevalence") +
  #theme_dark_black() +
   theme_ms() +
  #theme(legend.position = "bottom") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=12),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  geom_jitter(data=spec.uni,
              aes(y= ApicystisParasitismRate, x= Lat,
                  color = SiteScreened),
              width=0.05) +
  scale_color_gradient(low = "grey80", high = "grey20") +
  labs(color = "Screened individuals")+
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 6))+
  coord_cartesian(
    xlim = range(crithidia_lat_apis$Lat, na.rm = TRUE)
  )




lat_full <- ggarrange(p5, p7, p6, p8,
                         p2,p1,p3,p4,
                         labels = c("A", "B", "C","D", "E", "F", "G", "H"),
                         ncol = 4, nrow = 2, 
                         heights = 1, widths = c(3, 3, 3, 3))

ggsave(lat_full, file="figures/lat_full.pdf",
       height=6, width=14)

# lat_par <- ggarrange(p5,p7,p6,p8, 
#                            labels = c("A", "B", "C","D"),  
#                            ncol = 2, nrow = 2,  
#                            font.label = list(color = "white"))
# 
# ggsave(lat_par, file="figures/lat_par.jpg",
#        height=8, width=12)

################################################################################
## Plant diversity and bee diversity
################################################################################
# Load model for floral div
load(file="saved/parasiteFit_Bombus_CrithidiaPresenceApicystisSpp_bee_div_fd.Rdata")
fit.fd <- fit.parasite.bombus

## Generate newdata draws
cond.effects <- conditional_effects(fit.fd)


beediv_floraldiv <-
  cond.effects[["scaleNetBeeDiversity.scaleNetBeeDiversity_MeanFloralDiversity"]]

plantdiv_beediv <- ggplot(beediv_floraldiv, aes(x = MeanFloralDiversity, 
                                                y = estimate__)) +
  geom_line(aes(x = MeanFloralDiversity, y= estimate__), 
            size = 1.5, color = "#3182bd") +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), 
              alpha=0.4, fill = "#3182bd") +
  labs(y = "Bee Species Diversity", x = "Floral Diversity",
       fill = "Credible interval") +
  theme(legend.position = "bottom") +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        text = element_text(size=16)) +
  theme_ms() +
  ##theme_dark_black()+
  geom_point(data=spec.uni,
             aes(y= scale(Net_BeeDiversity), x= MeanFloralDiversity), cex=2)

ggsave(plantdiv_beediv, file="figures/Beediv_floraldiv.pdf",
       height=4, width=5)

