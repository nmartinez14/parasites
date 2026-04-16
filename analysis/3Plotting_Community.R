## Script for plotting all of the community level explanatory variables.
rm(list=ls())
setwd("C:/Users/na_ma")
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
source("src/plotMod.R")

site.screened <- spec.net %>%
  group_by(Site, Year, SampleRound) %>%
  summarise(SiteScreened = sum(!is.na(Apidae)))

spec.net <- left_join(spec.net, site.screened)
## data subsetted to unique values
spec.uni <- spec.net[spec.net$Weights ==1,]


## ***************************************************************************
# Load model for latitude
load(file="saved/parasiteFit_Bombus_CrithidiaPresenceApicystisSpp_lat_lat.Rdata")
fit.bombus.l <- fit.parasite.bombus


## Community level visuals
## ***********************************************************************
## bee community diversity and latitude
## ***********************************************************************

p1 <- plot_cond_effects(fit.bombus.l, data = spec.uni,
                        this.response = "scaleNetBeeDiversity",
                        this.effect = "Lat",
                        significance = "95",
                        dat.x = "Lat",
                        dat.y = "Net_BeeDiversity",
                        y.label = "Bee diversity",
                        x.label = "Latitude (log)",
                        parasite = FALSE,
                        scale.y = TRUE,
                        angle.x = TRUE,
                        x.axis.lab = FALSE)


################################################################################
## Plant community diversity and latitude
################################################################################
p2 <- plot_cond_effects(fit.bombus.l, data = spec.uni,
                        this.response = "scaleMeanFloralDiversity",
                        this.effect = "Lat",
                        significance = "97",
                        dat.x = "Lat",
                        dat.y = "MeanFloralDiversity",
                        y.label = "Mean floral diversity",
                        x.label = "Latitude (log)",
                        parasite = FALSE,
                        scale.y = TRUE,
                        angle.x = TRUE,
                        x.axis.lab = FALSE)



################################################################################
## Bee abundance and lat
################################################################################
p3 <- plot_cond_effects(fit.bombus.l, data = spec.uni,
                        this.response = "scaleNetBombusAbundance",
                        this.effect = "Lat",
                        significance = "97",
                        dat.x = "Lat",
                        dat.y = "Net_BombusAbundance",
                        y.label = expression(atop(italic("Bombus"), "abundance (log)")),
                        x.label = "Latitude (log)",
                        parasite = FALSE,
                        scale.y = TRUE,
                        angle.x = TRUE)




## Honeybee abundance
p4 <- plot_cond_effects(fit.bombus.l, data = spec.uni,
                        this.response = "scaleNetHBAbundance",
                        this.effect = "Lat",
                        significance = "97",
                        dat.x = "Lat",
                        dat.y = "Net_HBAbundance",
                        y.label = expression(atop(italic("Apis"), "abundance (log)")),
                        x.label = "Latitude (log)",
                        parasite = FALSE,
                        scale.y = TRUE,
                        angle.x = TRUE
                        )




# lat_community <- ggarrange(p2,p1,p3,p4, #plots that are going to be included in this multipanel figure
#                        labels = c("A", "B", "C","D"), #labels given each panel 
#                        ncol = 2, nrow = 2, #adjust plot space 
#                        font.label = list(color = "white"))

# ggsave(lat_community, file="figures/lat_community.pdf",
#        height=8, width=12)

## ***********************************************************************
## lat and crithidia
## ***********************************************************************
p5 <- plot_cond_effects(fit.bombus.l, data = spec.uni,
                        this.response = "CrithidiaPresence",
                        this.effect = "Lat",
                        significance = "ns",
                        dat.x = "Lat",
                        dat.y = "CrithidiaParasitismRate",
                        y.label = expression(atop(italic("Crithidia")~ "prevalence", "in"~italic("Bombus"))),
                        x.label = "Latitude (log)",
                        angle.x = TRUE,
                        x.axis.lab = FALSE)

## ***********************************************************************
## lat and apicystis
## ***********************************************************************
p6 <- plot_cond_effects(fit.bombus.l, data = spec.uni,
                        this.response = "ApicystisSpp",
                        this.effect = "Lat",
                        significance = "97",
                        dat.x = "Lat",
                        dat.y = "ApicystisParasitismRate",
                        y.label = expression(atop(italic("Apicystis")~ "prevalence", "in"~italic("Bombus"))),
                        x.label = "Latitude (log)",
                        angle.x = TRUE,
                        x.axis.lab = FALSE)

################################################################################
## crithida ~ lat in Apis
## ***************************************************************************
# Load model for lat
load(file="saved/parasiteFit_Apis_CrithidiaPresenceApicystisSpp_lat_lat.Rdata")
fit.apis.l <- fit.parasite.apis

p7 <- plot_cond_effects(fit.apis.l, data = spec.uni,
                        this.response = "CrithidiaPresence",
                        this.effect = "Lat",
                        significance = "ns",
                        dat.x = "Lat",
                        dat.y = "CrithidiaParasitismRate",
                        y.label = expression(atop(italic("Crithidia")~ "prevalence", "in"~italic("Apis"))),
                        x.label = "Latitude (log)",
                        angle.x = TRUE,
                        x.axis.lab = FALSE)



################################################################################
## Lat and apicystis in apis
################################################################################
p8 <- plot_cond_effects(fit.apis.l, data = spec.uni,
                        this.response = "ApicystisSpp",
                        this.effect = "Lat",
                        significance = "97",
                        dat.x = "Lat",
                        dat.y = "ApicystisParasitismRate",
                        y.label = expression(atop(italic("Apicystis")~ "prevalence", "in"~italic("Apis"))),
                        x.label = "Latitude (log)",
                        angle.x = TRUE,
                        x.axis.lab = FALSE)

lat_full <- ggarrange(p5, p7, p6, p8,
                         p2,p1,p3,p4,
                         labels = c("A", "B", "C","D", "E", "F", "G", "H"),
                         font.label = list(size = 12),
                         common.legend = TRUE,
                         legend = "bottom",
                         ncol = 2, nrow = 4, 
                         heights = c(3, 3, 3, 3), widths = 1)

ggsave(lat_full, file="figures/fig4_lat_full.pdf",
       height=22, width=18, units = "cm")

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

p9 <- plot_cond_effects(fit.fd, data = spec.uni,
                        this.response = "scaleNetBeeDiversity",
                        this.effect = "MeanFloralDiversity",
                        significance = "97",
                        dat.x = "MeanFloralDiversity",
                        dat.y = "Net_BeeDiversity",
                        y.label = "Bee diversity",
                        x.label = "Mean Floral Diversity",
                        parasite = FALSE,
                        scale.y = TRUE,
                        angle.x = TRUE)

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

ggsave(plantdiv_beediv, file="figures/figS3_Beediv_floraldiv.pdf",
       height=4, width=5)

