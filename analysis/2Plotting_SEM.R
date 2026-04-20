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
source("src/plotMod.R")


site.screened <- spec.net %>%
  group_by(Site, Year, SampleRound) %>%
  summarise(SiteScreened = sum(!is.na(Apidae)))

spec.net <- left_join(spec.net, site.screened)

## data subsetted to unique values
spec.uni <- spec.net[spec.net$Weights ==1,]


## ***************************************************************************
# Load model for bee diversity
load(file="saved/parasiteFit_Bombus_CrithidiaPresenceApicystisSpp_bee_div_cp.Rdata")
fit.bombus.bd <- fit.parasite.bombus

## ***************************************************************************
## Crithidia ~ bee diversity
p1 <- plot_cond_effects(fit.bombus.bd, data = spec.uni,
                        significance = "97",
                        x.axis.lab = FALSE,
                        y.label = expression(bolditalic("Crithidia") ~ bold("prevalence")),
                        x.label = "Bee diversity",
                        text.size = 14)




################################################################################
## Apicystis ~ bee diversity Bombus
################################################################################

p2 <- plot_cond_effects(fit.bombus.bd, data = spec.uni,
                        this.response = "ApicystisSpp",
                        dat.y = "ApicystisParasitismRate",
                        significance = "97",
                        y.label = expression(bolditalic("Apicystis") ~ bold("prevalence")),
                        x.label = "Bee diversity",
                        text.size = 14)


## ***************************************************************************
# Load model for floral diversity
load(file="saved/parasiteFit_bombus_CrithidiaPresenceApicystisSpp_floral_div_lat.Rdata")
fit.bombus.fd <- fit.parasite.bombus

## ***************************************************************************
## Crithidia ~ floral diversity

p3 <- plot_cond_effects(fit.bombus.fd, data = spec.uni,
                        this.effect = "MeanFloralDiversity",
                        dat.x = "MeanFloralDiversity",
                        significance = "97",
                        x.axis.lab = FALSE,
                        y.label = expression(bolditalic("Crithidia") ~ bold("prevalence")),
                        x.label = "Mean Floral diversity",
                        text.size = 14)




## ***************************************************************************
## Apicystis ~ floral diversity

p4 <- plot_cond_effects(fit.bombus.fd, data = spec.uni,
                        this.response = "ApicystisSpp",
                        this.effect = "MeanFloralDiversity",
                        dat.x = "MeanFloralDiversity", 
                        dat.y = "ApicystisParasitismRate",
                        significance = "97",
                        y.label = expression(bolditalic("Apicystis") ~ bold("prevalence")),
                        x.label = "Mean Floral diversity",
                        text.size = 14)



    
parasite.dilution <- ggarrange(p1, p3, p2, p4,
                            labels = c("A", "B", "C","D"),
                            font.label = list(size = 12),
                            ncol = 2, nrow = 2,
                            common.legend = TRUE,
                            legend = "bottom")

ggsave(parasite.dilution, file="figures/fig2_parasite_diversity.pdf", 
       height=6, width=10)

## ***************************************************************************
# Load model for bombus abundance
load(file="saved/parasiteFit_Bombus_CrithidiaPresenceApicystisSpp_bombus_abund_cp.Rdata")
fit.bombus.ba <- fit.parasite.bombus

## ***************************************************************************
## crithidia ~ bombus abundance

p5 <- plot_cond_effects(fit.bombus.ba, data = spec.uni,
                        this.effect = "Net_BombusAbundance",
                        dat.x = "Net_BombusAbundance",
                        x.axis.lab = FALSE,
                        significance = "ns",
                        y.label = expression(bolditalic("Crithidia") ~ bold("prevalence")),
                        x.label = expression(bolditalic("Bombus") ~ bold("abundance")),
                        text.size = 14)


################################################################################
## Apicysits ~ bombus abundance
################################################################################
p6 <- plot_cond_effects(fit.bombus.ba, data = spec.uni,
                        this.response = "ApicystisSpp",
                        this.effect = "Net_BombusAbundance",
                        dat.x = "Net_BombusAbundance", 
                        dat.y = "ApicystisParasitismRate",
                        significance = "97",
                        y.label = expression(bolditalic("Apicystis") ~ bold("prevalence")),
                        x.label = expression(bolditalic("Bombus") ~ bold("abundance")),
                        text.size = 14)



## ***************************************************************************
# Load model for apis abundance
load(file="saved/parasiteFit_Bombus_CrithidiaPresenceApicystisSpp_hb_abund_cp.Rdata")
fit.bombus.ha <- fit.parasite.bombus

################################################################################
## Crithidia ~ apis abundance
################################################################################

p7 <- plot_cond_effects(fit.bombus.ha, data = spec.uni,
                        this.effect = "Net_HBAbundance",
                        dat.x = "Net_HBAbundance",
                        significance = "97",
                        x.axis.lab = FALSE,
                        y.label = expression(bolditalic("Crithidia") ~ bold("prevalence")),
                        x.label = expression(bolditalic("Apis") ~ bold("abundance")),
                        text.size = 14)


################################################################################
## Apicysits ~ apis abundance
################################################################################

p8 <- plot_cond_effects(fit.bombus.ha, data = spec.uni,
                        this.response = "ApicystisSpp",
                        this.effect = "Net_HBAbundance",
                        dat.x = "Net_HBAbundance", 
                        dat.y = "ApicystisParasitismRate",
                        significance = "97",
                        y.label = expression(bolditalic("Apicystis") ~ bold("prevalence")),
                        x.label = expression(bolditalic("Apis") ~ bold("abundance")),
                        text.size = 14)


parasite.amplification <- ggarrange(p5, p7, p6, p8,
                                    nrow= 2, ncol = 2,
                                    font.label = list(size = 12),
                                    labels = c("A", "B", "C", "D"),
                                    common.legend = TRUE,
                                    legend = "bottom")

ggsave(parasite.amplification, file="figures/fig3_parasite_amplification.pdf",
       height=6, width=10)


