rm(list=ls())
source("lab_paths.R")
setwd(local.path)

ncores <- 6

setwd("parasites/analysis")
source("src/misc.R")
source("src/writeResultsTable.R")
source("src/makeMultiLevelData.R")
source("src/runParasiteModels.R")
source("src/standardize_weights.R")
source("src/pp_checks.R")

## "site" for community models as predictor variable or 
## "lat" to have site as random effect
site.or.lat <- "lat"

variables.to.log <- c("rare.degree", "Lat", "Area")

## some zeros in data
variables.to.log.1 <- c("Net_HBAbundance", "Net_BombusAbundance")

## loads specimen data
source("src/init.R")

## drop VC (Valles caldera) because it was more of a grassland, we
## only surveyed it one year)
print("Before dropping VC")
dim(spec.net)
spec.net <- filter(spec.net, Site != "VC")
print("After dropping VC")
dim(spec.net)
## because only Bombus and apis models converge, set the rest of
## the trait data to NA so that the variables scale properly
screened.bombus <- unique(spec.net$GenusSpecies[spec.net$Apidae == 1 &
                                                spec.net$Genus == "Bombus"])
screened.bombus <- screened.bombus[!is.na(screened.bombus)]


spec.net$rare.degree[!spec.net$GenusSpecies %in%
                     c("Apis mellifera", screened.bombus)] <- NA

dim(spec.net)


## Make SEM weights
spec.net <- prepDataSEM(spec.net, variables.to.log, variables.to.log.1,
                        standardize = FALSE)

spec.net$Site <- as.character(spec.net$Site)
## otherwise levels with no data are not properly dropped using subset
spec.net$Year <- as.character(spec.net$Year)
spec.net$GenusSpecies <- as.character(spec.net$GenusSpecies)

spec.net$BombusWeights <- ifelse(spec.net$Apidae == 1 & spec.net$Genus == "Bombus", 1, 0)
spec.net$ApisWeights <- ifelse(spec.net$Apidae == 1 & spec.net$Genus == "Apis", 1, 0)

## define all the formulas for the different parts of the models
source("src/plant_poll_models.R")

save(spec.net, file="saved/spec_weights.Rdata")

## Species that are not in the phylogeny are not used. brms is not
## allowing an incomplete phylogeny, to avoid the error we changed the
## species not present to one that is in the phylogeny.  We chose a
## species for which we did not do parasite screening and should not
## influence results.
not_in_phylo <- unique(spec.net$GenusSpecies[!spec.net$GenusSpecies
                                             %in%
                                             phylo$tip.label])

## only bombus model is multi species given the others do not converge
spec.net$GenusSpecies[spec.net$GenusSpecies %in%
                         not_in_phylo]<- "Agapostemon angelicus"


## **********************************************************
## Parasite models set up
## **********************************************************
## Bombus predictor variables
## phylogeny must be last in all xvar sets

xvars.fd <-  c("scale(MeanFloralDiversity)",
               "scale(Cumulative_Precip)",
               "scale(Lat)","scale(Area)",
               "(1|Site)",
               "(1|gr(GenusSpecies, cov = phylo_matrix))")

xvars.bd <-  c("scale(Net_BeeDiversity)",
               "scale(Cumulative_Precip)",
               "scale(Lat)","scale(Area)",
               "(1|Site)",
               "(1|gr(GenusSpecies, cov = phylo_matrix))")

xvars.ba <-  c("scale(Net_BombusAbundance)",
               "scale(MeanFloralDiversity)",
               "scale(Cumulative_Precip)",
               "scale(Lat)","scale(Area)",
               "(1|Site)",
               "(1|gr(GenusSpecies, cov = phylo_matrix))")


xvars.ha <-  c("scale(Net_HBAbundance)",
               "scale(MeanFloralDiversity)",
               "scale(Cumulative_Precip)",
               "scale(Lat)","scale(Area)",
               "(1|Site)",
               "(1|gr(GenusSpecies, cov = phylo_matrix))")

xvars.d <-  c("scale(rare.degree)",
              "scale(Cumulative_Precip)",
              "scale(Lat)","scale(Area)",
               "(1|Site)",
               "(1|gr(GenusSpecies, cov = phylo_matrix))")

xvars.l <-  c("scale(Lat)",
               "(1|Site)",
               "(1|gr(GenusSpecies, cov = phylo_matrix))")

xvars.a <-  c("scale(Area)",
               "(1|Site)",
               "(1|gr(GenusSpecies, cov = phylo_matrix))")


xvars.cp <-  c("scale(Cumulative_Precip)",
               "scale(Lat)",
               "(1|Site)",
               "(1|gr(GenusSpecies, cov = phylo_matrix))")


## **********************************************************
## Parasite presence
## **********************************************************
## Bombus
## **********************************************************
## floral diversity
mod.fd <- runCombinedParasiteModels(spec.data= spec.net,
                                       xvars=xvars.fd,
                                       ncores=ncores,
                                       data2= list(phylo_matrix=phylo_matrix),
                                       xvar.name="floral_div",
                                       top.level = "lat")


mod.fd <- runCombinedParasiteModels(spec.data= spec.net,
                                       xvars=xvars.fd,
                                       ncores=ncores,
                                       data2= list(phylo_matrix=phylo_matrix),
                                       xvar.name="floral_div",
                                       top.level = "cp")


mod.fd <- runCombinedParasiteModels(spec.data= spec.net,
                                       xvars=xvars.fd,
                                       ncores=ncores,
                                       data2= list(phylo_matrix=phylo_matrix),
                                       xvar.name="floral_div",
                                       top.level = "area")
## bee diversity

mod.bd <- runCombinedParasiteModels(spec.data= spec.net,
                                       xvars=xvars.bd,
                                       ncores=ncores,
                                       data2= list(phylo_matrix=phylo_matrix),
                                       xvar.name="bee_div",
                                       top.level = "lat")


mod.bd <- runCombinedParasiteModels(spec.data= spec.net,
                                       xvars=xvars.bd,
                                       ncores=ncores,
                                       data2= list(phylo_matrix=phylo_matrix),
                                       xvar.name="bee_div",
                                       top.level = "cp")

mod.bd <- runCombinedParasiteModels(spec.data= spec.net,
                                       xvars=xvars.bd,
                                       ncores=ncores,
                                       data2= list(phylo_matrix=phylo_matrix),
                                       xvar.name="bee_div",
                                       top.level = "area")

mod.bd <- runCombinedParasiteModels(spec.data= spec.net,
                                       xvars=xvars.bd,
                                       ncores=ncores,
                                       data2= list(phylo_matrix=phylo_matrix),
                                       xvar.name="bee_div",
                                       top.level = "fd")

## bombus abundance
mod.ba <- runCombinedParasiteModels(spec.data= spec.net,
                                         xvars=xvars.ba,
                                         ncores=ncores,
                                         data2= list(phylo_matrix=phylo_matrix),
                                         xvar.name="bombus_abund",
                                         top.level = "lat")

mod.ba <- runCombinedParasiteModels(spec.data= spec.net,
                                    xvars=xvars.ba,
                                    ncores=ncores,
                                    data2= list(phylo_matrix=phylo_matrix),
                                    xvar.name="bombus_abund",
                                    top.level = "area")

mod.ba <- runCombinedParasiteModels(spec.data= spec.net,
                                    xvars=xvars.ba,
                                    ncores=ncores,
                                    data2= list(phylo_matrix=phylo_matrix),
                                    xvar.name="bombus_abund",
                                    top.level = "cp")

mod.ba <- runCombinedParasiteModels(spec.data= spec.net,
                                    xvars=xvars.ba,
                                    ncores=ncores,
                                    data2= list(phylo_matrix=phylo_matrix),
                                    xvar.name="bombus_abund",
                                    top.level = "fd")
## honey bee abundance

mod.ha <- runCombinedParasiteModels(spec.data= spec.net,
                                       xvars=xvars.ha,
                                       ncores=ncores,
                                       data2= list(phylo_matrix=phylo_matrix),
                                       xvar.name="hb_abund",
                                       top.level = "lat")

mod.ha <- runCombinedParasiteModels(spec.data= spec.net,
                                    xvars=xvars.ha,
                                    ncores=ncores,
                                    data2= list(phylo_matrix=phylo_matrix),
                                    xvar.name="hb_abund",
                                    top.level = "area")

mod.ha <- runCombinedParasiteModels(spec.data= spec.net,
                                    xvars=xvars.ha,
                                    ncores=ncores,
                                    data2= list(phylo_matrix=phylo_matrix),
                                    xvar.name="hb_abund",
                                    top.level = "cp")

mod.ha <- runCombinedParasiteModels(spec.data= spec.net,
                                    xvars=xvars.ha,
                                    ncores=ncores,
                                    data2= list(phylo_matrix=phylo_matrix),
                                    xvar.name="hb_abund",
                                    top.level = "fd")
## diet breadth

mod.d <- runCombinedParasiteModels(spec.data= spec.net,
                                       xvars=xvars.d,
                                       ncores=ncores,
                                       data2= list(phylo_matrix=phylo_matrix),
                                       xvar.name="degree",
                                       top.level = "lat")

mod.d <- runCombinedParasiteModels(spec.data= spec.net,
                                      xvars=xvars.d,
                                      ncores=ncores,
                                      data2= list(phylo_matrix=phylo_matrix),
                                      xvar.name="degree",
                                      top.level = "area")

mod.d <- runCombinedParasiteModels(spec.data= spec.net,
                                      xvars=xvars.d,
                                      ncores=ncores,
                                      data2= list(phylo_matrix=phylo_matrix),
                                      xvar.name="degree",
                                      top.level = "cp")

mod.d <- runCombinedParasiteModels(spec.data= spec.net,
                                      xvars=xvars.d,
                                      ncores=ncores,
                                      data2= list(phylo_matrix=phylo_matrix),
                                      xvar.name="degree",
                                      top.level = "fd")

## Latitude

mod.l <- runCombinedParasiteModels(spec.data= spec.net,
                                       xvars=xvars.l,
                                       ncores=ncores,
                                       data2= list(phylo_matrix=phylo_matrix),
                                       xvar.name="lat",
                                      top.level = "lat")
## Area

mod.a <- runCombinedParasiteModels(spec.data= spec.net,
                                      xvars=xvars.a,
                                      ncores=ncores,
                                      data2= list(phylo_matrix=phylo_matrix),
                                      xvar.name="area",
                                      top.level = "area")
## Cumulative precipitation

mod.cp <- runCombinedParasiteModels(spec.data= spec.net,
                                      xvars=xvars.cp,
                                      ncores=ncores,
                                      data2= list(phylo_matrix=phylo_matrix),
                                      xvar.name="cp",
                                      top.level = "cp")



