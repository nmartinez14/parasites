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
source("src/runPlotFreqModelDiagnostics.R")

## site or lat as the geographic variable
site.or.lat <- "lat"

## all of the variables that are explanatory variables and thus need
## to be centered
vars_yearsr <- c("MeanFloralDiversity",
                 "Net_BeeDiversity",
                 "Net_BeeAbundance",
                 "Net_BombusAbundance",
                 "Net_HBAbundance",
                 "SRDoyPoly1",
                 "SRDoyPoly2"
                 )
vars_yearsrsp <- "rare.degree"
vars_site <- c("Lat", "Area")

variables.to.log <- c("rare.degree", "Lat", "Net_BeeAbundance", "Area")

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
## raw, non standardized data for plotting
spec.orig <- prepDataSEM(spec.net, variables.to.log, variables.to.log.1,
                         standardize=FALSE)

## Make SEM weights and standardize data.
spec.net <- prepDataSEM(spec.net, variables.to.log, variables.to.log.1,
                        vars_yearsr = vars_yearsr,
                        vars_yearsrsp = vars_yearsrsp,
                        vars_site=vars_site)

spec.net$Site <- as.character(spec.net$Site)
## otherwise levels with no data are not properly dropped using subset
spec.net$Year <- as.character(spec.net$Year)
spec.net$GenusSpecies <- as.character(spec.net$GenusSpecies)


## bombus only data
spec.bombus <- makeGenusSubset(spec.net, "Bombus")
## apis only data
spec.apis <- makeGenusSubset(spec.net, "Apis")
## can repeat for other genera but have deleted since the models do
## not converge

## define all the formulas for the different parts of the models
source("src/plant_poll_models.R")

save(spec.net, spec.orig, file="saved/spec_weights.Rdata")

## Species that are not in the phylogeny are not used. brms is not
## allowing an incomplete phylogeny, to avoid the error we changed the
## species not present to one that is in the phylogeny.  We chose a
## species for which we did not do parasite screening and should not
## influence results.
not_in_phylo <- unique(spec.net$GenusSpecies[!spec.net$GenusSpecies
                                             %in%
                                             phylo$tip.label])

## only bombus model is muti species given the others do not converge
spec.bombus$GenusSpecies[spec.bombus$GenusSpecies %in%
                         not_in_phylo]<- "Agapostemon angelicus"
## **********************************************************
## community model, check assumptions first before adding parasites
## **********************************************************

run_plot_freq_model_diagnostics(remove_subset_formula(formula.flower.div),
                                this_data=spec.net[spec.net$Weights == 1,],
                                this_family="students", site.lat=site.or.lat)

run_plot_freq_model_diagnostics(remove_subset_formula(formula.bombus.abund),
                                this_data=spec.net[spec.net$Weights == 1,],
                                this_family="students", site.lat=site.or.lat)

run_plot_freq_model_diagnostics(remove_subset_formula(formula.HB.abund),
                                this_data=spec.net[spec.net$Weights == 1,],
                                this_family="students", site.lat=site.or.lat)

run_plot_freq_model_diagnostics(remove_subset_formula(formula.bee.div),
                                this_data=spec.net[spec.net$Weights == 1,],
                                this_family="gaussian",
                                site.lat=site.or.lat)

## **********************************************************
## Parasite models set up
## **********************************************************
## Bombus predictor variables
## phylogeny must be last in all xvar sets

xvars.fd <-  c("MeanFloralDiversity",
               "Year",
               "SRDoy",
               "Lat","Area",
               "(1|Site)",
               "(1|gr(GenusSpecies, cov = phylo_matrix))")

xvars.bd <-  c("Net_BeeDiversity",
               "Year",
               "SRDoy",
               "Lat","Area",
               "(1|Site)",
               "(1|gr(GenusSpecies, cov = phylo_matrix))")

xvars.ba <-  c("Net_BombusAbundance",
               "SRDoy", "Year",
               "Lat","Area",
               "(1|Site)",
               "(1|gr(GenusSpecies, cov = phylo_matrix))")


xvars.ha <-  c("Net_HBAbundance",
               "Year",
               "SRDoy",
               "Lat","Area",
               "(1|Site)",
               "(1|gr(GenusSpecies, cov = phylo_matrix))")

xvars.d <-  c("rare.degree",
               "Year",
               "SRDoy",
               "Lat","Area",
               "(1|Site)",
               "(1|gr(GenusSpecies, cov = phylo_matrix))")

xvars.l <-  c("Lat",
               "(1|Site)",
               "(1|gr(GenusSpecies, cov = phylo_matrix))")

xvars.a <-  c("Area",
               "(1|Site)",
               "(1|gr(GenusSpecies, cov = phylo_matrix))")

xvars.sr <-  c("SRDoy",
               "(1|Site)",
               "(1|gr(GenusSpecies, cov = phylo_matrix))")

xvars.y <-  c("Year",
               "(1|Site)",
               "(1|gr(GenusSpecies, cov = phylo_matrix))")


## **********************************************************
## Parasite presence
## **********************************************************
## Bombus
## **********************************************************

bombus.fd <- runCombinedParasiteModels(spec.data= spec.bombus,
                                       species.group="bombus",
                                       xvars=xvars.fd,
                                       ncores=ncores,
                                       data2= list(phylo_matrix=phylo_matrix),
                                       xvar.name="floral_div")

bombus.bd <- runCombinedParasiteModels(spec.data= spec.bombus,
                                       species.group="bombus",
                                       xvars=xvars.bd,
                                       ncores=ncores,
                                       data2= list(phylo_matrix=phylo_matrix),
                                       xvar.name="bee_div")


bombus.ba <- runCombinedParasiteModels(spec.data= spec.bombus,
                                         species.group="bombus",
                                         xvars=xvars.ba,
                                         ncores=ncores,
                                         data2= list(phylo_matrix=phylo_matrix),
                                         xvar.name="bombus_abund")

bombus.ha <- runCombinedParasiteModels(spec.data= spec.bombus,
                                       species.group="bombus",
                                       xvars=xvars.ha,
                                       ncores=ncores,
                                       data2= list(phylo_matrix=phylo_matrix),
                                       xvar.name="hb_abund")

bombus.d <- runCombinedParasiteModels(spec.data= spec.bombus,
                                       species.group="bombus",
                                       xvars=xvars.d,
                                       ncores=ncores,
                                       data2= list(phylo_matrix=phylo_matrix),
                                       xvar.name="degree")

bombus.l <- runCombinedParasiteModels(spec.data= spec.bombus,
                                       species.group="bombus",
                                       xvars=xvars.l,
                                       ncores=ncores,
                                       data2= list(phylo_matrix=phylo_matrix),
                                       xvar.name="lat")

bombus.a <- runCombinedParasiteModels(spec.data= spec.bombus,
                                      species.group="bombus",
                                      xvars=xvars.a,
                                      ncores=ncores,
                                      data2= list(phylo_matrix=phylo_matrix),
                                      xvar.name="area")

bombus.sr <- runCombinedParasiteModels(spec.data= spec.bombus,
                                      species.group="bombus",
                                      xvars=xvars.sr,
                                      ncores=ncores,
                                      data2= list(phylo_matrix=phylo_matrix),
                                      xvar.name="SRDOY")

bombus.y <- runCombinedParasiteModels(spec.data= spec.bombus,
                                      species.group="bombus",
                                      xvars=xvars.y,
                                      ncores=ncores,
                                      data2= list(phylo_matrix=phylo_matrix),
                                      xvar.name="year")


## **********************************************************
## Parasite presence
## **********************************************************
## Apis
## **********************************************************

apis.fd <- runCombinedParasiteModels(spec.data= spec.apis,
                                     species.group="apis",
                                     xvars=xvars.fd[-length(xvars.fd)],
                                     ncores=ncores,
                                     xvar.name= "floral_div")


apis.bd <- runCombinedParasiteModels(spec.data= spec.apis,
                                     species.group="apis",
                                     xvars=xvars.bd[-length(xvars.bd)],
                                     ncores=ncores,
                                     xvar.name= "bee_div")

apis.ba <- runCombinedParasiteModels(spec.data= spec.apis,
                                     species.group="apis",
                                     xvars=xvars.ba[-length(xvars.ba)],
                                     ncores=ncores,
                                     xvar.name= "bombus_abund")

apis.ha <- runCombinedParasiteModels(spec.data= spec.apis,
                                     species.group="apis",
                                     xvars=xvars.ha[-length(xvars.ha)],
                                     ncores=ncores,
                                     xvar.name= "ha_abund")

apis.d <- runCombinedParasiteModels(spec.data= spec.apis,
                                     species.group="apis",
                                     xvars=xvars.d[-length(xvars.d)],
                                     ncores=ncores,
                                     xvar.name= "degree")

apis.l <- runCombinedParasiteModels(spec.data= spec.apis,
                                     species.group="apis",
                                     xvars=xvars.l[-length(xvars.l)],
                                     ncores=ncores,
                                     xvar.name= "lat")

apis.a <- runCombinedParasiteModels(spec.data= spec.apis,
                                     species.group="apis",
                                     xvars=xvars.a[-length(xvars.a)],
                                     ncores=ncores,
                                     xvar.name= "area")

apis.sr <- runCombinedParasiteModels(spec.data= spec.apis,
                                     species.group="apis",
                                     xvars=xvars.sr[-length(xvars.sr)],
                                     ncores=ncores,
                                     xvar.name= "SRDOY")

apis.y <- runCombinedParasiteModels(spec.data= spec.apis,
                                     species.group="apis",
                                     xvars=xvars.y[-length(xvars.y)],
                                     ncores=ncores,
                                     xvar.name= "year")

