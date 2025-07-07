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
## Parasite models set up
## **********************************************************
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

xvars.ba.2 <-  c("Net_BombusAbundance",
               "SRDoy", 
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
## Parasite presence
## **********************************************************
## Bombus
## **********************************************************
table(spec.bombus$CrithidiaPresence[spec.bombus$WeightsPar == 1])
table(spec.bombus$ApicystisSpp[spec.bombus$WeightsPar ==1 ])

## model running function also runs frequentist models and
## diagnostics. Bayesian models are beta binomial, frequentist are
## binomial because there are no packages that include beta binomial


bombus.fd <- runCombinedParasiteModels(spec.data= spec.bombus,
                                       species.group="bombus",
                                       xvars=xvars.fd,
                                       ncores=ncores,
                                       data2= list(phylo_matrix=phylo_matrix),
                                       xvar.name="floral_div", 
                                       site.lat = "lat")

bombus.bd <- runCombinedParasiteModels(spec.data= spec.bombus,
                                       species.group="bombus",
                                       xvars=xvars.bd,
                                       ncores=ncores,
                                       data2= list(phylo_matrix=phylo_matrix),
                                       xvar.name="bee_div", 
                                       site.lat = "lat")

## For Critihidia: bombus abundance is collinear with year. I've ran individual 
## models for each parasite to account for this and remove year in the Crithidia.

bombus.ba.a <- runCombinedParasiteModels(spec.data= spec.bombus,
                                         species.group="bombus",
                                         parasites= "ApicystisSpp",
                                         xvars=xvars.ba,
                                         ncores=ncores,
                                         data2= list(phylo_matrix=phylo_matrix),
                                         xvar.name="bombus_abund", 
                                         site.lat = "lat")

bombus.ba.c <- runCombinedParasiteModels(spec.data= spec.bombus,
                                       species.group="bombus",
                                       parasites = "CrithidiaPresence",
                                       xvars=xvars.ba.2,
                                       ncores=ncores,
                                       data2= list(phylo_matrix=phylo_matrix),
                                       xvar.name="bombus_abund", 
                                       site.lat = "lat")


bombus.ha <- runCombinedParasiteModels(spec.data= spec.bombus,
                                       species.group="bombus",
                                       xvars=xvars.ha,
                                       ncores=ncores,
                                       data2= list(phylo_matrix=phylo_matrix),
                                       xvar.name="hb_abund", 
                                       site.lat = "lat")

bombus.d <- runCombinedParasiteModels(spec.data= spec.bombus,
                                       species.group="bombus",
                                       xvars=xvars.d,
                                       ncores=ncores,
                                       data2= list(phylo_matrix=phylo_matrix),
                                       xvar.name="degree", 
                                      site.lat = "lat")

bombus.l <- runCombinedParasiteModels(spec.data= spec.bombus,
                                       species.group="bombus",
                                       xvars=xvars.l,
                                       ncores=ncores,
                                       data2= list(phylo_matrix=phylo_matrix),
                                       xvar.name="lat", 
                                      site.lat = "lat")

bombus.a <- runCombinedParasiteModels(spec.data= spec.bombus,
                                      species.group="bombus",
                                      xvars=xvars.a,
                                      ncores=ncores,
                                      data2= list(phylo_matrix=phylo_matrix),
                                      xvar.name="area", 
                                      site.lat = "lat")

bombus.sr <- runCombinedParasiteModels(spec.data= spec.bombus,
                                      species.group="bombus",
                                      xvars=xvars.sr,
                                      ncores=ncores,
                                      data2= list(phylo_matrix=phylo_matrix),
                                      xvar.name="SRDOY", 
                                      site.lat = "lat")

bombus.y <- runCombinedParasiteModels(spec.data= spec.bombus,
                                      species.group="bombus",
                                      xvars=xvars.y,
                                      ncores=ncores,
                                      data2= list(phylo_matrix=phylo_matrix),
                                      xvar.name="year", 
                                      site.lat = "lat")

## **********************************************************
## Bombus loo summaries
## **********************************************************
## crithidia
## **********************************************************
## not including ss because HB abundance and bombus abundance are colinear

load("saved/parasiteFit_bombus_CrithidiaPresenceApicystisSpp_ba_noSRDoy.Rdata")
bombus.ba <- fit.parasite
load("saved/parasiteFit_bombus_CrithidiaPresenceApicystisSpp_apis_abundance_noSRDoy.Rdata")
bombus.ha <- fit.parasite
load("saved/parasiteFit_bombus_CrithidiaPresenceApicystisSpp_lat_diversity.Rdata")
bombus.div <- fit.parasite

loo.crithidia.ba <- loo(bombus.ba, resp="CrithidiaPresence")
loo.crithidia.ha <- loo(bombus.ha, resp="CrithidiaPresence")
loo.crithidia.div <- loo(bombus.div, resp="CrithidiaPresence")
loo_compare(loo.crithidia.ba, loo.crithidia.ha, loo.crithidia.div)



## **********************************************************
## apicystis
## **********************************************************

load("saved/parasiteFit_bombus_CrithidiaPresenceApicystisSpp_lat_social_species.Rdata")
bombus.ss <- fit.parasite
load("saved/parasiteFit_bombus_CrithidiaPresenceApicystisSpp_bombus_abundance.Rdata")
bombus.ba <- fit.parasite
load("saved/parasiteFit_bombus_CrithidiaPresenceApicystisSpp_apis_abundance.Rdata")
bombus.ha <- fit.parasite

loo.apicystis.ba <- loo(bombus.ba, resp="ApicystisSpp")
loo.apicystis.ss <- loo(bombus.ss, resp="ApicystisSpp")
loo.apicystis.ha <- loo(bombus.ha, resp="ApicystisSpp")

loo_compare(loo.apicystis.ss,loo.apicystis.ba, loo.apicystis.ha)
## The best fit is the abundance of bombus and apis together, which in
## this model are not colinear.

## **********************************************************
## Honey bees
## **********************************************************
table(spec.apis$CrithidiaPresence[spec.apis$WeightsPar == 1])
table(spec.apis$ApicystisSpp[spec.apis$WeightsPar ==1 ])

xvars.fd <-  c("MeanFloralDiversity",
               "Year",
               "SRDoy",
               "Lat","Area",
               "(1|Site)")

xvars.fd.2 <-  c("MeanFloralDiversity",
               "Year",
               "Area",
               "(1|Site)")

xvars.bd <-  c("Net_BeeDiversity",
               "Area",
               "(1|Site)")

xvars.bd.2 <-  c("Net_BeeDiversity",
               "Year",
               "SRDoy",
               "Area",
               "(1|Site)")

xvars.ba <-  c("Net_BombusAbundance",
               "Year",
               "SRDoy",
               "Area",
               "(1|Site)")

xvars.ha <-  c("Net_HBAbundance",
               "Year",
               "SRDoy",
               "Lat","Area",
               "(1|Site)")

xvars.ha.2 <-  c("Net_HBAbundance",
               "Year",
               "SRDoy",
               "Area",
               "(1|Site)")

xvars.d <-  c("rare.degree",
              "Year",
              "SRDoy",
              "Lat","Area",
              "(1|Site)")

xvars.d.2 <-  c("rare.degree",
              "SRDoy","Area",
              "(1|Site)")

xvars.l <-  c("Lat",
              "(1|Site)")

xvars.a <-  c("Area",
              "(1|Site)")

xvars.sr <-  c("SRDoy",
               "(1|Site)")

xvars.y <-  c("Year",
              "(1|Site)")

## Ran split the parasite models based on appropriate noncolinear variables

apis.fd.a <- runCombinedParasiteModels(spec.data= spec.apis,
                                     species.group="apis",
                                     parasites = "ApicystisSpp",
                                     xvars=xvars.fd,
                                     ncores=ncores,
                                     site.lat = site.or.lat,
                                     xvar.name= "floral_div")

apis.fd.c <- runCombinedParasiteModels(spec.data= spec.apis,
                                     species.group="apis",
                                     parasites = "CrithidiaPresence",
                                     xvars=xvars.fd.2,
                                     ncores=ncores,
                                     site.lat = site.or.lat,
                                     xvar.name= "floral_div")

apis.bd.a <- runCombinedParasiteModels(spec.data= spec.apis,
                                     species.group="apis",
                                     parasites= "ApicystisSpp",
                                     xvars=xvars.bd,
                                     ncores=ncores,
                                     site.lat = site.or.lat,
                                     xvar.name= "bee_div")

apis.bd.c <- runCombinedParasiteModels(spec.data= spec.apis,
                                     species.group="apis",
                                     parasites = "CrithidiaPresence",
                                     xvars=xvars.bd.2,
                                     ncores=ncores,
                                     site.lat = site.or.lat,
                                     xvar.name= "bee_div")

apis.ba <- runCombinedParasiteModels(spec.data= spec.apis,
                                     species.group="apis",
                                     xvars=xvars.ba,
                                     ncores=ncores,
                                     site.lat = site.or.lat,
                                     xvar.name= "bombus_abund")

apis.ha.a <- runCombinedParasiteModels(spec.data= spec.apis,
                                     species.group="apis",
                                     parasites = "ApicystisSpp",
                                     xvars=xvars.ha,
                                     ncores=ncores,
                                     site.lat = site.or.lat,
                                     xvar.name= "ha_abund")

apis.ha.c <- runCombinedParasiteModels(spec.data= spec.apis,
                                     species.group="apis",
                                     parasites = "CrithidiaPresence",
                                     xvars=xvars.ha.2,
                                     ncores=ncores,
                                     site.lat = site.or.lat,
                                     xvar.name= "ha_abund")

apis.d.a <- runCombinedParasiteModels(spec.data= spec.apis,
                                     species.group="apis",
                                     parasites = "ApicystisSpp",
                                     xvars=xvars.d,
                                     ncores=ncores,
                                     site.lat = site.or.lat,
                                     xvar.name= "degree")

apis.d.c <- runCombinedParasiteModels(spec.data= spec.apis,
                                    species.group="apis",
                                    parasites = "CrithidiaPresence",
                                    xvars=xvars.d.2,
                                    ncores=ncores,
                                    site.lat = site.or.lat,
                                    xvar.name= "degree")

apis.l <- runCombinedParasiteModels(spec.data= spec.apis,
                                     species.group="apis",
                                     xvars=xvars.l,
                                     ncores=ncores,
                                     site.lat = site.or.lat,
                                     xvar.name= "lat")

apis.a <- runCombinedParasiteModels(spec.data= spec.apis,
                                     species.group="apis",
                                     xvars=xvars.a,
                                     ncores=ncores,
                                     site.lat = site.or.lat,
                                     xvar.name= "area")

apis.sr <- runCombinedParasiteModels(spec.data= spec.apis,
                                     species.group="apis",
                                     xvars=xvars.sr,
                                     ncores=ncores,
                                     site.lat = site.or.lat,
                                     xvar.name= "SRDOY")

apis.y <- runCombinedParasiteModels(spec.data= spec.apis,
                                     species.group="apis",
                                     xvars=xvars.y,
                                     ncores=ncores,
                                     site.lat = site.or.lat,
                                     xvar.name= "year")


## **********************************************************
## Apis loo summaries
## **********************************************************
## crithidia
## **********************************************************
load("saved/parasiteFit_apis_CrithidiaPresenceApicystisSpp_apis_abundance_noFloralDiv.Rdata")
apis.ha <- fit.parasite
load("saved/parasiteFit_apis_CrithidiaPresenceApicystisSpp_diversity_noSRDoy.Rdata")
apis.div <- fit.parasite


loo.crithidia.ha2 <- loo(apis.ha2, resp="CrithidiaPresence")
loo.crithidia.div <- loo(apis.div, resp="CrithidiaPresence")

loo_compare(loo.crithidia.ha2, loo.crithidia.div)

## **********************************************************
## apicystis
## **********************************************************

load("saved/parasiteFit_apis_CrithidiaPresenceApicystisSpp_apis_abundance_noFloralDiv.Rdata")
apis.ha <- fit.parasite
load("saved/parasiteFit_apis_CrithidiaPresenceApicystisSpp_ba_noArea.Rdata")
apis.ba <- fit.parasite
load("saved/parasiteFit_apis_CrithidiaPresenceApicystisSpp_diversity.Rdata")
apis.div <- fit.parasite

loo.apicystis.ha2 <- loo(apis.ha2, resp="ApicystisSpp")
loo.apicystis.div <- loo(apis.div, resp="ApicystisSpp")
loo.apicystis.ba <- loo(apis.ba, resp="ApicystisSpp")

loo_compare(loo.apicystis.ha2, loo.apicystis.div, loo.apicystis.ba)
