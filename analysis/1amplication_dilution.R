rm(list=ls())
source("lab_paths.R")
setwd(local.path)

ncores <- 3

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
## social species abundances
xvars.ss <-  c("Net_BeeDiversity",
               "Net_BombusAbundance",
               "Net_HBAbundance",
               "rare.degree",
               "MeanFloralDiversity",
               #"SRDoy",
               "Lat","Area",
               "(1|Site)",
               "(1|gr(GenusSpecies, cov = phylo_matrix))")

## bumble bee abundance only (vif sometimes indicates HB abundance and
## bombus abundance are colinear)
xvars.ba <- xvars.ss[xvars.ss != "Net_HBAbundance"]

## honey bee abundance only (vif sometimes indicates floral diversity and
## bombus abundance are colinear)
xvars.ha <- xvars.ss[xvars.ss != "Net_BombusAbundance"]

xvars.ha.2 <- xvars.ha[xvars.ha != "MeanFloralDiversity"]

## diversity only
xvars.div <- xvars.ss[!xvars.ss %in%
                        c("Net_HBAbundance", "Net_BombusAbundance")]

## **********************************************************
## community model, check assumptions first before adding parasites
## **********************************************************

run_plot_freq_model_diagnostics(remove_subset_formula(formula.flower.div),
                                this_data=spec.net[spec.net$Weights == 1,],
                                this_family="students", site.lat=site.or.lat)

run_plot_freq_model_diagnostics(remove_subset_formula(formula.bee.abund),
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

xvar.order <- c("social_species",
                "bombus_abundance",
                "apis_abundance",
                "diversity")

## social species abundance as x var
bombus.ss <- runCombinedParasiteModels(spec.data= spec.bombus,
                                       species.group="bombus",
                                       xvars=xvars.ss,
                                       ncores=ncores,
                                       data2= list(phylo_matrix=phylo_matrix),
                                       site.lat=site.or.lat,
                                       xvar.name=xvar.order[1])


## HB and bombus abundance are colinear in crithidia model, so not a valid model

## bombus abundance only
bombus.ba <- runCombinedParasiteModels(spec.data= spec.bombus,
                                       species.group="bombus",
                                       xvars=xvars.ba,
                                       ncores=ncores,
                                       data2= list(phylo_matrix=phylo_matrix),
                                       site.lat=site.or.lat,
                                       xvar.name=xvar.order[2])

## apis abundance only
bombus.ha <- runCombinedParasiteModels(spec.data= spec.bombus,
                                       species.group="bombus",
                                       xvars=xvars.ha,
                                       ncores=ncores,
                                       data2= list(phylo_matrix=phylo_matrix),
                                       site.lat=site.or.lat,
                                       xvar.name=xvar.order[3])
## For Crithidia bee div collinear with SRDOY(not quite VIF > 5)


## no abundances, just bee diversity 
bombus.div <- runCombinedParasiteModels(spec.data= spec.bombus,
                                       species.group="bombus",
                                       xvars=xvars.div,
                                       ncores=ncores,
                                       data2= list(phylo_matrix=phylo_matrix),
                                       site.lat=site.or.lat,
                                       xvar.name=xvar.order[4])



## **********************************************************
## Bombus loo summaries
## **********************************************************
## crithidia
## **********************************************************
## not including ss because HB abundance and bombus abundance are colinear

load("saved/parasiteFit_bombus_CrithidiaPresenceApicystisSpp_lat_bombus_abundance.Rdata")
bombus.ba <- fit.parasite
load("saved/parasiteFit_bombus_CrithidiaPresenceApicystisSpp_lat_social_species.Rdata")
bombus.ss <- fit.parasite
load("saved/parasiteFit_bombus_CrithidiaPresenceApicystisSpp_lat_apis_abundance.Rdata")
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
## Including ba (bombus abundance) since HB and bombus abundance
## abundance are not colinear in these models

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


xvar.order <- c("social_species",
                "apis_abundance",
                "bombus_abundance",
                "diversity")

## social species abundance as x var
apis.ss <- runCombinedParasiteModels(spec.data= spec.apis,
                                     species.group="apis",
                                     xvars=xvars.ss[-length(xvars.ss)],
                                     ncores=ncores,
                                     site.lat=site.or.lat,
                                     xvar.name=xvar.order[1])
## Multiple variables are collinear for both parasites, so not a valid model

## apis abundance
apis.ha <- runCombinedParasiteModels(spec.data= spec.apis,
                                     species.group="apis",
                                     xvars=xvars.ha[-length(xvars.ha)],
                                     ncores=ncores,
                                     site.lat=site.or.lat,
                                     xvar.name=xvar.order[2])
## Multiple variables are collinear for both parasites, so not a valid model

## apis abundance - floral diversity
apis.ha2 <- runCombinedParasiteModels(spec.data= spec.apis,
                                     species.group="apis",
                                     xvars=xvars.ha.2[-length(xvars.ha.2)],
                                     ncores=ncores,
                                     site.lat=site.or.lat,
                                     xvar.name=xvar.order[2])

## bombus abundance
apis.ba <- runCombinedParasiteModels(spec.data= spec.apis,
                                      species.group="apis",
                                      xvars=xvars.ba[-length(xvars.ba)],
                                      ncores=ncores,
                                      site.lat=site.or.lat,
                                      xvar.name=xvar.order[3])
## All variables collinear for crithidia. 
## For apicystis lat, area and bombus abundance collinear


## no abundances, just bee diversity 
apis.div <- runCombinedParasiteModels(spec.data= spec.apis,
                                      species.group="apis",
                                      xvars=xvars.div[-length(xvars.div)],
                                      ncores=ncores,
                                      data2= list(phylo_matrix=phylo_matrix),
                                      site.lat=site.or.lat,
                                      xvar.name=xvar.order[4])
## For Crithida model floral diversity/Lat/SRDoy collinear.

## **********************************************************
## Apis loo summaries
## **********************************************************
## crithidia
## **********************************************************
load("saved/parasiteFit_apis_CrithidiaPresenceApicystisSpp_lat_apis_abundance.Rdata")
apis.ha2 <- fit.parasite
load("saved/parasiteFit_apis_CrithidiaPresenceApicystisSpp_lat_diversity.Rdata")
apis.div <- fit.parasite
load("saved/parasiteFit_apis_CrithidiaPresenceApicystisSpp_lat_bombus_abundance.Rdata")
apis.ba <- fit.parasite

loo.crithidia.ha2 <- loo(apis.ha2, resp="CrithidiaPresence")
loo.crithidia.div <- loo(apis.div, resp="CrithidiaPresence")

loo_compare(loo.crithidia.ha2, loo.crithidia.div)

## **********************************************************
## apicystis
## **********************************************************

loo.apicystis.ha2 <- loo(apis.ha2, resp="ApicystisSpp")
loo.apicystis.div <- loo(apis.div, resp="ApicystisSpp")
loo.apicystis.ba <- loo(apis.ba, resp="ApicystisSpp")

loo_compare(loo.apicystis.ha2, loo.apicystis.div, loo.apicystis.ba)

## **********************************************************
## save models and loo results
## **********************************************************

write.csv(rbind(loo.bombus.crithidia,
                loo.apis.crithidia,
                loo.bombus.apicystis,
                loo.apis.apicystis),
          row.names=FALSE,
          file="saved/loo.csv")

write.table(rbind(loo.bombus.crithidia,
                  loo.apis.crithidia,
                  loo.bombus.apicystis,
                  loo.apis.apicystis),
            row.names=FALSE,
            file="saved/loo.txt", sep= " & ")

save(loo.apis.apicystis,
     loo.apis.crithidia,
     loo.bombus.apicystis,
     loo.bombus.crithidia,
     fit.apis.all,
     fit.bombus.all,
     fit.apis.ss,
     fit.bombus.ss,
     fit.apis.ha,
     fit.bombus.ba,
     file="saved/all_models_loo.R")
