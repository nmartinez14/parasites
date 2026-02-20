library(ggplot2)
library(brms)
library(dplyr)
library(tidyr)
library(tidybayes)
library(ggthemes)
library(bayestestR)
library(car)
library(DHARMa)
library(cmdstanr)

dir.create(path="saved", showWarnings = FALSE)
dir.create(path="saved/tables", showWarnings = FALSE)

load('../parasitesData.Rdata')

parasites <- c(
    "AscosphaeraSpp",
    "ApicystisSpp",
    "CrithidiaExpoeki",
    "CrithidiaBombi",
    "CrithidiaSpp",
    "CrithidiaMellificae",
    "NosemaBombi",
    "NosemaCeranae")


#instructions for installing cmdstanr
# options(repos = c(
# stan = "https://stan-dev.r-universe.dev",
# CRAN = "https://cran.rstudio.com"
# ))
# install.packages("cmdstanr")
# library(cmdstanr)
# #will take a while:
# cmdstanr::install_cmdstan(overwrite=TRUE)
# #check if installed
# cmdstanr::cmdstan_version()

# Rebuild the precompiled headers if update laptop
#rebuild_cmdstan()



