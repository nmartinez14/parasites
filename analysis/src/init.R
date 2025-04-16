library(ggplot2)
library(brms)
library(dplyr)
library(tidyr)
library(tidybayes)
library(ggthemes)
library(bayestestR)
library(car)
library(DHARMa)

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






