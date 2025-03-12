# parasites
Analysis for understanding how natural diversity gradients influence parasite 
prevalence in bees.

## Repository organization  

packages.sh  
- This shell script includes code to install all of the required R packages to 
run analysis and figure-generating scripts. 

parasitesData.Rdata  
- This Rdata file includes all of the necessary pieces of data to run the 
scripts in this repository. These include:
  - spec.net  
    - This data frame includes all of the specimen-level data and relevant 
    explanatory variables for all analyses. 
  - phylo_matrix  
    - This matrix array is the phylogenetic covariance matrix of bee species 
    included in this manuscript. The original published supermatrix was 
    downloaded and trimmed using our species list. 
    (Henriquez Piskulich, Patricia Andrea; Hugall, Andrew F.; Stuart-Fox, Devi (2023). 
    A supermatrix phylogeny of the world’s bees (Hymenoptera: Anthophila) 
    [Dataset]. Dryad. https://doi.org/10.5061/dryad.80gb5mkw1
    )  
  - sites_sf  
    - This sf spatial object includes spatial data from our sites used for 
    generating the map figure.  
  
analysis  
  - This folder includes the scripts to run the two main analyses for 
  this manuscript, including:  
    - SEM  
      - 1_amplification_dilution.R  
        - This script runs the structural equation modeling script for each 
        host genus and parasite.
      - 2_Plotting_SEM.R  
        - This script creates the scatterplots for the results from the parasite
        prevalence models of the SEM.
      - 3_Plotting_Community.R  
        - This script creates the scatterplots for the results from the 
        community models portion of the SEM.
    
figures  
  - This folder includes the scripts to generate the site map.  

## To use this repository  
  The R scripts in the analysis folder should be run in order 
(i.e., 0, 1, 2) to reproduce the analyses and figures. All analyses were run 
using R (version 4.4.1). 