rm(list=ls())
source("lab_paths.R")
setwd(local.path)

load("parasites/parasitesData.Rdata")

library(tidyverse)
library(patchwork)  

# Vector of species to analyze
species_vector <- c("Bombus centralis", "Bombus huntii", 
                    "Apis mellifera", "Bombus bifarius")

# Set bootstrap parameters
n_bootstrap <- 1000
sample_size <- 5
set.seed(123)

# Initialize list to store ALL plots across all species
all_plots <- list()

# Loop over species
for (chosen_species in species_vector) {
  
  message(paste0("Processing species: ", chosen_species))
  
  # Filter data to chosen species
  species_data <- spec.net %>%
    filter(GenusSpecies == chosen_species & Apidae == 1) %>% 
    select(Site, SampleRound, Year, GenusSpecies, ParasitePresence)
  
  # Get unique sites for this species
  sites <- unique(species_data$Site)
  
  # Initialize lists for this species
  all_results <- list()
  all_summaries <- list()
  
  # Loop over all sites
  for (site in sites) {
    
    # Get unique sample rounds for this site
    sample_rounds <- species_data %>%
      filter(Site == site) %>%
      pull(SampleRound) %>%
      unique()
    
    # Loop over sample rounds within each site
    for (round in sample_rounds) {
      
      # Get unique years for this combination
      years <- species_data %>%
        filter(Site == site, SampleRound == round) %>%
        pull(Year) %>%
        unique()
      
      for (year in years) {
        
        # Get data for this site-round-year combination
        combo_data <- species_data %>%
          filter(Site == site, SampleRound == round, Year == year)
        
        # Need at least 6 individuals to compare
        if (nrow(combo_data) < 6) {
          message(paste0("Skipping ", chosen_species, " - ", site, " - Round ", round, 
                         " - Year ", year, ": only ", nrow(combo_data), " individuals"))
          next
        }
        
        # Calculate true prevalence (all individuals)
        true_prev <- mean(combo_data$ParasitePresence, na.rm = TRUE)
        n_total <- nrow(combo_data)
        
        # Bootstrap across different sample sizes
        sample_sizes <- 1:min(n_total, 15)
        
        bootstrap_results <- purrr::map_dfr(sample_sizes, function(n) {
          
          size_results <- purrr::map_dfr(1:n_bootstrap, function(i) {
            sampled <- combo_data %>%
              sample_n(n, replace = FALSE)
            
            tibble(
              iteration = i,
              sample_size = n,
              prevalence = mean(sampled$ParasitePresence, na.rm = TRUE)
            )
          })
          
          return(size_results)
        })
        
        # Add metadata
        bootstrap_results <- bootstrap_results %>%
          mutate(
            Site = site,
            SampleRound = round,
            Year = year,
            Species = chosen_species,
            true_prevalence = true_prev,
            n_total = n_total,
            abs_diff = abs(prevalence - true_prev)
          )
        
        # Store results
        combo_name <- paste0(chosen_species, "_", site, "_Round", round, "_Year", year)
        all_results[[combo_name]] <- bootstrap_results
        
        # Calculate summary by sample size
        summary_stats <- bootstrap_results %>%
          group_by(sample_size) %>%
          summarise(
            Site = first(Site),
            SampleRound = first(SampleRound),
            Year = first(Year),
            Species = first(Species),
            n_total = first(n_total),
            true_prev = first(true_prevalence),
            mean_prevalence = mean(prevalence),
            sd_prevalence = sd(prevalence),
            mean_abs_diff = mean(abs_diff),
            .groups = "drop"
          )
        
        all_summaries[[combo_name]] <- summary_stats
      }
    }
  }
  
  # Combine all results for this species
  all_bootstrap_data <- bind_rows(all_results)
  all_summary_data <- bind_rows(all_summaries)
  
  # Calculate bootstrap mean prevalence for each site-round-year-sample_size combination
  bootstrap_means <- all_bootstrap_data %>%
    group_by(Site, SampleRound, Year, sample_size) %>%
    summarise(
      bootstrap_mean_prev = mean(prevalence),
      .groups = "drop"
    )
  
  # Calculate observed prevalence (from all individuals) and join with bootstrap means
  observed_data <- species_data %>%
    filter(GenusSpecies == chosen_species) %>%
    group_by(Site, SampleRound, Year) %>%
    summarise(
      observed_prevalence = mean(ParasitePresence, na.rm = TRUE),
      true_sample_size = n(),
      .groups = "drop"
    )
  
  # Join and calculate prevalence difference for each bootstrap iteration
  plot_data <- all_bootstrap_data %>%
    left_join(bootstrap_means, by = c("Site", "SampleRound", "Year", "sample_size")) %>%
    left_join(observed_data, by = c("Site", "SampleRound", "Year")) %>%
    mutate(
      prevalence_diff = observed_prevalence - bootstrap_mean_prev
    )
  
  # Create boxplot
  p <- ggplot(plot_data, aes(x = factor(true_sample_size), y = prevalence_diff)) +
    geom_boxplot() +
    labs(
      title = paste0("Species: ", chosen_species),
      x = "True Sample Size (total individuals sampled)",
      y = "Prevalence Difference\n(Observed - Bootstrap Mean)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "italic"),
      axis.text.x = element_text(angle = 0, hjust = 0.5)
    )
  
  # Store plot
  all_plots[[chosen_species]] <- p
}

# Combine all plots into a panel and save as PDF

combined_plot <- wrap_plots(all_plots, ncol = 2,
                            labels = c("A", "B", "C", "D"))

print(combined_plot)
ggsave(combined_plot, file="parasites/figures/bootstrap_diff_spp.pdf",
       height=10, width=12)

dev.off()


  

  