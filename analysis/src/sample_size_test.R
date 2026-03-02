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
target_n <- 5  # fixed sample size per the analysis description
set.seed(123)


# Initialize storage for all species-level bootstrap data
all_species_bootstrap_data <- list()

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
  
  # Initialize list for this species
  all_results <- list()
  
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
        
        # Must have more than 5 individuals to perform the resampling test
        if (nrow(combo_data) <= target_n) {
          message(paste0("Skipping ", chosen_species, " - ", site, " - Round ", round, 
                         " - Year ", year, ": only ", nrow(combo_data), " individuals (need > ", target_n, ")"))
          next
        }
        
        # Compute p_full from ALL individuals in this stratum
        p_full <- mean(combo_data$ParasitePresence, na.rm = TRUE)
        n_total <- nrow(combo_data)
        
        # Repeatedly draw exactly target_n individuals without replacement
        bootstrap_results <- purrr::map_dfr(1:n_bootstrap, function(i) {
          sampled <- combo_data %>%
            sample_n(target_n, replace = FALSE)
          
          p_5 <- mean(sampled$ParasitePresence, na.rm = TRUE)
          
          tibble(
            iteration     = i,
            p_5           = p_5,
            p_full        = p_full,
            abs_diff      = abs(p_full - p_5),
            signed_diff   = p_full - p_5,  # for checking centering around 0
            Site          = site,
            SampleRound   = round,
            Year          = year,
            Species       = chosen_species,
            n_total       = n_total
          )
        })
        
        # Store results
        combo_name <- paste0(chosen_species, "_", site, "_Round", round, "_Year", year)
        all_results[[combo_name]] <- bootstrap_results
      }
    }
  }
  
  # Combine all results for this species
  all_bootstrap_data <- bind_rows(all_results)
  
  
  if (nrow(all_bootstrap_data) == 0) {
    message(paste0("No valid strata for species: ", chosen_species, ". Skipping plot."))
    next
  }
  
  # Create a stratum label for grouping on the x-axis
  all_bootstrap_data <- all_bootstrap_data %>%
    mutate(stratum = paste0(Site, "\nRound ", SampleRound, "\n", Year,
                            "\n(n=", n_total, ")"))
  
  # Store species-level tibble
  all_species_bootstrap_data[[chosen_species]] <- all_bootstrap_data
  
  # Plot: distribution of signed differences (p_full - p_5) per stratum
  # Centered around 0 = no meaningful bias from sampling 5 individuals
  p <- ggplot(all_bootstrap_data, aes(x = stratum, y = signed_diff)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.7) +
    geom_boxplot(outlier.size = 0.8, fill = "grey85") +
    labs(
      title    = paste0("Species: ", chosen_species),
      subtitle = paste0("Prevalence stability: ", n_bootstrap, 
                        " draws of n=", target_n, " vs. p_full (all individuals)"),
      x        = "Stratum (Site · Round · Year · total n)",
      y        = expression(p[full] - p[5])
    ) +
    theme_minimal() +
    theme(
      plot.title    = element_text(size = 12, face = "italic"),
      plot.subtitle = element_text(size = 9, color = "grey40"),
      axis.text.x   = element_text(angle = 45, hjust = 1, size = 7)
    )
  
  # Store plot
  all_plots[[chosen_species]] <- p
  
}

final_bootstrap_tibble <- bind_rows(all_species_bootstrap_data, .id = "SpeciesName")

summary <- final_bootstrap_tibble %>%
  mutate(mean_abs_error_total = mean(abs_diff, na.rm = TRUE),
         mean_signed_error_total = mean(signed_diff, na.rm = TRUE)) %>% 
  group_by(Species) %>% 
  summarize(
    mean_abs_error = mean(abs_diff),
    mean_signed_error = mean(signed_diff),
    n_combo = n_distinct(paste(Site, SampleRound, Year))
  ) 
  
# Combine all plots into a panel and save as PDF

combined_plot <- wrap_plots(all_plots, ncol = 2,
                            labels = c("A", "B", "C", "D"))

print(combined_plot)
ggsave(combined_plot, file="parasites/figures/bootstrap_diff_spp.jpg",
       height=10, width=12)

dev.off()


  

  