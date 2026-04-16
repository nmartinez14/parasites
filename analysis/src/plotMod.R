## Plotting function for the brms SEM
## This function takes the brms model output. You specify the predictors (x) and
## it generates a posterior prediction distribution using the conditional effects func
## Additional details provided include significance, color, and axis labels. 

plot_cond_effects <- function(model.fit, # brms model
                              this.response = "CrithidiaPresence", #response varible
                              this.effect = "Net_BeeDiversity", #predictor of interest
                              color.fill = "grey", #significance color, defaults to no significant
                              significance = "ns", # significance of the effect size, default not significant
                              x.label = "Bee Diversity", #characters for the figure labels
                              y.label = "Crithidia prevalence",
                              x.axis.lab = TRUE, #defaults TRUE, change if you want to mute the x axis label
                              data, #data frame for plotting the points
                              dat.x = "Net_BeeDiversity", #column name for plot our data for x
                              dat.y = "CrithidiaParasitismRate", #column name for plot our data for y
                              parasite = TRUE, # defaults to TRUE, if FALSE it won't include sample size of screened indiv
                              text.size = 12, # for formatting text size in fig
                              transform.x = NULL,# add a transformation when needed
                              angle.x = FALSE,
                              scale.y = FALSE){
  # Derive fill color and linetype from significance once, at the top
  fill.color <- switch(significance,
                       "ns"  = "grey60",
                       "95"  = "grey60",
                       "97"  = "#3182bd")
  
  line.type <- ifelse(significance %in% c("95", "97"), "solid", "dashed")
  
  # Extract CRI levels
  ce_95 <- conditional_effects(model.fit, effects = this.effect, prob = 0.95)
  ce_80 <- conditional_effects(model.fit, effects = this.effect, prob = 0.80)
  ce_50 <- conditional_effects(model.fit, effects = this.effect, prob = 0.50)
  
  effect.key <- paste(this.response, ".", this.response, "_", this.effect, sep = "")
  
  plot_data <- bind_rows(
    ce_95[[effect.key]] %>% mutate(level = "95%"),
    ce_80[[effect.key]] %>% mutate(level = "80%"),
    ce_50[[effect.key]] %>% mutate(level = "50%")
  )
  
  # Apply x transformation if provided
  if (!is.null(transform.x)) {
    plot_data <- plot_data %>%
      mutate(!!this.effect := transform.x(.data[[this.effect]]))
  }
  # Build base plot
  plot.obj <- ggplot(plot_data, aes(x = .data[[this.effect]], y = estimate__)) +
    geom_ribbon(aes(ymin = lower__, ymax = upper__, alpha = level), fill = fill.color) +
    scale_alpha_manual(values = c("95%" = 0.2, "80%" = 0.35, "50%" = 0.55),
                       name = "Credible interval") +
    guides(alpha = "none") + 
    geom_line(linewidth = 1.5, linetype = line.type) +
    labs(x = x.label, y = y.label) +
    theme_ms() +
    theme(axis.title.y = element_text(size = text.size),
          text = element_text(size = text.size),
          axis.title.x = if (x.axis.lab) element_text(size = text.size) else element_blank(),
          axis.text.x = if (angle.x) element_text(angle = 45, vjust = 1, hjust = 1, size = text.size) 
          else element_text(size = text.size))
  
  # Handle parasite toggling here
  # For parasite = TRUE
  if (parasite) {
    x_points <- if (!is.null(transform.x)) transform.x(data[[dat.x]]) else data[[dat.x]]
    plot.obj <- plot.obj +
      geom_jitter(data = data,
                  aes(y = .data[[dat.y]], x = x_points, color = SiteScreened),
                  width = 0.05) +
      scale_colour_gradient(low = "grey80", high = "grey20") +
      labs(color = "Screened individuals")+
      scale_x_continuous(
        breaks = scales::pretty_breaks(n = 6))+
      coord_cartesian(xlim = range(data[[dat.x]], na.rm = TRUE))
  } else {
    x_points <- if (!is.null(transform.x)) transform.x(data[[dat.x]]) else data[[dat.x]]
    y_points <- if (scale.y) scale(data[[dat.y]]) else data[[dat.y]]  
    jitter_data <- data %>% mutate(.x_points = x_points, .y_points = y_points)
    plot.obj <- plot.obj +
      geom_jitter(data = jitter_data,
                  aes(y = .y_points, x = .x_points),
                  cex = 2, color = "grey30") +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
      coord_cartesian(xlim = range(data[[dat.x]], na.rm = TRUE))
  }
  
  plot.obj
}
