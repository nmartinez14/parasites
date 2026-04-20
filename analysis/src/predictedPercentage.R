get_prevalence_range <- function(model.fit,
                                 this.response = "CrithidiaPresence",
                                 this.effect = "Net_BeeDiversity",
                                 data,
                                 dat.x = "Net_BeeDiversity") {
  
  # Build prediction grid at just min and max of predictor
  pred_grid <- model.fit$data %>%
    summarise(
      across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
      across(where(is.factor), ~ names(sort(table(.x), decreasing = TRUE))[1]),
      across(where(is.character), ~ names(sort(table(.x), decreasing = TRUE))[1])
    ) %>%
    slice(rep(1, 2)) %>%  # just 2 rows - min and max
    mutate(!!dat.x := c(min(data[[dat.x]], na.rm = TRUE),
                        max(data[[dat.x]], na.rm = TRUE)),
           label = c("min", "max"))
  
  # Draw from posterior
  draws <- posterior_epred(model.fit,
                           newdata = pred_grid,
                           resp = this.response,
                           re_formula = NA)
  
  # Compute median and 95% CRI at min and max
  results <- data.frame(
    predictor = c("min", "max"),
    predictor_value = c(min(data[[dat.x]], na.rm = TRUE),
                        max(data[[dat.x]], na.rm = TRUE)),
    estimate  = apply(draws, 2, median) * 100,  # convert to percentage
    lower_95  = apply(draws, 2, quantile, probs = 0.025) * 100,
    upper_95  = apply(draws, 2, quantile, probs = 0.975) * 100
  )
  
  # Print in a sentence-ready format
  cat(sprintf(
    "Predicted prevalence at minimum %s (%.2f): %.1f%% [95%% CRI: %.1f%%-%.1f%%]\n",
    dat.x, results$predictor_value[1],
    results$estimate[1], results$lower_95[1], results$upper_95[1]
  ))
  cat(sprintf(
    "Predicted prevalence at maximum %s (%.2f): %.1f%% [95%% CRI: %.1f%%-%.1f%%]\n",
    dat.x, results$predictor_value[2],
    results$estimate[2], results$lower_95[2], results$upper_95[2]
  ))
  
  invisible(results)  # also return the data frame silently in case you need the numbers
}

get_prevalence_range(fit.apis.l,
                     this.response = "ApicystisSpp",
                     this.effect = "Lat",
                     data = spec.uni,
                     dat.x = "Lat")
