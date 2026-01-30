## Function takes the models and creates a pdf with all the posterior predictive
## check plots. 

pp_check_all <- function(
    fit_bombus,
    fit_apis,
    xvar.name,
    top.level,
    save.dir = "figures/diagnostics/ppcheck",
    beta_binomial_responses = c("CrithidiaPresence", "ApicystisSpp"),
    ndraws = 50,
    density_type = "dens_overlay"
) {
  
  run_ppc <- function(fit, pollinator) {
    
    responses <- c('scaleMeanFloralDiversity', 'scaleNetBombusAbundance', 
                   'scaleNetHBAbundance', 'scaleNetBeeDiversity', 
                   'CrithidiaPresence', 'ApicystisSpp')
    
    file <- sprintf(
      "%s/ppcheck_%s_%s_%s.pdf",
      save.dir,
      pollinator,
      xvar.name,
      top.level
    )
    
    pdf(file, width = 7, height = 6)
    
    for (resp in responses) {
      
      message("Running pp_check for ", pollinator, " | response: ", resp)
      
      if (resp %in% beta_binomial_responses) {
        ppc_type <- "bars"
      } else {
        ppc_type <- density_type
      }
      
      p <- brms::pp_check(
        fit,
        resp = resp,
        type = ppc_type,
        ndraws = ndraws
      ) +
        ggplot2::ggtitle(
          paste0(
            pollinator, " SEM – Posterior predictive check: ", resp
          )
        )
      
      print(p)
    }
    
    dev.off()
    
    message("Saved PPCs to: ", file)
  }
  
  run_ppc(fit_bombus, pollinator = "Bombus")
  run_ppc(fit_apis,   pollinator = "Apis")
}

