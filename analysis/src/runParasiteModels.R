## This function creates and runs the models for the parasites.
runCombinedParasiteModels <- function(spec.data,## data
                                      parasites =
                                        c("CrithidiaPresence",
                                          "ApicystisSpp"),
                                      ## name of the parasites for the model
                                      xvars,## explanatory variables for the parasite model
                                      ncores=ncores, ## number of cores
                                      iter = 10^4,
                                      chains = 3,
                                      thin=1,
                                      init=0,
                                      data2 = NULL, ## Data for the Phylogeny
                                      SEM = TRUE,
                                      top.level,
                                      xvar.name,
                                      ...){
  
  ## Ensure weight columns exist and are properly formatted
  if(!"BombusWeights" %in% names(spec.data) | !"ApisWeights" %in% names(spec.data)){
    stop("BombusWeights and ApisWeights columns must exist in spec.data")
  }
  
  ## Convert to numeric and ensure no NAs in weight columns
  spec.data$BombusWeights <- as.numeric(spec.data$BombusWeights)
  spec.data$ApisWeights <- as.numeric(spec.data$ApisWeights)
  
  ## Replace any NAs with 0
  spec.data$BombusWeights[is.na(spec.data$BombusWeights)] <- 0
  spec.data$ApisWeights[is.na(spec.data$ApisWeights)] <- 0
  
  ## Print summary of weights
  print(paste("Bombus weights - 0s:", sum(spec.data$BombusWeights == 0), 
              "1s:", sum(spec.data$BombusWeights == 1)))
  print(paste("Apis weights - 0s:", sum(spec.data$ApisWeights == 0), 
              "1s:", sum(spec.data$ApisWeights == 1)))
  
  ## Create separate xvars for Bombus (all variables) and Apis (exclude last variable)
  xvars.bombus <- xvars
  xvars.apis <- xvars[-length(xvars)]  ## Remove the last element (phylogeny term)
  
  ## Create separate formula lists for Bombus and Apis weights
  bf.parasite.formulas.bombus <- vector(mode="list", length=length(parasites))
  bf.parasite.formulas.apis <- vector(mode="list", length=length(parasites))
  names(bf.parasite.formulas.bombus) <- paste0(parasites, "_Bombus")
  names(bf.parasite.formulas.apis) <- paste0(parasites, "_Apis")
  
  ## Create the models for parasites with BombusWeights (using all xvars)
  for(parasite in parasites){
    formula.bombus <- as.formula(paste(
      paste(parasite, "| subset(BombusWeights) + trials(1)"),
      paste(xvars.bombus, collapse=" + "),
      sep=" ~ "))
    bf.parasite.formulas.bombus[[paste0(parasite, "_Bombus")]] <- 
      bf(formula.bombus, family="beta_binomial")
    print(paste("Bombus formula for", parasite, ":"))
    print(formula.bombus)
  }
  
  ## Create the models for parasites with ApisWeights (excluding last xvar)
  for(parasite in parasites){
    formula.apis <- as.formula(paste(
      paste(parasite, "| subset(ApisWeights) + trials(1)"),
      paste(xvars.apis, collapse=" + "),
      sep=" ~ "))
    bf.parasite.formulas.apis[[paste0(parasite, "_Apis")]] <- 
      bf(formula.apis, family="beta_binomial")
    print(paste("Apis formula for", parasite, ":"))
    print(formula.apis)
  }
  
  ## Build separate SEM models for Bombus and Apis
  if(SEM){
    ## Bombus SEM model
    if(top.level == "lat"){
      print("lat - Bombus model")
      bform.bombus <- bf.fdiv.lat +
        bf.bombusabund.lat + bf.HBabund.lat +
        bf.bdiv.lat +
        bf.parasite.formulas.bombus[[1]] +
        bf.parasite.formulas.bombus[[2]] +
        set_rescor(FALSE)
      
      print("lat - Apis model")
      bform.apis <- bf.fdiv.lat +
        bf.bombusabund.lat + bf.HBabund.lat +
        bf.bdiv.lat +
        bf.parasite.formulas.apis[[1]] +
        bf.parasite.formulas.apis[[2]] +
        set_rescor(FALSE)
    } else if(top.level == "cp"){
      print("cum precip - Bombus model")
      bform.bombus <- bf.fdiv.cp +
        bf.bombusabund.cp + bf.HBabund.cp +
        bf.bdiv.cp +
        bf.parasite.formulas.bombus[[1]] +
        bf.parasite.formulas.bombus[[2]] +
        set_rescor(FALSE)
      
      print("cum precip - Apis model")
      bform.apis <- bf.fdiv.cp +
        bf.bombusabund.cp + bf.HBabund.cp +
        bf.bdiv.cp +
        bf.parasite.formulas.apis[[1]] +
        bf.parasite.formulas.apis[[2]] +
        set_rescor(FALSE)
    } else if(top.level == "area"){
      print("area - Bombus model")
      bform.bombus <- bf.fdiv.a +
        bf.bombusabund.a + bf.HBabund.a +
        bf.bdiv.a +
        bf.parasite.formulas.bombus[[1]] +
        bf.parasite.formulas.bombus[[2]] +
        set_rescor(FALSE)
      
      print("area - Apis model")
      bform.apis <- bf.fdiv.a +
        bf.bombusabund.a + bf.HBabund.a +
        bf.bdiv.a +
        bf.parasite.formulas.apis[[1]] +
        bf.parasite.formulas.apis[[2]] +
        set_rescor(FALSE)
    } else if(top.level == "fd"){
      print("floral div - Bombus model")
      bform.bombus <- bf.fdiv.a +
        bf.bombusabund.fd + bf.HBabund.fd +
        bf.bdiv.fd +
        bf.parasite.formulas.bombus[[1]] +
        bf.parasite.formulas.bombus[[2]] +
        set_rescor(FALSE)
      
      print("floral div - Apis model")
      bform.apis <- bf.fdiv.a +
        bf.bombusabund.fd + bf.HBabund.fd +
        bf.bdiv.fd +
        bf.parasite.formulas.apis[[1]] +
        bf.parasite.formulas.apis[[2]] +
        set_rescor(FALSE)
    }
  } else {
    ## Non-SEM models
    bform.bombus <- bf.parasite.formulas.bombus[[1]] +
      bf.parasite.formulas.bombus[[2]] +
      set_rescor(FALSE)
    
    bform.apis <- bf.parasite.formulas.apis[[1]] +
      bf.parasite.formulas.apis[[2]] +
      set_rescor(FALSE)
  }
  
  ## Fit Bombus model
  print("Fitting Bombus model...")
  fit.parasite.bombus <- brm(bform.bombus, spec.data,
                             cores= ncores,
                             backend = "cmdstanr",
                             iter = iter,
                             chains = chains,
                             thin=thin,
                             init=init,
                             control = list(adapt_delta = 0.999999999999, max_treedepth=12),
                             save_pars = save_pars(all = TRUE),
                             data2 = data2,
                             drop_unused_levels = TRUE,
                             ...)
  
  ## Fit Apis model
  print("Fitting Apis model...")
  fit.parasite.apis <- brm(bform.apis, spec.data,
                           cores= ncores,
                           backend = "cmdstanr",
                           iter = iter,
                           chains = chains,
                           thin=thin,
                           init=init,
                           control = list(adapt_delta = 0.999999999999, max_treedepth=12),
                           save_pars = save_pars(all = TRUE),
                           data2 = list(phylo_matrix = NULL),  ## No phylogeny for Apis
                           drop_unused_levels = TRUE,
                           ...)
  
  ## Calculate r2 values
  r2.bombus <- bayes_R2(fit.parasite.bombus)
  r2.apis <- bayes_R2(fit.parasite.apis)
  print("Bombus R2:")
  print(round(r2.bombus, 2))
  print("Apis R2:")
  print(round(r2.apis, 2))
  
  
  ## Save the Bombus model results
  save(fit.parasite.bombus, spec.data, r2.bombus,
       file=sprintf("saved/parasiteFit_Bombus_%s_%s_%s.Rdata",
                    paste(parasites, collapse=""),
                    xvar.name, top.level))
  
  ## Save the Apis model results
  save(fit.parasite.apis, spec.data, r2.apis,
       file=sprintf("saved/parasiteFit_Apis_%s_%s_%s.Rdata",
                    paste(parasites, collapse=""),
                    xvar.name, top.level))
  
  ## Plot the residuals for Bombus
  plot.res(fit.parasite.bombus, sprintf("Bombus_%s_%s_%s",
                                        paste(parasites, collapse=""),
                                        xvar.name, top.level))
  
  ## Plot the residuals for Apis
  plot.res(fit.parasite.apis, sprintf("Apis_%s_%s_%s",
                                      paste(parasites, collapse=""),
                                      xvar.name, top.level))
  

  ## Create tables with the results for Bombus
  write.ms.table(fit.parasite.bombus,
                 sprintf("parasitism_Bombus_%s_%s_%s",
                         paste(parasites, collapse=""),
                         xvar.name, top.level))
  
  ## Create tables with the results for Apis
  write.ms.table(fit.parasite.apis,
                 sprintf("parasitism_Apis_%s_%s_%s",
                         paste(parasites, collapse=""),
                         xvar.name, top.level))
  
  ## Make pp_check plots
  pp_check_all(
    fit_bombus = fit.parasite.bombus,
    fit_apis   = fit.parasite.apis,
    xvar.name  = xvar.name,
    top.level  = top.level
  )
  
  return(list(fit.bombus=fit.parasite.bombus,
              fit.apis=fit.parasite.apis,
              r2.bombus=r2.bombus,
              r2.apis=r2.apis,
              formulas_bombus=bf.parasite.formulas.bombus,
              formulas_apis=bf.parasite.formulas.apis))
}