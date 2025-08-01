################################################################################
#                     Recherche d'initialisations viables                      #
################################################################################

# Librairies utiles ------------------------------------------------------------

library(nimble)
# library(sf)
# library(tidyverse)

# fonction pour déterminer des valeurs initiales -------------------------------

valeurs_inits <- function(code, constants, data, params, inits) {
  # MCMC settings
  nc <- 2
  nburn <- 10000 
  ni <- nburn + 30000 
  nt <- 1
  
  model <- nimbleModel(
    code = code,
    constants = constants,
    data = data,
    inits = inits()
  )
  
  model$initializeInfo()
  
  comp <- 0
  res <- -Inf
  while ((comp < 100) & is.infinite(res)) {
    print(comp)
    comp <- comp + 1
    model$simulate()
    res <- model$calculate()
  }
  
  inits <- NULL
  for (p in params) {
    inits[[p]] <- model[[p]]
  }
  
  # inits <- list(
  #   beta0 = model$beta0,
  #   beta = model$beta, 
  #   alpha0 = model$alpha0,
  #   alpha1 = model$alpha1,
  #   tau = model$tau # à mettre ?!
  # )
  
  return(inits)
}


# Pour une liste de valeurs initiales différentes ------------------------------

list_inits <- function(code, constants, data, params, inits) {
  # MCMC settings
  nc <- 2
  nburn <- 10000
  ni <- nburn + 30000
  nt <- 1

  model <- nimbleModel(
    code = code,
    constants = constants,
    data = data,
    inits = inits()
  )

  inits <- list()
  inits[[1]] <- list()
  inits[[2]] <- list()
  
  for (i in 1:2) {
    comp <- 0
    res <- -Inf
    while ((comp < 500) & is.infinite(res)) {
      print(comp)
      comp <- comp + 1
      model$simulate()
      res <- model$calculate()
    }

    for (p in params) {
      inits[[i]][[p]] <- model[[p]]
    }
    
    # inits[[i]] <- list(
    #   beta0 = model$beta0,
    #   beta = model$beta,
    #   alpha0 = model$alpha0,
    #   alpha1 = model$alpha1,
    #   tau = model$tau # à mettre ?!
    }

  return(inits)
}

# # application au gros modèle sur la période ------------------------------------
# 
# donnees_utiles <- extract_data(grid_sf, periode) 
# 
# liste_inits <- list_inits(
#   code = code2,
#   constants = donnees_utiles$constants,
#   data = donnees_utiles$data,
#   params = params(2),
#   inits = inits2
# )
# 
# save(liste_inits, file = "RData/Inits/list_inits_modele2_50km2_2010-2024.RData")
