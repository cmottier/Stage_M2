################################################################################
#          Ajustement du modèle de régression - approche bayésienne            #
################################################################################

# /!\ Modèle à sélectionner : effort voulu et nombre d'observations par cellule
# 
# Lasso bayésien implémenté via le choix du prior des coefficients beta


# Librairies utiles ------------------------------------------------------------

library(nimble)
library(sf)
library(tidyverse)
library(MCMCvis)


# Chargement et travail préalable sur grid -------------------------------------

# Chargement de la grille complète
load("grid_sf_5km2.RData")

# # on garde les cellules qui intersectent à au moins 10%
# grid_selec <- grid_selec %>%
#   filter(as.numeric(area)>0.1*max(as.numeric(area)))

# # on renumérote les cellules de la sélection
# grid_selec$grid_id = 1:lengths(grid_selec)[1] 


# Codes Nimble -----------------------------------------------------------------

# Une détection seule par cellule

code_uni <- nimbleCode({
  # intensité et effort pour toutes les cellules
  for(pixel in 1:npixel){
    # intensité
    log(lambda[pixel]) <- beta[1] +
      beta[2] * x_1[pixel] +
      beta[3] * x_2[pixel] +
      beta[4] * x_3[pixel] + 
      beta[5] * x_4[pixel] + 
      beta[6] * x_5[pixel] +
      cell_area[pixel] 
    # effort
    logit(b[pixel]) <-  alpha[1] + alpha[2] * h_1[pixel]
  }
  
  # loi jointe, pour les m cellules contenant une détection (po_pixel)
  po_denominator <- inprod(lambda[1:npixel], b[1:npixel]) / m
  
  for(po in 1:m){
    ones[po] ~ dbern(
      exp(
        log(lambda[po_pixel[po]] * b[po_pixel[po]]) -
          po_denominator)
      / CONSTANT) 
  }
  
  # Priors 
  beta[1] ~ dnorm(0, sd = 2) # intercept
  for(i in 2:6){
    beta[i] ~ ddexp(0, tau) 
  }
  tau ~ dunif(0.001,10)
  
  # for(i in 1:6){
  #   beta[i] ~ dnorm(0, sd = 2)
  # }

  for(j in 1:2){
    alpha[j] ~ dnorm(0, sd = 2)
  }
  
})


# Multiples détections par cellule

code_multi <- nimbleCode({
  # intensité et effort pour toutes les cellules
  for(pixel in 1:npixel){
    # intensité
    log(lambda[pixel]) <- beta[1] +
      beta[2] * x_1[pixel] +
      beta[3] * x_2[pixel] +
      beta[4] * x_3[pixel] + 
      beta[5] * x_4[pixel] + 
      beta[6] * x_5[pixel] + 
      cell_area[pixel]
    # effort
    logit(b[pixel]) <-  alpha[1] + alpha[2] * h_1[pixel]
  }
  
  # pour les nobs observations (obs_pixel)
  obs_denominator <- inprod(lambda[1:npixel], b[1:npixel]) / nobs
  
  for(obs in 1:nobs){
    ones[obs] ~ dbern(
      exp(
        log(lambda[obs_pixel[obs]] * b[obs_pixel[obs]]) -
          obs_denominator)   
      / CONSTANT) 
  }
  
  # Priors 
  beta[1] ~ dnorm(0, sd = 2) # intercept
  for(i in 2:6){
    beta[i] ~ ddexp(0, tau) 
  }
  tau ~ dunif(0.001,10)
  
  # for(i in 1:6){
  #   beta[i] ~ dnorm(0, sd = 2)
  # }

  for(j in 1:2){
    alpha[j] ~ dnorm(0, sd = 2)
  }
  
})



# Estimation -------------------------------------------------------------------

## Fonction pour choisir le modèle  ##################

#' Estimation des paramètres (IPP)
#'
#' @param grid # grille à utiliser
#' @param modele # 1 : une détection, 2 : multiples détections
#' @param effort # 'prox' : proximité aux routes, 'gbif': densité d'observations GBIF
#' @param annee 
#'
estim_param <- function(grid, modele, effort, annee) {
  
  # Nombre de pixels
  npix <- nrow(grid)
  
  # Aire des pixels
  s.area <- as.numeric(units::set_units(grid$area,"km^2"))
  logarea <- log(s.area)
  
  # ID des cellules où il y a au moins une occurrence
  pixel.id.det <- grid$grid_id[grid[[paste0("nnutria", annee)]] > 0]
  
  # Troncature des observations
  nb_observations <- grid[[paste0("nnutria", annee)]] 
  nb_observations[nb_observations > 50] <- 50 # valeur arbitraire à définir
  
  # nombre total d'observations prises en compte
  nobs = sum(nb_observations)
  
  # pixel associé aux observations (avec répétition)
  obs_pixel <- NULL
  for (i in 1:length(pixel.id.det)){
    obs_pixel <- c(obs_pixel, rep(pixel.id.det[i], nb_observations[pixel.id.det[i]]))
  }
  
  # Variables 
  data <- list(cell_area = logarea,
               x_1 = scale(grid$dist_eau)[,1],
               x_2 = scale(grid$logdensity)[,1],
               x_3 = scale(grid$agri_cover)[,1],
               x_4 = scale(grid[[paste0("pcum_", annee)]])[,1],
               x_5 = scale(grid[[paste0("tmin_", annee)]])[,1]
               )
  
  if (effort == 'prox') {
    data$h_1 <- scale(grid$dist_acces)[, 1]
  } else {
    data$h_1 <- scale(grid[[paste0("dgbif_", annee)]])[, 1]
  }
  
  if (modele == 1) { 
    data$ones <- rep(1, length(pixel.id.det)) 
    } else {
    data$ones <- rep(1, nobs)
  }

  
  constants <- list(
    npixel = npix,
    nobs = nobs,
    m = length(pixel.id.det), 
    CONSTANT = 50000,
    po_pixel = pixel.id.det,
    obs_pixel = obs_pixel
    ) 
  
  # Initialisation
  inits <- function(){
    list(
      beta = rnorm(6, 0, 1), 
      alpha = rnorm(2, 0, 1)
    )
  }
  
  # Paramètres à suivre
  params <- c("alpha", "beta") 
  
  # MCMC settings
  nc <- 2
  nburn <- 10000 
  ni <- nburn + 30000 
  nt <- 1
  
  # MCMC
  if (modele == 1) {code <- code_uni} else {code <- code_multi}
  
  out <- nimbleMCMC(
    code = code,
    constants = constants,
    data = data,
    inits = inits(),
    monitors = params,
    niter = ni,
    nburnin = nburn,
    nchains = nc,
    thin = nt,
    # WAIC = TRUE
  )
  
  # sortie
  return(out)
}


## Lancement et sauvegarde #################

periode = 2016:2016

for (annee in periode) {
  print(annee)
  set.seed(123)
  # assign(
  #   x = paste0("outMCMC_", annee),
  #   value = estim_param(
  #     grid = grid_sf,
  #     modele = 1,
  #     effort = "gbif",
  #     annee = annee
  #   )
  # )
  out <- estim_param(
    grid = grid_sf,
    modele = 2,
    effort = "gbif",
    annee = annee
  )
  save(out, file = paste0("out_multi_gbif_", annee, ".RData"))
}

# save.image(file = "out_uni_gbif_2021_2021.RData")
