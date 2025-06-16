################################################################################
#                           Modèle multi-années                                #
################################################################################

# Librairies utiles ------------------------------------------------------------

library(nimble)
library(sf)
library(tidyverse)


# Construction des données utiles aux modèles ----------------------------------

extract_data <- function(grid, periode) {
  # Nombre d'années
  nb_annees <- length(periode)
  
  # Nombre de pixels
  npix <- nrow(grid)
  
  # Aire des pixels
  s.area <- as.numeric(units::set_units(grid$area,"km^2"))
  logarea <- log(s.area)
  
  # Observations annuelles
  pixel.id.det <- NULL
  nb_observations <- NULL
  nobs <- NULL
  obs_pixel <- NULL
  
  for (a in 1:nb_annees) {
    # ID des cellules où il y a au moins une occurrence
    pixel.id.det[[a]] <- grid$grid_id[grid[[paste0("nnutria", periode[a])]] > 0]
    
    # Troncature des observations
    nb_observations[[a]] <- grid[[paste0("nnutria", periode[a])]] 
    nb_observations[[a]][nb_observations[[a]] > 50] <- 50 # valeur arbitraire à définir
    
    # nombre total d'observations prises en compte
    nobs[a] <- sum(nb_observations[[a]])
  }
  
  # pixel associé aux observations (avec répétition)
  obs_pixel <- matrix(NA, nrow = max(nobs), ncol = nb_annees)
  for (a in 1:nb_annees) {
    s <- 0
    for (i in (1:length(pixel.id.det[[a]]))){
      obs_pixel[(s+1):(s + nb_observations[[a]][pixel.id.det[[a]][i]]), a] <- pixel.id.det[[a]][i]
      s <- s + nb_observations[[a]][pixel.id.det[[a]][i]]
    }
  }
  
  # matrice de 1 correspondant aux observations (one's trick)
  ones <- matrix(NA, nrow = max(nobs), ncol = nb_annees)
  for (a in 1:nb_annees) {
    ones[1:nobs[a], a] <- 1
  }
  
  pcum_periode <- grid_sf %>% select(paste0("pcum_" , periode)) %>% st_drop_geometry %>% as.matrix()
  tmin_periode <- grid_sf %>% select(paste0("tmin_" , periode)) %>% st_drop_geometry %>% as.matrix()
  dgbif_periode <- grid_sf %>% select(paste0("dgbif_" , periode)) %>% st_drop_geometry %>% as.matrix()
  
  # Variables 
  data <- list(cell_area = logarea,
               x_1 = scale(grid$dist_eau)[,1],
               x_2 = scale(grid$logdensity)[,1],
               x_3 = scale(grid$agri_cover)[,1],
               x_4 = (pcum_periode - mean(pcum_periode))/sd(pcum_periode),
               x_5 = (tmin_periode - mean(tmin_periode))/sd(tmin_periode),
               h_1 = (dgbif_periode - mean(dgbif_periode))/sd(dgbif_periode),
               ones = ones
  )
  
  # str(data,1)
  
  constants <- list(
    npixel = npix,
    nb_annees = nb_annees,
    nobs = nobs,
    CONSTANT = 50000,
    obs_pixel = obs_pixel
  ) 
  
  return(list(data = data, constants = constants))
}  



# Code Nimble ------------------------------------------------------------------
# Multiples détections par cellule

## Modèle (M0) avec tout constant ##############################################

code0 <- nimbleCode({
  # intensité et effort pour toutes les cellules
  for(pixel in 1:npixel){
    for (a in 1:nb_annees){
      # intensité
      log(lambda[pixel,a]) <- beta0 +
        beta[1] * x_1[pixel] +
        beta[2] * x_2[pixel] +
        beta[3] * x_3[pixel] + 
        beta[4] * x_4[pixel,a] + 
        beta[5] * x_5[pixel,a] + 
        cell_area[pixel]
      # effort
      logit(b[pixel,a]) <-  alpha0[a] + alpha1 * h_1[pixel,a]
    }
  }
  
  # pour les nobs observations (obs_pixel)
  for (a in 1:nb_annees) {
    obs_denominator[a] <- inprod(lambda[1:npixel,a], b[1:npixel,a]) / nobs[a]
    
    for(obs in 1:nobs[a]){
      ones[obs, a] ~ dbern(
        exp(
          log(lambda[obs_pixel[obs,a],a] * b[obs_pixel[obs,a],a]) -
            obs_denominator[a])   
        / CONSTANT) 
    }
  }
  
  # Priors 
  beta0  ~ dnorm(0, sd = 2) 
  
  for(i in 1:5){
    beta[i] ~ ddexp(0, tau) 
  }
  tau ~ dunif(0.001,10)
  
  for (a in 1:nb_annees) {
    alpha0[a] ~ dnorm(mu.a, sd = sd.a)
  }
  mu.a ~ dnorm(0, sd = 2)
  sd.a ~ dunif(0, 5)
  
  alpha1 ~ dnorm(0, sd = 2)
  
})



## Modèle avec intercepts annuels indépendants #################################
# intercepts libres (beta0)

code1 <- nimbleCode({
  # intensité et effort pour toutes les cellules
  for(pixel in 1:npixel){
    for (a in 1:nb_annees){
      # intensité
      log(lambda[pixel,a]) <- beta0[a] +
        beta[1] * x_1[pixel] +
        beta[2] * x_2[pixel] +
        beta[3] * x_3[pixel] + 
        beta[4] * x_4[pixel,a] + 
        beta[5] * x_5[pixel,a] + 
        cell_area[pixel]
      # effort
      logit(b[pixel,a]) <-  alpha0[a] + alpha1 * h_1[pixel,a]
    }
  }
  
  # pour les nobs observations (obs_pixel)
  for (a in 1:nb_annees) {
    obs_denominator[a] <- inprod(lambda[1:npixel,a], b[1:npixel,a]) / nobs[a]
    
    for(obs in 1:nobs[a]){
      ones[obs, a] ~ dbern(
        exp(
          log(lambda[obs_pixel[obs,a],a] * b[obs_pixel[obs,a],a]) -
            obs_denominator[a])   
        / CONSTANT) 
    }
  }
  
  # Priors 
  for (a in 1:nb_annees) {
    beta0[a] ~ dnorm(0, sd = 2) 
  }
  
  for(i in 1:5){
    beta[i] ~ ddexp(0, tau) 
  }
  tau ~ dunif(0.001,10)
  
  for (a in 1:nb_annees) {
    alpha0[a] ~ dnorm(mu.a, sd = sd.a)
  }
  mu.a ~ dnorm(0, sd = 2)
  sd.a ~ dunif(0, 5)
  
  alpha1 ~ dnorm(0, sd = 2)
  
})


## Modèle avec marche aléatoire sur intercepts #################################
# Cf Outhwaite
# Intercepts dépendants des années précédentes

code2 <- nimbleCode({
  # intensité et effort pour toutes les cellules
  for(pixel in 1:npixel){
    for (a in 1:nb_annees){
      # intensité
      log(lambda[pixel,a]) <- beta0[a] +
        beta[1] * x_1[pixel] +
        beta[2] * x_2[pixel] +
        beta[3] * x_3[pixel] + 
        beta[4] * x_4[pixel,a] + 
        beta[5] * x_5[pixel,a] + 
        cell_area[pixel]
      # effort
      logit(b[pixel,a]) <-  alpha0[a] + alpha1 * h_1[pixel,a]
    }
  }
  
  # pour les nobs observations (obs_pixel)
  for (a in 1:nb_annees) {
    obs_denominator[a] <- inprod(lambda[1:npixel,a], b[1:npixel,a]) / nobs[a]
    
    for(obs in 1:nobs[a]){
      ones[obs, a] ~ dbern(
        exp(
          log(lambda[obs_pixel[obs,a],a] * b[obs_pixel[obs,a],a]) -
            obs_denominator[a])   
        / CONSTANT) 
    }
  }
  
  # Priors intercepts
  beta0[1] ~ dnorm(mu.b, sd = 2) 
  
  for (a in 2:nb_annees) {
    beta0[a] ~ dnorm(beta0[a-1], sd = sd.b) 
  }
  
  mu.b ~ dnorm(0, sd = 2)
  sd.b ~ dunif(0, 5)
  
  # Priors pentes
  for(i in 1:5){
    beta[i] ~ ddexp(0, tau) 
  }
  tau ~ dunif(0.001,10)
  
  for (a in 1:nb_annees) {
    alpha0[a] ~ dnorm(mu.a, sd = sd.a)
  }
  mu.a ~ dnorm(0, sd = 2)
  sd.a ~ dunif(0, 5)
  
  alpha1 ~ dnorm(0, sd = 2)
  
})



## Modèle avec tendance linéaire #################################

code3 <- nimbleCode({
  # intensité et effort pour toutes les cellules
  for(pixel in 1:npixel){
    for (a in 1:nb_annees){
      # intensité
      log(lambda[pixel,a]) <- beta0 + c*a +
        beta[1] * x_1[pixel] +
        beta[2] * x_2[pixel] +
        beta[3] * x_3[pixel] + 
        beta[4] * x_4[pixel,a] + 
        beta[5] * x_5[pixel,a] + 
        cell_area[pixel]
      # effort
      logit(b[pixel,a]) <-  alpha0[a] + alpha1 * h_1[pixel,a]
    }
  }
  
  # pour les nobs observations (obs_pixel)
  for (a in 1:nb_annees) {
    obs_denominator[a] <- inprod(lambda[1:npixel,a], b[1:npixel,a]) / nobs[a]
    
    for(obs in 1:nobs[a]){
      ones[obs, a] ~ dbern(
        exp(
          log(lambda[obs_pixel[obs,a],a] * b[obs_pixel[obs,a],a]) -
            obs_denominator[a])   
        / CONSTANT) 
    }
  }
  
  # Priors 
  beta0 ~ dnorm(0, sd = 2) 
  c ~ dnorm(0, sd = 2)
  
  for(i in 1:5){
    beta[i] ~ ddexp(0, tau) 
  }
  tau ~ dunif(0.001,10)
  
  for (a in 1:nb_annees) {
    alpha0[a] ~ dnorm(mu.a, sd = sd.a)
  }
  mu.a ~ dnorm(0, sd = 2)
  sd.a ~ dunif(0, 5)
  
  alpha1 ~ dnorm(0, sd = 2)
  
})



# Estimation -------------------------------------------------------------------

params <- function(modele) {
  if (modele == 0 | modele == 1) {
    return(c("alpha0", "alpha1", "beta0", "beta", "mu.a", "sd.a"))
  }
  if (modele == 2) {
    return(c("alpha0", "alpha1", "beta0", "beta", "mu.a", "mu.b", "sd.a", "sd.b"))
  }
  if (modele == 3) {
    return(c("alpha0", "alpha1", "beta0", "c", "beta", "mu.a", "sd.a"))
  }
}

choix_code <- function(modele) {
  return(get(paste0("code", modele)))
}

#' Estimation des coefficients du modèle sélectionné
#'
#' @param grid # grille à utiliser
#' @param modele # numéro du modèle à utiliser
#' @param periode # vecteur des années voulues
#' @param inits # valeurs initiales à utiliser (déterminées séparément)
#'
estim_param <- function(grid, modele, periode, inits) {
  
  donnees_utiles <- extract_data(grid, periode)
  
  # Paramètres à suivre
  params <- params(modele)
  params2 <- c("lambda", "b")
  
  # MCMC settings
  nc <- 2
  nburn <- 10000 
  ni <- nburn + 30000 
  nt <- 1
  
  # MCMC
  model <- nimbleModel(
    code = choix_code(modele),
    constants = donnees_utiles$constants,
    data = donnees_utiles$data,
    inits = inits
  )
  
  Cmodel <- compileNimble(model)
  
  config <- configureMCMC(Cmodel, 
                          monitors = params,
                          thin = nt,
                          monitors2 = params2,
                          thin2 = 1000, 
                          enableWAIC = TRUE
  )
  
  Rmcmc <- buildMCMC(config)
  Cmcmc <- compileNimble(Rmcmc, project = Cmodel)
  
  out <- runMCMC(
    Cmcmc,
    niter = ni,
    nburnin = nburn,
    nchains = nc,
    setSeed = TRUE,
    # inits = inits_list
    WAIC = TRUE
  )
  
  # sortie
  return(out)
}

# Initiations utilisées --------------------------------------------------------

# Pour le modèle 0
inits0 <- function(){
  list(
    beta0 = rnorm(1, 0, 1),
    beta = rnorm(5, 0, 1), # à mettre ?!
    tau = runif(1, 0.001,10), # à mettre ?!
    alpha0 = rnorm(nb_annees, 0, 1), # à mettre ?!
    mu.a = rnorm(1, 0, 1),
    sd.a = 2,
    alpha1 = rnorm(1, 0, 1)
  )
}

# Pour le modèle 1
inits1 <- function(){
  list(
    beta0 = rnorm(nb_annees, 0, 1),
    beta = rnorm(5, 0, 1), # à mettre ?!
    tau = runif(1, 0.001,10), # à mettre ?!
    alpha0 = rnorm(nb_annees, 0, 1), # à mettre ?!
    mu.a = rnorm(1, 0, 1),
    sd.a = 2,
    alpha1 = rnorm(1, 0, 1)
  )
}

# Pour le modèle 2
inits2 <- function(){
  list(
    beta0 = rnorm(nb_annees, 0, 1),
    mu.b = 0,
    sd.b = 2,
    beta = rnorm(5, 0, 1),
    tau = runif(1, 0.001,10), # à mettre ?!
    alpha0 = rnorm(nb_annees, 0, 1),
    mu.a = rnorm(1, 0, 1),
    sd.a = 2,
    alpha1 = rnorm(1, 0, 1)
  )
}

# Pour le modèle 3
inits3 <- function(){
  list(
    beta0 = rnorm(1, 0, 1),
    c = rnorm(1, 0, 1),
    beta = rnorm(5, 0, 1), # à mettre ?!
    tau = runif(1, 0.001,10), # à mettre ?!
    alpha0 = rnorm(nb_annees, 0, 1), # à mettre ?!
    mu.a = rnorm(1, 0, 1),
    sd.a = 2,
    alpha1 = rnorm(1, 0, 1)
  )
}
