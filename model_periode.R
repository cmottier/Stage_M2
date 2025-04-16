# Librairies utiles ------------------------------------------------------------

library(nimble)
library(MCMCvis)
library(ggplot2)
library(sf)
library(tidyverse)


# Chargement et travail préalable sur grid -------------------------------------

# Chargement de la grille complète
# load("RData/grid_sf_5km2_periode.RData")

# on enlève si nécessaire les cellules sans intersection avec Occitanie 
grid_selec <- grid_sf %>%
  filter(as.numeric(area)!=0)

# # on garde les cellules qui intersectent à au moins 10%
# grid_selec <- grid_selec %>%
#   filter(as.numeric(area)>0.1*max(as.numeric(area)))

# on renumérote les cellules de la sélection
grid_selec$grid_id = 1:lengths(grid_selec)[1] 

# gestion des NA et -Inf de grid_sf (à paramétrer)
summary(grid_selec$logdensity)
grid_selec$logdensity[is.infinite(grid_selec$logdensity)] <- -50 # valeur artificielle à déterminer 
grid_selec$agri_cover[is.na(grid_selec$agri_cover)] <- 0


# Codes Nimble -----------------------------------------------------------------

# Une détection par cellule
code_uni <- nimbleCode({
  # pour toutes les cellules
  for(pixel in 1:npixel){
    # intensité
    log(lambda[pixel]) <- beta[1] +
      beta[2] * x_1[pixel] +
      beta[3] * x_2[pixel] +
      beta[4] * x_3[pixel] + 
      beta[5] * x_4[pixel] + 
      beta[6] * x_5[pixel] +
      beta[7] * x_6[pixel] +
      beta[8] * x_7[pixel] +
      cell_area[pixel] 
    # effort
    logit(b[pixel]) <-  alpha[1] + alpha[2] * h_1[pixel]
  }
  
  # pour les m cellules contenant une détection
  po_denominator <- inprod(lambda[1:npixel], b[1:npixel]) / m
  
  for(po in 1:m){
    ones[po] ~ dbern(
      exp(
        log(lambda[po_pixel[po]] * b[po_pixel[po]]) -
          po_denominator)
      / CONSTANT) 
  }
  
  # Priors 
  for(i in 1:8){
    beta[i] ~ dnorm(0, sd = 2)
  }

  for(j in 1:2){
    alpha[j] ~ dnorm(0, sd = 2)
  }
})


# Multiples détections par cellule
code_multi <- nimbleCode({
  # pour toutes les cellules
  for(pixel in 1:npixel){
    # intensité
    log(lambda[pixel]) <- beta[1] +
      beta[2] * x_1[pixel] +
      beta[3] * x_2[pixel] +
      beta[4] * x_3[pixel] + 
      beta[5] * x_4[pixel] + 
      beta[6] * x_5[pixel] + 
      beta[7] * x_6[pixel] +
      beta[8] * x_7[pixel] +
      cell_area[pixel]
    # effort
    logit(b[pixel]) <-  alpha[1] + alpha[2] * h_1[pixel]
  }
  
  # pour les nobs observations
  obs_denominator <- inprod(lambda[1:npixel], b[1:npixel]) / nobs
  
  for(obs in 1:nobs){
    ones[obs] ~ dbern(
      exp(
        log(lambda[obs_pixel[obs]] * b[obs_pixel[obs]]) -
          obs_denominator)   
      / CONSTANT) 
  }
  
  # Priors 
  for(i in 1:8){
    beta[i] ~ dnorm(0, sd = 2)
  }

  for(j in 1:2){
    alpha[j] ~ dnorm(0, sd = 2)
  }
})



# Estimation -------------------------------------------------------------------

#' Estimation des paramètres (IPP)
#'
#' @param grid # grille à utiliser
#' @param modele # une détection (1) ou multiples détections (2)
#' @param effort # proximité aux routes ('prox') ou GBIF ('gbif')
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
               x_1 = scale(grid_selec$dist_eau)[,1],
               x_2 = scale(grid_selec$logdensity)[,1],
               x_3 = scale(grid_selec$agri_cover)[,1],
               x_4 = scale(grid_selec[[paste0("pcum_", annee)]])[,1],
               x_5 = scale(grid_selec[[paste0("tmin_", annee)]])[,1],
               x_6 = scale(grid_selec[[paste0("tmax_", annee)]])[,1],
               x_7 = scale(grid_selec[[paste0("tmean_", annee)]])[,1])
  
  if (effort == 'prox') {
    data$h_1 <- scale(grid_selec$dist_acces)[, 1]
  }
  else {
    data$h_1 <- scale(grid_selec[[paste0("dgbif_", annee)]])[, 1]
  }
  
  if (modele == 1) { data$ones <- rep(1, length(pixel.id.det)) }
  else {
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
      beta = rnorm(8, 0, 1), 
      alpha = rnorm(2, 0, 1)
    )
  }
  
  # Paramètres à suivre
  params <- c("alpha", "beta") 
  
  # MCMC settings
  nc <- 2
  nburn <- 8000 
  ni <- nburn + 10000 
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

set.seed(123)
out_2010 <- estim_param(grid = grid_selec, modele = 1, effort = "prox", annee = 2010)
MCMCsummary(out_2010, param=c("alpha", "beta"))


# On sauve
# out_prox_mult <- out
# save(out, file = "out_2010.RData")

# Vérifications
MCMCsummary(out, param=c("alpha", "beta"))
MCMCtrace(out, pdf = FALSE, ind = TRUE, params = c("alpha", "beta"))
MCMCplot(out, params = "beta")

# Résultats
res <- rbind(out$chain1, out$chain2)

## Retour à l'Occitanie ###########

# # paramètre à tracer : "lambda", "prob", "^b[^e]"
# param = "lambda"
# 
# # estimateurs
# mask <- str_detect(colnames(res), param)
# res_param <- res[,mask]
# paramestim <- apply(res_param, 2, median)
# parammoy <- apply(res_param, 2, mean)
# 
# # plot
# p <- ggplot() +
#   geom_sf(data = st_intersection(grid_selec, occitanie), color = NA, aes(fill = parammoy)) +
#   labs(fill = "Intensité") + # à modifier
#   scale_fill_viridis_c(begin = 0, end = 1) +
#   geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
#   # geom_sf(data = nutria) +
#   theme_light()
# p
# # ggsave(plot = p, "Images/map_int_env_5km2_dist.png", dpi = 600)


# select alpha
mask <- str_detect(colnames(res), "alpha")
res_alpha <- res[,mask]
alphaestim <- apply(res_alpha, 2, median)
alphamoy <- apply(res_alpha, 2, mean)

# select beta
mask <- str_detect(colnames(res), "beta")
res_beta <- res[,mask]
betaestim <- apply(res_beta, 2, median)
betamoy <- apply(res_beta, 2, mean)
# 
# # select lambda
# mask <- str_detect(colnames(res), "lambda")
# res_lambda <- res[,mask]
# lambdaestim <- apply(res_lambda, 2, median)
# lambdamoy <- apply(res_lambda, 2, mean)
# 
# # select b
# mask <- str_detect(colnames(res), "^b[^e]") # pour ne pas sélectionner beta et lambda...
# res_b <- res[,mask]
# bestim <- apply(res_b, 2, median)
# bmoy <- apply(res_b, 2, mean)
# 
# # select prob
# mask <- str_detect(colnames(res), "prob")
# res_prob <- res[,mask]
# probestim <- apply(res_prob, 2, median)
# probmoy <- apply(res_prob, 2, mean)
# 
# lambda et b
grid_selec$lambda <- exp(betaestim[1] +
                           betaestim[2] * data$x_1 +
                           betaestim[3] * data$x_2 +
                           betaestim[4] * data$x_3 +
                           betaestim[5] * data$x_4 +
                           betaestim[6] * data$x_5 +
                           betaestim[7] * data$x_6 +
                           betaestim[8] * data$x_7 +
                           + data$cell_area)
grid_selec$b <- plogis(alphaestim[1] +
                         alphaestim[2] * data$h_1)

# grid_selec <- grid_selec %>% st_transform(st_crs(occitanie))
# 
# plot
p_lambda <- ggplot() +
  geom_sf(data = st_intersection(grid_selec, occitanie), color = NA, aes(fill = lambda)) +
  labs(fill = "Intensité") +
  scale_fill_viridis_c(begin = 0, end = 1) +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  # geom_sf(data = nutria) +
  theme_light()
p_lambda
# ggsave(plot = p_lambda, "Images/map_int_env_5km2_dist.png", dpi = 600)

p_b <- ggplot() +
  geom_sf(data = st_intersection(grid_selec, occitanie), color = NA, aes(fill = b)) +
  labs(fill = "Effort") +
  scale_fill_viridis_c(begin = 0, end = 1) +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  # geom_sf(data = nutria) +
  theme_light()
p_b
# ggsave(plot = p_b, "Images/map_eff_env_5km2_dist.png", dpi = 600)

p_p <- ggplot() +
  geom_sf(data = st_intersection(grid_selec, occitanie), color= NA, aes(fill = 1-exp(-lambda))) +
  labs(fill = "Présence potentielle estimée du ragondin") +
  scale_fill_viridis_c(begin = 0, end = 1) +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  # geom_sf(data = nutria) +
  theme_light()
p_p
# # ggsave(plot = p_p, "Images/map_pres_env_5km2_dist.png", dpi = 600)
