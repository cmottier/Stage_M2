################################################################################
#               Extraction de la population - Meso@LR                          #
################################################################################

# /!\ à lancer sur le cluster

library(tidyverse)
library(sf)
# library(ggplot2)

# # contours des départements d'Occitanie
# dpts_occitanie <- st_read("Data/departements-d-occitanie.shp") 
# 
# # contours de la région
# occitanie <- dpts_occitanie %>% st_union()
# 
# rm(dpts_occitanie)
# 
# # Box
# occitanie_box <- st_bbox(c(
#   xmin = -0.4,
#   xmax = 5,
#   ymax = 45.1,
#   ymin = 42
# ), crs = st_crs(4326)) %>%
#   st_as_sfc() %>%
#   st_transform(crs = st_crs(occitanie))
# 
# # Extrcation de la population utile
# pop <- st_read("Data/grid_1km_surf.gpkg") %>%
#   st_transform(crs = st_crs(occitanie)) %>%
#   st_crop(occitanie_box) %>%
#   st_intersection(occitanie)
# 
# # sauvegarde de la sortie
# save(pop, file = "pop.RData")

# # plot
# ggplot() +
#   geom_sf(data = pop, aes(fill = log(TOT_P_2021))) +
#   scale_fill_viridis_c()


load("pop.RData")
load("Data/grid_sf_5km2.RData")


# # Comparaison
# print("calcul de difference")
# Diff <- st_sym_difference(agri12, agri18)
# 
# save(Diff, file = "Diff.RData")

# Avec la grille
grid_pop <- pop %>%
  st_transform(crs = st_crs(grid_sf)) %>%
  select(TOT_P_2011, TOT_P_2018, TOT_P_2021, geometry) %>%
  st_intersection(grid_sf) %>%
  group_by(grid_id) %>%
  summarise(hab_2021 = sum(TOT_P_2021), hab_2011 = sum(TOT_P_2011), hab_2018 = sum(TOT_P_2018)) %>%
  as_tibble() %>%
  select(-geometry)

  
# Densité 
grid_pop_evol <- grid_sf %>% 
  select(grid_id, area, grid) %>%
  full_join(grid_pop, by = "grid_id") %>%
  mutate(logdensity_2021 = log(as.numeric(hab_2021/area) + 10^(-20))) %>% 
  mutate(logdensity_2011 = log(as.numeric(hab_2011/area) + 10^(-20))) %>% 
  mutate(logdensity_2018 = log(as.numeric(hab_2018/area) + 10^(-20)))


save(grid_pop_evol, file = "grid_pop.RData")

# Etude des différences de la grille
grid_pop_evol <- grid_pop_evol %>%
  mutate(diff_11_21 = abs(grid_pop_evol$logdensity_2021 - grid_pop_evol$logdensity_2011))

sum(grid_pop_evol$diff_11_21 != 0) # la plupart des cellules contient un changement
summary(grid_pop_evol$diff_11_21)



ggplot() +
  geom_sf(data = grid_pop_evol, aes(fill = diff_11_21)) +
  scale_fill_viridis_c() +
  theme_void()

ggplot() +
  geom_sf(data = grid_pop_evol, aes(fill = as.numeric(diff_11_21>1))) +
  scale_fill_viridis_c() +
  theme_void()


# Différences entre les modèles obtenus

library(nimble)
library(MCMCvis)

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
  for(i in 1:6){
    beta[i] ~ dnorm(0, sd = 2)
  }
  
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
               x_2 = scale(grid_pop_evol$logdensity_2021)[,1],
               x_3 = scale(grid$agri_cover)[,1],
               x_4 = scale(grid[[paste0("pcum_", annee)]])[,1],
               x_5 = scale(grid[[paste0("tmin_", annee)]])[,1]
  )
  
  if (effort == 'prox') {
    data$h_1 <- scale(grid$dist_acces)[, 1]
  }
  else {
    data$h_1 <- scale(grid[[paste0("dgbif_", annee)]])[, 1]
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

set.seed(123)
out_pop_2011 <- estim_param(grid_sf, 1, "prox", 2018)
set.seed(123)
out_pop_2018 <- estim_param(grid_sf, 1, "prox", 2018)
set.seed(123)
out_pop_2021 <- estim_param(grid_sf, 1, "prox", 2018)

MCMCsummary(out_pop_2011)      
MCMCsummary(out_pop_2018)
MCMCsummary(out_pop_2021)





