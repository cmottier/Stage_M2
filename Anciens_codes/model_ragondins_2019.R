# Librairies utiles ------------------------------------------------------------

library(nimble)
library(MCMCvis)
library(ggplot2)
library(sf)
# library(stringr)
library(tidyverse)
#
# /!\ Choisir l'effort et le modèle (1 détection ou plusieurs) avant lancement des MCMC
#

# Chargement et travail préalable sur grid -----------------------------------------------

# Chargement de la grille complète
load("RData/grid_sf_5km2.RData")

# on enlève les cellules sans intersection avec Occitanie 
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





# Paramètres communs à tous les modèles ----------------------------------------

# Nombre de pixels
npix <- nrow(grid_selec)

# Aire des pixels
s.area <- as.numeric(units::set_units(grid_selec$area,"km^2"))
logarea <- log(s.area)

# ID des cellules où il y a au moins une occurrence
pixel.id.det <- grid_selec$grid_id[grid_selec$nnutria > 0] 

# Variables environnementales
data <- list(
  cell_area = logarea,
  x_1 = scale(grid_selec$dist_eau)[,1],
  x_2 = scale(grid_selec$logdensity)[,1],
  x_3 = scale(grid_selec$agri_cover)[,1],
  x_4 = scale(grid_selec$prec_cum)[,1],
  x_5 = scale(grid_selec$temp_min)[,1],
  x_6 = scale(grid_selec$temp_max)[,1],
  x_7 = scale(grid_selec$temp_mean)[,1]
  )




# Choix de l'effort ------------------------------------------------------------

data$h_1 <- scale(grid_selec$dist_acces)[,1] # proxy
# data$h_1 <- scale(grid_selec$GBIF)[,1]
# data$h_1 <- scale(grid_selec$logGBIF)[,1]




# Modèle à une détection par cellule -------------------------------------------

## Code ###############

# Bayesian version of the Koshkina (2017) model.
code <- nimbleCode({
  # latent-state model
  for(pixel in 1:npixel){
    # latent state linear predictor
    #
    # x_s  = covariates for latent state
    # beta = latent state model regression coefficients
    # cell_area = log area of grid cell
    #
    
    log(lambda[pixel]) <- beta[1] +
      beta[2] * x_1[pixel] +
      beta[3] * x_2[pixel] +
      beta[4] * x_3[pixel] + 
      beta[5] * x_4[pixel] + 
      beta[6] * x_5[pixel] +
      beta[7] * x_6[pixel] +
      beta[8] * x_7[pixel] +
      cell_area[pixel] # prise en compte de la surface intersectée
    # Species presence in a gridcell as a Bernoulli trial
    # z[pixel] ~ dbern(1 - exp(-lambda[pixel]))
    # presence only thinning prob linear predictor
    #
    # h_s = covariates for thinning probability
    # alpha  = presence-only data model regression coefficients
    #
    logit(b[pixel]) <-  alpha[1] + alpha[2] * h_1[pixel]
  }
  # The presence only data model.
  #
  # This part of the model just uses the
  #  what we have calculated above (lambda
  #  and b). The denominator of this likelihood
  #  is actually a scalar so we can calculate it
  #  outside of a for loop. Let's do that first.
  #
  # The presence_only data model denominator, which
  #  is the thinned poisson process across the
  #  whole region (divided by the total number of
  #  data points because it has to be
  #  evaluated for each data point).
  # m is the number of presence-only data points
  po_denominator <- inprod(lambda[1:npixel], b[1:npixel]) / m
  #
  # Loop through each presence-only data point
  #  using Bernoulli one's trick. The numerator
  #  is just the thinned poisson process for
  #  the ith data point.
  #  po_pixel denotes the grid cell of the ith presence only data point
  for(po in 1:m){
    ones[po] ~ dbern(
      exp(
        log(lambda[po_pixel[po]] * b[po_pixel[po]]) -
          po_denominator)
      / CONSTANT) # attention, voir issue https://github.com/mfidino/integrated-occupancy-model/issues/1
  }
  # Priors for latent state model
  for(i in 1:8){
    beta[i] ~ dnorm(0, sd = 2)
  }
  # Priors for presence-only data model
  for(j in 1:2){
    alpha[j] ~ dnorm(0, sd = 2)
  }
  # Derived parameter, the number of cells occupied
  # zsum <- sum(z[1:npixel])
  # probabilité de présence
  prob[1:npixel] <- 1-exp(-lambda[1:npixel])
})


## data et constantes particulières ##########

data$ones =  rep(1, length(pixel.id.det))

constants <- list(
  npixel = npix,
  m = length(pixel.id.det), 
  CONSTANT = 50000,
  po_pixel = pixel.id.det) 

# Initialisation
inits <- function(){
  list(
    beta = rnorm(8, 0, 1), 
    alpha = rnorm(2, 0, 1)
  )
}

# Paramètres à suivre
params <- c("alpha", "beta") #, "lambda", "b", "prob") # très long si on suit tout !



# Modèle à multiples détections par cellule ------------------------------------

## Code ############

# Bayesian version of the Koshkina (2017) model.
code <- nimbleCode({
  for(pixel in 1:npixel){
    log(lambda[pixel]) <- beta[1] +
      beta[2] * x_1[pixel] +
      beta[3] * x_2[pixel] +
      beta[4] * x_3[pixel] + 
      beta[5] * x_4[pixel] + 
      beta[6] * x_5[pixel] + 
      beta[7] * x_6[pixel] +
      beta[8] * x_7[pixel] +
      cell_area[pixel]

    logit(b[pixel]) <-  alpha[1] + alpha[2] * h_1[pixel]
  }
  
  obs_denominator <- inprod(lambda[1:npixel], b[1:npixel]) / nobs
  
  # on fait ici une boucle sur toutes les observations 
  for(obs in 1:nobs){
    ones[obs] ~ dbern(
      exp(
        log(lambda[po_pixel[obs]] * b[po_pixel[obs]]) -
          obs_denominator)   
      / CONSTANT) 
  }
  
  # Priors for latent state model
  for(i in 1:8){
    beta[i] ~ dnorm(0, sd = 2)
  }
  # Priors for presence-only data model
  for(j in 1:2){
    alpha[j] ~ dnorm(0, sd = 2)
  }
  # Derived parameter, the number of cells occupied
  # zsum <- sum(z[1:npixel])
})



## Data et constantes particulières ########

# nombre d'observations
# troncature pour éviter cellule particulière
grid_selec$nnutria[grid_selec$nnutria > 50] <- 50 # valeur arbitraire à définir

# nombre total d'observations prises en compte
nobs = sum(grid_selec$nnutria)

# pixel associé aux observations (avec répétition)
obs_pixel <- NULL
for (i in 1:length(pixel.id.det)){
  obs_pixel <- c(obs_pixel, rep(pixel.id.det[i], grid_selec$nnutria[pixel.id.det[i]]))
}

constants <- list(
  npixel = npix,
  nobs = nobs,
  CONSTANT = 50000, 
  po_pixel = obs_pixel) 

# data
data$ones <- rep(1, nobs)

# initialisation
inits <- function(){
  list(
    beta = rnorm(8, 0, 1),
    alpha = rnorm(2, 0, 1)
  )
}

# à suivre
params <- c("alpha", "beta")



# MCMC -------------------------------------------------------------------------

## Exécution du modèle #################

# MCMC settings
nc <- 2
nburn <- 8000 
ni <- nburn + 10000 
nt <- 1

# On lance
set.seed(123)
start <- Sys.time()
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
end <- Sys.time()
end - start

# On sauve
# out_prox_mult <- out
# save(out, file = "out_prox_5km2_ens_mult.RData")

# Vérifications
MCMCsummary(out, param=c("alpha", "beta"))
MCMCtrace(out, pdf = FALSE, ind = TRUE, params = c("alpha", "beta"))
MCMCplot(out, params = "beta")

# Résultats
res <- rbind(out$chain1, out$chain2)

## Retour à l'Occitanie ###########

# paramètre à tracer : "lambda", "prob", "^b[^e]"
param = "lambda"

# estimateurs
mask <- str_detect(colnames(res), param)
res_param <- res[,mask]
paramestim <- apply(res_param, 2, median)
parammoy <- apply(res_param, 2, mean)

# plot
p <- ggplot() +
  geom_sf(data = st_intersection(grid_selec, occitanie), color = NA, aes(fill = parammoy)) +
  labs(fill = "Intensité") + # à modifier
  scale_fill_viridis_c(begin = 0, end = 1) +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  # geom_sf(data = nutria) +
  theme_light()
p
# ggsave(plot = p, "Images/map_int_env_5km2_dist.png", dpi = 600)


# # select alpha
# mask <- str_detect(colnames(res), "alpha")
# res_alpha <- res[,mask]
# alphaestim <- apply(res_alpha, 2, median)
# alphamoy <- apply(res_alpha, 2, mean)
# 
# # select beta
# mask <- str_detect(colnames(res), "beta")
# res_beta <- res[,mask]
# betaestim <- apply(res_beta, 2, median)
# betamoy <- apply(res_beta, 2, mean)
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
# # # lambda et b
# # grid_selec$lambda <- exp(betaestim[1] +
# #                            betaestim[2] * data$x_1 +
# #                            betaestim[3] * data$x_2 +
# #                            betaestim[4] * data$x_3 +
# #                            betaestim[5] * data$x_4 +
# #                            betaestim[6] * data$x_5 +
# #                            betaestim[7] * data$x_6 +
# #                            betaestim[8] * data$x_7 +
# #                            + data$cell_area)
# # grid_selec$b <- plogis(alphaestim[1] +
# #                          alphaestim[2] * data$h_1)
# # 
# grid_selec <- grid_selec %>% st_transform(st_crs(occitanie))
# 
# # plot
# p_lambda <- ggplot() +
#   geom_sf(data = st_intersection(grid_selec, occitanie), color = NA, aes(fill = lambdamoy)) +
#   labs(fill = "Intensité") +
#   scale_fill_viridis_c(begin = 0, end = 1) +
#   geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
#   # geom_sf(data = nutria) +
#   theme_light()
# p_lambda
# # ggsave(plot = p_lambda, "Images/map_int_env_5km2_dist.png", dpi = 600)
# 
# p_b <- ggplot() +
#   geom_sf(data = st_intersection(grid_selec, occitanie), color = NA, aes(fill = bmoy)) +
#   labs(fill = "Effort") +
#   scale_fill_viridis_c(begin = 0, end = 1) +
#   geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
#   # geom_sf(data = nutria) +
#   theme_light()
# p_b
# # ggsave(plot = p_b, "Images/map_eff_env_5km2_dist.png", dpi = 600)
# 
# p_p <- ggplot() +
#   geom_sf(data = st_intersection(grid_selec, occitanie), color= NA, aes(fill = probmoy)) +
#   labs(fill = "Présence potentielle estimée du ragondin") +
#   scale_fill_viridis_c(begin = 0, end = 1) +
#   geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
#   # geom_sf(data = nutria) +
#   theme_light()
# p_p
# # ggsave(plot = p_p, "Images/map_pres_env_5km2_dist.png", dpi = 600)
