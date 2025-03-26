library(nimble)
library(MCMCvis)
library(ggplot2)
library(sf)
library(stringr)
library(tidyverse)

load("grid_sf.RData")

###################################
# avec 1 détection par cellule...
###################################

# gestion des NA et -Inf de grid_sf (à vérifier...)
grid_sf$surface_en_eau[is.na(grid_sf$surface_en_eau)] <- 0
grid_sf$logdensity[is.na(grid_sf$logdensity)] <- 0
grid_sf$logdensity[is.infinite(grid_sf$logdensity)] <- -100 # valeur artificielle à déterminer ? 
grid_sf$agri_cover[is.na(grid_sf$agri_cover)] <- 0

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
      beta[9] * x_8[pixel] + cell_area
    # Species presence in a gridcell as a Bernoulli trial
    z[pixel] ~ dbern(1 - exp(-lambda[pixel]))
    # presence only thinning prob linear predictor
    #
    # h_s = covariates for thinning probability
    # alpha  = presence-only data model regression coefficients
    #
    logit(b[pixel]) <-  alpha[1] + alpha[2] * h_1[pixel] + alpha[3] * h_2[pixel]
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
  for(i in 1:9){
    beta[i] ~ dnorm(0, sd = 2)
  }
  # Priors for presence-only data model
  for(j in 1:3){
    alpha[j] ~ dnorm(0, sd = 2)
  }
  # Derived parameter, the number of cells occupied
  # zsum <- sum(z[1:npixel])
})

head(grid_sf$grid_id) # les ID des cellules

pixel.id.det <- grid_sf$grid_id[grid_sf$nnutria > 0] # les ID des cellules où il y a au moins une occurrence
head(pixel.id.det)

npix <- nrow(grid_sf)
s.area <- as.numeric(units::set_units(st_area(grid_sf)[1],"km^2"))
(logarea <- log(s.area / npix)) # erreur dans le code d'origine ? -> s.area déjà par cellule...
# logarea <- log(s.area) # pose problème pour logProb...

data <- list(
  cell_area = logarea,
  x_1 = scale(grid_sf$surface_en_eau)[,1],
  x_2 = scale(grid_sf$logdensity)[,1],
  x_3 = scale(grid_sf$agri_cover)[,1],
  x_4 = scale(grid_sf$temp_min)[,1],
  x_5 = scale(grid_sf$temp_max)[,1],
  x_6 = scale(grid_sf$temp_mean)[,1],
  x_7 = scale(grid_sf$prec_cum)[,1],
  x_8 = scale(grid_sf$lgr_rivieres)[,1],
  h_1 = scale(grid_sf$lgr_chemins)[,1],
  h_2 = scale(grid_sf$lgr_routes)[,1],
  ones = rep(1, length(pixel.id.det)))

constants <- list(
  npixel = npix,
  m = length(pixel.id.det), 
  CONSTANT = 50000,
  po_pixel = pixel.id.det) 

zinit <- numeric(npix)
zinit[pixel.id.det] <- 1
inits <- function(){
  list(
    beta = rnorm(9, 0, 1), 
    alpha = rnorm(3, 0, 1),
    z = zinit # pourquoi en commentaire dans le code initial ?
  )
}

params <- c("alpha", "beta", "z")

# MCMC settings (pour tester...)
nc <- 2
nburn <- 5000 #5000
ni <- nburn + 10000 #30000
nt <- 1

# run the model!
set.seed(10) # graine pour voir apparaître -Inf dans logProb_ones
set.seed(123) # graine pour ne pas voir apparaître -Inf dans logProb_ones
model <-  nimbleModel(
  code = code,
  constants = constants,
  data = data,
  inits = inits())
model$logProb_ones 
# si on corrige logarea, tous les essais ont donné -Inf et le MCMC ne marche pas...

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
  thin = nt)
end <- Sys.time()
end - start

MCMCsummary(out)

MCMCtrace(out, pdf = FALSE, ind = TRUE, params = "alpha")
MCMCplot(out, params = "beta")

res <- rbind(out$chain1, out$chain2)

# select z
mask <- str_detect(colnames(res), "z")
res_z <- res[,mask]
grid_sf$zestim <- apply(res_z, 2, median)
grid_sf$zmoy <- apply(res_z, 2, mean)

# viz
ggplot() +
  geom_sf(data = grid_sf, lwd = 0.1, aes(fill = as_factor(zestim))) +
  labs(fill = "Présence potentielle estimée du ragondin") +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  # geom_sf(data = nutria) +
  theme_void()

ggplot() +
  geom_sf(data = grid_sf, lwd = 0.1, aes(fill = zmoy)) +
  labs(fill = "Présence potentielle estimée du ragondin") +
  scale_fill_viridis_c() +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  # geom_sf(data = nutria) +
  theme_void()


###################################
# avec multiples détections par cellules -> ne marche pas... 
###################################

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
      beta[9] * x_8[pixel] + cell_area
    # Species presence in a gridcell as a Bernoulli trial
    z[pixel] ~ dbern(1 - exp(-lambda[pixel]))
    # presence only thinning prob linear predictor
    #
    # h_s = covariates for thinning probability
    # alpha  = presence-only data model regression coefficients
    #
    logit(b[pixel]) <-  alpha[1] + alpha[2] * h_1[pixel] + alpha[3] * h_2[pixel]
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
        nobs_pixel[po]*log(lambda[po_pixel[po]] * b[po_pixel[po]]) -
          po_denominator)   # modif ici pour prendre en compte le nombre d'observations par pixel : nobs_pixel[po]
      / CONSTANT) # attention, voir issue https://github.com/mfidino/integrated-occupancy-model/issues/1
  }
  # Priors for latent state model
  for(i in 1:9){
    beta[i] ~ dnorm(0, sd = 2)
  }
  # Priors for presence-only data model
  for(j in 1:3){
    alpha[j] ~ dnorm(0, sd = 2)
  }
  # Derived parameter, the number of cells occupied
  # zsum <- sum(z[1:npixel])
})

head(grid_sf$grid_id) # ID des cellules

pixel.id.det <- grid_sf$grid_id[grid_sf$nnutria > 0] # les ID des cellules où il y a au moins une occurrence
head(pixel.id.det)

npix <- nrow(grid_sf)
s.area <- as.numeric(units::set_units(st_area(grid_sf)[1],"km^2")) # si les cellules sont d'aires identiques
logarea <- log(s.area / npix)

data <- list(
  cell_area = logarea,
  x_1 = scale(grid_sf$surface_en_eau)[,1],
  x_2 = scale(grid_sf$logdensity)[,1],
  x_3 = scale(grid_sf$agri_cover)[,1],
  x_4 = scale(grid_sf$temp_min)[,1],
  x_5 = scale(grid_sf$temp_max)[,1],
  x_6 = scale(grid_sf$temp_mean)[,1],
  x_7 = scale(grid_sf$prec_cum)[,1],
  x_8 = scale(grid_sf$lgr_rivieres)[,1],
  h_1 = scale(grid_sf$lgr_chemins)[,1],
  h_2 = scale(grid_sf$lgr_routes)[,1],
  ones = rep(1, length(pixel.id.det)))

constants <- list(
  npixel = npix,
  m = length(pixel.id.det),
  CONSTANT = 1000000, # problème : la constante à choisir pour guarantir proba ! 
  po_pixel = pixel.id.det,
  nobs_pixel = grid_sf$nnutria[grid_sf$nnutria > 0]) # modifié pour prendre en compte toutes les détections

zinit <- numeric(npix)
zinit[pixel.id.det] <- 1
inits <- function(){
  list(
    beta = rnorm(9, 0, 1),
    alpha = rnorm(3, 0, 1),
    z = zinit # pourquoi mis en commentaire dans code initial ?
  )
}

params <- c("alpha", "beta", "z")

# MCMC settings (pour tester...)
nc <- 2
nburn <- 5000 #5000
ni <- nburn + 10000 #30000
nt <- 1

# pour comprendre le problème du modèle, il faut comprendre ça : 
model <- nimbleModel(code, constants = constants, data = data, inits = inits())
model$logProb_ones

# run the model!
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
  thin = nt)
end <- Sys.time()
end - start

MCMCsummary(out)

res <- rbind(out$chain1, out$chain2)


# # select z
# mask <- str_detect(colnames(res), "z")
# res_z <- res[,mask]
# grid_sf$zestim <- apply(res_z, 2, median)
# grid_sf$zmoy <- apply(res_z, 2, mean)
# 
# # viz
# ggplot() +
#   geom_sf(data = grid_sf, lwd = 0.1, aes(fill = as_factor(zestim))) +
#   labs(fill = "Présence potentielle estimée du ragondin") +
#   geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
#   geom_sf(data = nutria) +
#   theme_void()
# 
# ggplot() +
#   geom_sf(data = grid_sf, lwd = 0.1, aes(fill = zmoy)) +
#   labs(fill = "Présence potentielle estimée du ragondin") +
#   scale_fill_viridis_c() +
#   geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
#   geom_sf(data = nutria) +
#   theme_void()
