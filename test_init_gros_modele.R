################################################################################
#                           Modèle multi-année                                 #
################################################################################

# Librairies utiles ------------------------------------------------------------

library(nimble)
library(sf)
library(tidyverse)
# library(MCMCvis)
# library(patchwork)



# Chargement  de la grille -----------------------------------------------------

load("RData/50km2/grid_sf_50km2.RData")

# Code Nimble -----------------------------------------------------------------
# Multiples détections par cellule

# avec un effet aléatoire sur intercepts seuls

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
    alpha0[a] ~ dnorm(0, sd = 2)
  }
  for(i in 1:5){
    beta[i] ~ ddexp(0, tau) 
  }
  tau ~ dunif(0.001,10)
  alpha1 ~ dnorm(0, sd = 2)
  
  # probabilité de présence
  # p[1:npixel,periode] <- 1-exp(-lambda[1:npixel,periode])
})



# Estimation -------------------------------------------------------------------

grid <- grid_sf
periode <- 2010:2024

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


# Variables 
data <- list(cell_area = logarea,
             x_1 = scale(grid$dist_eau)[,1],
             x_2 = scale(grid$logdensity)[,1],
             x_3 = scale(grid$agri_cover)[,1],
             x_4 = apply(grid_sf %>% select(paste0("pcum_" , periode)) %>% st_drop_geometry, 2, function(x){return((x-mean(x))/sd(x))}),
             x_5 = apply(grid_sf %>% select(paste0("tmin_" , periode)) %>% st_drop_geometry, 2, function(x){return((x-mean(x))/sd(x))}),
             h_1 = apply(grid_sf %>% select(paste0("dgbif_" , periode)) %>% st_drop_geometry, 2, function(x){return((x-mean(x))/sd(x))}),
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

# Initialisation
inits <- function(){
  list(
    beta0 = rnorm(nb_annees, 0, 1),
    beta = rnorm(5, 0, 1), 
    alpha0 = rnorm(nb_annees, 0, 1),
    alpha1 = rnorm(1, 0, 1),
    tau = runif(1, 0.001,10) # à mettre ?!
  )
}

# Paramètres à suivre
params <- c("alpha0", "alpha1", "beta0", "beta") 
params2 <- c("lambda", "b")

# MCMC settings
nc <- 2
nburn <- 10000 
ni <- nburn + 30000 
nt <- 1


model <- nimbleModel(
  code = code1,
  constants = constants,
  data = data,
  inits = inits()
)

model$initializeInfo()

comp <- 0
res <- -Inf
while ((comp < 1000) & is.infinite(res)) {
  print(comp)
  comp <- comp + 1
  model$simulate()
  res <- model$calculate()
}

inits <- list(
    beta0 = model$beta0,
    beta = model$beta, 
    alpha0 = model$alpha0,
    alpha1 = model$alpha1,
    tau = model$tau # à mettre ?!
  )

save(inits, file = "RData/Inits/init_50km2_2010-2024")

model <- nimbleModel(
  code = code1,
  constants = constants,
  data = data,
  inits = inits
)

Cmodel <- compileNimble(model)

conf <- configureMCMC(Cmodel, 
                      monitors = params,
                      thin = nt,
                      monitors2 = params2,
                      thin2 = 1000
)

Rmcmc <- buildMCMC(conf)
Cmcmc <- compileNimble(Rmcmc, project = Cmodel)

out <- runMCMC(
  Cmcmc,
  niter = ni,
  nburnin = nburn,
  nchains = nc,
  setSeed = TRUE
  # WAIC = TRUE
)
  

save(out, file = "Resultats_MCMC/Gros_modele/50km2/modele1_avec_IE_10-24_50km2.RData")
