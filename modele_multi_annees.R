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

load("grid_sf_50km2.RData")

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

## Fonction pour choisir le modèle  ##################

#' Estimation des paramètres (IPP)
#'
#' @param grid # grille à utiliser
#' @param modele # modèle à utiliser
#' @param periode # période de temps, numérotée à partir de 1
#'
estim_param <- function(grid, modele, periode) {
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
  # params2 <- c("lambda", "b")
  
  # MCMC settings
  nc <- 2
  nburn <- 10000 
  ni <- nburn + 30000 
  nt <- 1
  
  # MCMC
  if (modele == 1) {code <- code1} 
  
  model <- nimbleModel(
    code = code,
    constants = constants,
    data = data,
    inits = inits()
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
  
  # sortie
  return(out)
}


## Lancement et sauvegarde #################

periode = 2010:2012

set.seed(123)

out <- estim_param(
  grid = grid_sf,
  modele = 1,
  periode = periode
)

save(out, file = "modele1_avec_IE_10-12_5km2.RData")


# ## Plot ############################
#
(resume <- MCMCsummary(out$samples, param=c("alpha0", "alpha1", "beta0", "beta")))
MCMCtrace(out$samples, pdf = FALSE, ind = TRUE, params = c("alpha0", "alpha1", "beta0", "beta"))
MCMCplot(out$samples)

library(patchwork)


int_eff <- grid_sf %>%
  select(grid)

# out <- get(paste0("outMCMC_", annee))
# load(paste0("Resultats_MCMC/5km2_avec_lasso/out_multi_gbif_l_iep", annee, ".RData"))
res <- rbind(out$samples2$chain1, out$samples2$chain2)
mask <- str_detect(colnames(res), "lambda")
res_lambda <- res[,mask]
res_lambda <- apply(res_lambda, 2, median)
res_lambda <- matrix(res_lambda, ncol = length(periode), byrow = FALSE)
for (a in 1:length(periode)) {
  int_eff[[paste0("lambda_med", periode[a]) ]] <- res_lambda[,a]
}
mask <- str_detect(colnames(res), "b[^d]")
res_b <- res[,mask]
res_b <- apply(res_b, 2, median)
res_b <- matrix(res_b, ncol = length(periode), byrow = FALSE)
for (a in 1:length(periode)) {
  int_eff[[paste0("b_med", periode[a]) ]] <- res_b[,a]
}
#
#
#
# Intensité
plot_l <- do.call(
  wrap_plots,
  lapply(periode,
         function(x) {
           ggplot() +
             geom_sf(data = int_eff, color = NA, aes(fill = get(paste0("lambda_med", x)))) +
             labs(fill = "Intensité") +
             scale_fill_viridis_c(begin = 0, end = 1, limits = c(0, 10)) + # pour fixer couleurs
           labs(title = x) +
             theme_light()
         })
)

plot_l

# effort
plot_b <- do.call(
  wrap_plots,
  lapply(periode,
         function(x) {
           ggplot() +
             geom_sf(data = int_eff, color = NA, aes(fill = get(paste0("b_med", x)))) +
             labs(fill = "Effort") +
             scale_fill_viridis_c(begin = 0, end = 1, limits = c(0, 1)) + # pour fixer couleurs
             labs(title = x) +
             theme_light()
         })
)

plot_b

#' #' # Code Nimble -----------------------------------------------------------------
#' #' # Multiples détections par cellule
#' #' 
#' #' # avec un effet aléatoire sur intercepts seuls
#' #' 
#' #' code1 <- nimbleCode({
#' #'   # intensité et effort pour toutes les cellules
#' #'   for(pixel in 1:npixel){
#' #'     for (a in 1:nb_annees){
#' #'       # intensité
#' #'       log(lambda[pixel,a]) <- beta0[a] +
#' #'         beta[1] * x_1[pixel] +
#' #'         beta[2] * x_2[pixel] +
#' #'         beta[3] * x_3[pixel] + 
#' #'         beta[4] * x_4[pixel,a] + 
#' #'         beta[5] * x_5[pixel,a] + 
#' #'         cell_area[pixel]
#' #'       # effort
#' #'       logit(b[pixel,a]) <-  alpha0[a] + alpha1 * h_1[pixel,a]
#' #'     }
#' #'   }
#' #'   
#' #'   # pour les nobs observations (obs_pixel)
#' #'   for (a in 1:nb_annees) {
#' #'     obs_denominator[a] <- inprod(lambda[1:npixel,a], b[1:npixel,a]) / nobs[a]
#' #'     
#' #'     for(obs in 1:nobs[a]){
#' #'       ones[[a]][obs] ~ dbern(
#' #'         exp(
#' #'           log(lambda[obs_pixel[[a]][obs]] * b[obs_pixel[[a]][obs]]) -
#' #'             obs_denominator[a])   
#' #'         / CONSTANT) 
#' #'     }
#' #'   }
#' #'   
#' #'   # Priors 
#' #'   for (a in 1:nb_annees) {
#' #'     beta0[a] ~ dnorm(0, sd = 2) 
#' #'     alpha0[a] ~ dnorm(0, sd = 2)
#' #'     }
#' #'   for(i in 1:5){
#' #'     beta[i] ~ ddexp(0, tau) 
#' #'   }
#' #'   tau ~ dunif(0.001,10)
#' #'   alpha1 ~ dnorm(0, sd = 2)
#' #' 
#' #'   # probabilité de présence
#' #'   # p[1:npixel,periode] <- 1-exp(-lambda[1:npixel,periode])
#' #' })
#' #' 
#' #' 
#' #' 
#' #' # Estimation -------------------------------------------------------------------
#' #' 
#' #' ## Fonction pour choisir le modèle  ##################
#' #' 
#' #' #' Estimation des paramètres (IPP)
#' #' #'
#' #' #' @param grid # grille à utiliser
#' #' #' @param modele # modèle à utiliser
#' #' #' @param periode # période de temps, numérotée à partir de 1
#' #' #'
#' #' estim_param <- function(grid, modele, periode) {
#' #'   # Nombre d'années
#' #'   nb_annees <- length(periode)
#' #'   
#' #'   # Nombre de pixels
#' #'   npix <- nrow(grid)
#' #'   
#' #'   # Aire des pixels
#' #'   s.area <- as.numeric(units::set_units(grid$area,"km^2"))
#' #'   logarea <- log(s.area)
#' #'   
#' #'   # Observations annuelles
#' #'   pixel.id.det <- NULL
#' #'   nb_observations <- NULL
#' #'   nobs <- NULL
#' #'   ones <- NULL
#' #'   obs_pixel <- NULL
#' #'   
#' #'   for (a in 1:nb_annees) {
#' #'     # ID des cellules où il y a au moins une occurrence
#' #'     pixel.id.det[[a]] <- grid$grid_id[grid[[paste0("nnutria", periode[a])]] > 0]
#' #'     
#' #'     # Troncature des observations
#' #'     nb_observations[[a]] <- grid[[paste0("nnutria", periode[a])]] 
#' #'     nb_observations[[a]][nb_observations[[a]] > 50] <- 50 # valeur arbitraire à définir
#' #'     
#' #'     # nombre total d'observations prises en compte
#' #'     nobs[a] <- sum(nb_observations[[a]])
#' #'     
#' #'     # pixel associé aux observations (avec répétition)
#' #'     obs_pixel[[a]] <- rep(NA, nobs[a])
#' #'     s <- 0
#' #'     for (i in (1:length(pixel.id.det[[a]]))){
#' #'       obs_pixel[[a]][(s+1):(s + nb_observations[[a]][pixel.id.det[[a]][i]])] <- pixel.id.det[[a]][i]
#' #'       s <- s + nb_observations[[a]][pixel.id.det[[a]][i]]
#' #'     }
#' #'     
#' #'     # Liste de vecteur de 1
#' #'     ones[[a]] <- rep(1, nobs[a])
#' #'   }
#' #'   
#' #'   # Variables 
#' #'   data <- list(cell_area = logarea,
#' #'                x_1 = scale(grid$dist_eau)[,1],
#' #'                x_2 = scale(grid$logdensity)[,1],
#' #'                x_3 = scale(grid$agri_cover)[,1],
#' #'                x_4 = apply(grid_sf %>% select(paste0("pcum_" , periode)) %>% st_drop_geometry, 2, function(x){return((x-mean(x))/sd(x))}),
#' #'                x_5 = apply(grid_sf %>% select(paste0("tmin_" , periode)) %>% st_drop_geometry, 2, function(x){return((x-mean(x))/sd(x))}),
#' #'                h_1 = apply(grid_sf %>% select(paste0("dgbif_" , periode)) %>% st_drop_geometry, 2, function(x){return((x-mean(x))/sd(x))}),
#' #'                ones = ones
#' #'                )
#' #'   
#' #'   constants <- list(
#' #'     npixel = npix,
#' #'     nobs = nobs,
#' #'     CONSTANT = 50000,
#' #'     obs_pixel = obs_pixel
#' #'     ) 
#' #'   
#' #'   # Initialisation
#' #'   inits <- function(){
#' #'     list(
#' #'       beta0 = rnorm(nb_annees, 0, 1),
#' #'       beta = rnorm(5, 0, 1), 
#' #'       alpha0 = rnorm(1, 0, 1),
#' #'       alpha1 = rnorm(1, 0, 1),
#' #'     )
#' #'   }
#' #'   
#' #'   # Paramètres à suivre
#' #'   params <- c("alpha0", "alpha1", "beta0", "beta") 
#' #'   # params2 <- c("lambda", "b") 
#' #'   
#' #'   # MCMC settings
#' #'   nc <- 2
#' #'   nburn <- 10000 
#' #'   ni <- nburn + 30000 
#' #'   nt <- 1
#' #'   
#' #'   # MCMC
#' #'   if (modele == 1) {code <- code1} 
#' #'   
#' #'   model <- nimbleModel(
#' #'     code = code,
#' #'     constants = constants,
#' #'     data = data,
#' #'     inits = inits()
#' #'   )
#' #'   
#' #'   Cmodel <- compileNimble(model)
#' #'   
#' #'   conf <- configureMCMC(Cmodel, 
#' #'                         monitors = params,
#' #'                         thin = nt
#' #'                         # monitors2 = params2, 
#' #'                         # thin2 = 1000
#' #'                         )
#' #'   
#' #'   Rmcmc <- buildMCMC(conf)
#' #'   Cmcmc <- compileNimble(Rmcmc, project = Cmodel)
#' #'   
#' #'   out <- runMCMC(
#' #'     Cmcmc,
#' #'     niter = ni,
#' #'     nburnin = nburn,
#' #'     nchains = nc
#' #'     # WAIC = TRUE
#' #'   )
#' #'   
#' #'   # sortie
#' #'   return(out)
#' #' }
#' #' 
#' #' 
#' #' ## Lancement et sauvegarde #################
#' #' 
#' #' periode = 2010:2012
#' #' 
#' #' out <- estim_param(
#' #'     grid = grid_sf,
#' #'     modele = 1,
#' #'     periode = periode
#' #' )
#' #' 
