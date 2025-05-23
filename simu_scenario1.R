################################################################################
#                                Simulation : scenario 1                                   #
################################################################################

# On utilise les covariables construites par simDataDK
# On considère la même variable pour l'intensité et l'effort : x = w
# effort : b = plogis(-5w)
# intensité : l = exp(6+x)

# Librairies utiles ------------------------------------------------------------
library(spatstat)
library(nimble)
library(raster)

# library(ggplot2)
# library(tidyverse)
# library(sf)
# library(MCMCvis)
# library(plot.matrix)
# library(patchwork)
# library(mvtnorm)
# library(rgl)

# import de la fonction simDataDk telle que dans AHMbook
source("simdatadk_script.R") 


# Codes des modèles avec/sans effort -------------------------------------------

## Modèle avec effort ##########################
code_avec <- nimbleCode({
  # intensité et effort pour toutes les cellules
  for(pixel in 1:npixel){
    # intensité
    log(lambda[pixel]) <- beta[1] +
      beta[2] * x[pixel] +
      log_area[pixel]
    # effort
    logit(b[pixel]) <-  alpha[1] + alpha[2] * w[pixel]
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
  beta[2] ~ ddexp(0, tau) 
  tau ~ dunif(0.001,10)
  
  for(j in 1:2){
    alpha[j] ~ dnorm(0, sd = 2)
  }
  
})


## Modèle sans effort #######################

code_sans <- nimbleCode({
  # intensité pour toutes les cellules
  for(pixel in 1:npixel){
    # intensité
    log(lambda[pixel]) <- beta[1] +
      beta[2] * x[pixel] +
      log_area[pixel]
  }
  
  # pour les nobs observations (obs_pixel)
  obs_denominator <- sum(lambda[1:npixel]) / nobs
  
  for(obs in 1:nobs){
    ones[obs] ~ dbern(
      exp(
        log(lambda[obs_pixel[obs]]) -
          obs_denominator)   
      / CONSTANT) 
  }
  
  # Priors 
  beta[1] ~ dnorm(0, sd = 2) # intercept
  beta[2] ~ ddexp(0, tau) 
  tau ~ dunif(0.001,10)

})

## Paramètres communs #########################

# MCMC settings
nc <- 2
nburn <- 10000 
ni <- nburn + 30000 
nt <- 1
# nt2 <- 100

inits <- function(){
  list(
    beta = rnorm(2, 0, 1),
    alpha = rnorm(2, 0, 1)
  )
}


# Coefficients choisis ---------------------------------------------------------

alpha <- c(0, -5) # effort : b = plogis(...+...*w)
beta <- c(6, 1) # intensité : l = exp(...+...*x)


# Construction des données -----------------------------------------------------

# variable d'environnement x fournie par simDataDK

set.seed(123)

dim <- 50

dat <- simDataDK(
  sqrt.npix = dim, 
  alpha = alpha, 
  beta = beta, 
  drop.out.prop.pb = 0, 
  quadrat.size = 2,
  show.plot = FALSE
) # donne x et w centrées reduites

# surface des cellules
logarea <- log(dat$s.area / dat$npix)

# On change les formats (spatstat) en remplaçant w par x 

X <- matrix(dat$xcov, nrow = dim)
X <- (X-mean(X))/sd(X)
X.im <- as.im(X, square(dim))

W <- X
W.im <- X.im

# Construction de dat 

dat <- list()
dat$npix <- dim*dim
dat$s.loc <- as.data.frame(X.im)[,1:2]
dat$xcov <- as.data.frame(X.im)$value
dat$wcov <- as.data.frame(W.im)$value

# Vraies valeurs de lambda et b 
lambda_v <- as.im(exp(logarea + beta[1]+beta[2]*X), square(dim))
b_v <- as.im(plogis(alpha[1]+alpha[2]*W), square(dim))

dat$lambda <- exp(logarea+beta[1]+beta[2]*dat$xcov)
dat$b <- plogis(alpha[1]+alpha[2]*dat$wcov)


# Réalisation des simulations et estimations -----------------------------------

# simulation des observations (fonction)

simu <- function(dat) {
  # simulation de l'IPPP
  pp <- rpoispp(b_v*lambda_v)
  # nombre d'observations
  dat$N.det <- pp$n
  # localisation des observations
  dat$loc.det <- data.frame(x = pp$x, y = pp$y)
  # identifiant des cellules des observations
  dat$pixel.id.det <- NULL
  for (i in 1:dat$N.det) {
    dat$pixel.id.det[i] <- which(dat$s.loc[1] == floor(pp$x[i]) + 0.5 & dat$s.loc[2] == floor(pp$y[i]) + 0.5)
  }
  return(dat)
}

# exécution des MCMC avec et sans effort (fonction)

estim <- function(dat) {
  data <- list(
    log_area = rep(logarea, dat$npix),
    x = dat$xcov,
    w = dat$wcov,
    ones = rep(1, dat$N.det)
  )
  
  constants <- list(
    npixel = dat$npix,
    nobs = dat$N.det,
    CONSTANT = 50000,
    obs_pixel = dat$pixel.id.det
  )
  
  out <- NULL
  
  for (type in c("avec", "sans")) {
    
    if (type == "avec") {
      code <- code_avec
      monitors <- c("beta", "alpha")
    } else {
      code <- code_sans
      monitors <- c("beta")
    }
    
    out[[type]] <- nimbleMCMC(
      code = code,
      constants = constants,
      data = data,
      inits = inits(),
      monitors = monitors,
      niter = ni,
      nburnin = nburn,
      nchains = nc,
      thin = nt
    )
  }
  
  return(out)
}

# Extraction des éléments utiles 

resume <- function(out) {
  res <- NULL 
  out_a <- rbind(out$avec$chain1, out$avec$chain2)
  out_s <- rbind(out$sans$chain1, out$sans$chain2)
  
  res$avec <- apply(out_a, 2, quantile, probs = c(0.25, 0.5, 0.95))
  res$sans <- apply(out_s, 2, quantile, probs = c(0.25, 0.5, 0.95))
  
  return(res)
}

# Lancement des n simulations

n <- 50
resultat <- NULL 

for (i in 1:n) {
  print(i)
  dat <- simu(dat)
  out <- estim(dat)
  resultat[[i]] <- resume(out)
}

save(resultat, file = "simulation_scenario1.RData")

