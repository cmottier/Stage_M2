################################################################################
#                                Simulation : scenario 2                                   #
################################################################################

# Dans ce scenario, on utilise des covariables déterministes peu corrélées
# Exemple d'un cas très favorable avec effort (effort élevé)
# effort : b = plogis(-2+2*w)
# intensité : l = exp(-4+4*x1+1*x2)

# Librairies utiles ------------------------------------------------------------
# install.packages("spatstat")
library(spatstat)
library(nimble)
library(raster)

# Codes des modèles avec/sans effort -------------------------------------------

## Modèle avec effort ##########################
code_avec <- nimbleCode({
  # intensité et effort pour toutes les cellules
  for(pixel in 1:npixel){
    # intensité
    log(lambda[pixel]) <- beta[1] +
      beta[2] * x[pixel] +
      beta[3] * x2[pixel] +
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
  beta[3] ~ ddexp(0, tau) 
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
      beta[3] * x2[pixel] +
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
  beta[3] ~ ddexp(0, tau) 
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
    beta = rnorm(3, 0, 1),
    alpha = rnorm(2, 0, 1)
  )
}


# Coefficients choisis ---------------------------------------------------------

alpha <- c(-2, 2) # effort : b = plogis(...+...*w)
beta <- c(-4, 4, 1) # intensité : l = exp(...+...*x)


# Construction des données -----------------------------------------------------

dim <- 50
X.im <- as.im(function(x,y){sqrt(x+y)}, owin(xrange = c(0,dim), yrange = c(0,dim)), dimyx = dim)
X <- as.matrix.im(X.im)
X <- (X-mean(X))/sd(X)
X.im <- as.im(X, square(dim))

X2.im <- as.im(function(x,y){sqrt(100-2*y)}, owin(xrange = c(0,dim), yrange = c(0,dim)), dimyx = dim)
X2 <- as.matrix.im(X2.im)
X2 <- (X2-mean(X2))/sd(X2)
X2.im <- as.im(X2, square(dim))

W.im <- as.im(function(x,y){exp(-(y-dim/2)^2/50)}, owin(xrange = c(0,dim), yrange = c(0,dim)), dimyx = dim)
W <- as.matrix.im(W.im)
W <- (W-mean(W))/sd(W)
W.im <- as.im(W, square(dim))

dat <- list()
dat$npix <- dim*dim
dat$s.loc <- as.data.frame(X.im)[,1:2]
dat$xcov <- as.data.frame(X.im)$value
dat$x2cov <- as.data.frame(X2.im)$value
dat$wcov <- as.data.frame(W.im)$value

# aire des cellules
area <- 1
logarea <- log(area)
dat$logarea <- logarea

# Valeurs de lambda et b 
lambda_v <- as.im(area*exp(beta[1]+beta[2]*X+beta[3]*X2), square(dim))
b_v <- as.im(plogis(alpha[1]+alpha[2]*W), square(dim))



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
    x2 = dat$x2cov,
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

save(resultat, file = "simulation_scenario2.RData")

