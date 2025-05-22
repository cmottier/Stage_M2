################################################################################
#                                Simulations                                   #
################################################################################

# Librairies utiles ------------------------------------------------------------
# library(ggplot2)
library(tidyverse)
library(spatstat)
library(nimble)
library(MCMCvis)
# library(plot.matrix)
# library(raster)
library(patchwork)
library(mvtnorm)

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
  
  # probabilité de présence
  p[1:npixel] <- 1-exp(-lambda[1:npixel])
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
  
  # probabilité de présence
  p[1:npixel] <- 1-exp(-lambda[1:npixel])
})

## Paramètres communs #########################

# MCMC settings
nc <- 2
nburn <- 10000 
ni <- nburn + 30000 
nt <- 1
nt2 <- 100

inits <- function() {
  list(beta = rnorm(3, 0, 1), alpha = rnorm(2, 0, 1))
}



# Simulation des observations avec spatstat ------------------------------------

## Avec une loi multinormale (Fithian et al) ############################

simu_data <- function(dim, sigma) {
  # simulation des covariables centrées réduites
  covariables <- rmvnorm(dim^2, sigma = sigma)
  covariables <- apply(covariables, 2, function(x) {
    return((x - mean(x)) / sd(x))
  })
  
  # format qui va bien
  X <- matrix(covariables[, 1], nrow = dim)
  X.im <- as.im(X, square(dim))
  
  X2 <- matrix(covariables[, 2], nrow = dim)
  X2.im <- as.im(X2, square(dim))
  
  W <- matrix(covariables[, 3], nrow = dim)
  W.im <- as.im(W, square(dim))
  
  # aire des cellules
  area <- 1
  logarea <- log(area)
  
  # Objet dat contenant toutes les informations nécessaires
  dat <- list()
  dat$npix <- dim * dim # nombre de pixels
  dat$logarea <- logarea
  dat$s.loc <- as.data.frame(X.im)[, 1:2] 
  dat$xcov <- as.data.frame(X.im)$value
  dat$x2cov <- as.data.frame(X2.im)$value
  dat$wcov <- as.data.frame(W.im)$value
  dat$lambda <- exp(logarea + beta[1] + beta[2] * dat$xcov + beta[3] * dat$x2cov)
  dat$b <- plogis(alpha[1] + alpha[2] * dat$wcov)
  
  # simulation de l'IPPP
  pp <- rpoispp(as.im(matrix(dat$lambda*dat$b, nrow = dim), square(dim)))
  
  # nombre d'observations
  dat$N.det <- pp$n 
  
  # localisation des observations
  dat$loc.det <- data.frame(x = pp$x, y = pp$y) 
  
  # identifiant des cellules des observations
  dat$pixel.id.det <- NULL
  for (i in 1:dat$N.det) {
    dat$pixel.id.det[i] <- which(dat$s.loc[1] == floor(pp$x[i]) + 0.5 &
                                   dat$s.loc[2] == floor(pp$y[i]) + 0.5)
  }

  return(dat)
}

## Avec des fonctions construites ####################################

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

# Valeurs de lambda et b 
lambda_v <- as.im(area*exp(beta[1]+beta[2]*X+beta[3]*X2), square(dim))
b_v <- as.im(plogis(alpha[1]+alpha[2]*W), square(dim))

# simulation de l'IPPP
pp <- rpoispp(b_v*lambda_v)

(dat$N.det <- pp$n)
dat$loc.det <- data.frame(x = pp$x, y = pp$y)

# identifiant des cellules des observations
dat$pixel.id.det <- NULL
for (i in 1:dat$N.det) {
  dat$pixel.id.det[i] <- which(dat$s.loc[1] == floor(pp$x[i]) + 0.5 & dat$s.loc[2] == floor(pp$y[i]) + 0.5)
}

# intensité et effort
dat$lambda <- exp(logarea+beta[1]+beta[2]*dat$xcov+beta[3]*dat$x2cov)
dat$b <- plogis(alpha[1]+alpha[2]*dat$wcov)



# Identifiabilité : matrice de Fisher ------------------------------------------

identif <- function(dat) {
  X <- cbind(1, dat$xcov, dat$x2cov)
  W <- cbind(1, dat$wcov)
  
  Hmat = matrix(nrow=5, ncol=5)
  
  #  beta partials
  for (i in 1:3) {
    for (j in 1:i) {
      Hmat[i,j] = sum(X[,i] * X[,j] * dat$lambda * dat$b)
      Hmat[j,i] = Hmat[i,j]
    }
  }
  
  # alpha partials
  for (i in 1:2) {
    for (j in 1:i) {
      Hmat[3+i, 3+j] = sum(W[,i] * W[,j] * dat$lambda * dat$b * ((1-dat$b)^3) * (1 - exp(2 * W %*% alpha)) ) + sum(W[,i] * W[,j] * dat$b * (1-dat$b))
      Hmat[3+j, 3+i] = Hmat[3+i, 3+j]
    }
  }
  
  # alpha-beta partials
  for (i in 1:2) {
    for (j in 1:3) {
      Hmat[3+i, j] = sum(X[,j] * W[,i] * dat$lambda * dat$b * (1-dat$b))
      Hmat[j, 3+i] = Hmat[3+i, j]
    }
  }
  
  return(min(abs(eigen(Hmat)$values))/ max(abs(eigen(Hmat)$values)))
} 



# Illustrations des données simulées -------------------------------------------

var = "xcov"
var = "x2cov"
var = "wcov"
var = "lambda"
var = "b"

ggplot() +
  geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = dat[[var]]))+
  geom_point(data = dat$loc.det, aes(x = x, y = y), col = "white") +
  scale_fill_viridis_c() +
  labs(x = "", y = "", fill =var)

fill = dat$b*dat$lambda

ggplot() +
  geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = fill))+
  geom_point(data = dat$loc.det, aes(x = x, y = y), col = "white") +
  scale_fill_viridis_c() +
  labs(x = "", y = "", fill = "lambda*b")



# Estimation -------------------------------------------------------------------

estim <- function(dat) {
  data <- list(
    log_area = rep(dat$logarea, dat$npix),
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
      monitors2 <- c("lambda", "b", "p")
    } else {
      code <- code_sans
      monitors <- c("beta")
      monitors2 <- c("lambda", "p")
    }
    
    model <- nimbleModel(
      code = code,
      constants = constants,
      data = data,
      inits = inits()
    )
    
    Cmodel <- compileNimble(model)
    
    conf <- configureMCMC(
      Cmodel,
      monitors = monitors,
      monitors2 = monitors2,
      thin = nt,
      thin2 = nt2
    )
    
    Rmcmc <- buildMCMC(conf)
    Cmcmc <- compileNimble(Rmcmc, project = Cmodel)
    
    out[[type]] <- runMCMC(
        Cmcmc,
        niter = ni,
        nburnin = nburn,
        nchains = nc
        # WAIC = TRUE
      )
  }
  
  return(out)
}


# Extraction des éléments utiles -----------------------------------------------

resume <- function(out) {
  res <- NULL 
  out_a <- rbind(out$avec$samples2$chain1, out$avec$samples2$chain2)
  out_s <- rbind(out$sans$samples2$chain1, out$sans$samples2$chain2)
  
  # lambdas fittés
  mask <- str_detect(colnames(out_a), "lambda")
  res_lambda_a <- out_a[,mask]
  res$lambda_avec  <- apply(res_lambda_a, 2, median)
  
  mask <- str_detect(colnames(out_s), "lambda")
  res_lambda_s <- out_s[,mask]
  res$lambda_sans<- apply(res_lambda_s, 2, median)
  
  # effort fitté
  mask <- str_detect(colnames(out_a), "b[^ed]")
  res_b_a <- out_a[,mask]
  res$b_avec <- apply(res_b_a, 2, median)
  
  # coeff fittés
  res$coeff_avec <- MCMCsummary(out$avec$samples, probs = c(0.5), Rhat = FALSE, n.eff = FALSE)[,3]
  res$coeff_sans <- MCMCsummary(out$sans$samples, probs = c(0.5), Rhat = FALSE, n.eff = FALSE)[,3]
  
  return(res)
}

# # tests
# dim <- 50
# alpha <- c(-3, 1) # effort : b = plogis(...+...*w)
# beta <- c(-2, 1, 1) # intensité : l = exp(...+...*x)
# sigma <- matrix(c(1,0,-0.3,0,1,0,-0.3,0,1), nrow = 3, byrow = TRUE)
# dat <- simu_data(dim, sigma)
# dat$N.det
# identif(dat)
# out <- estim(dat)
# res <- resume(out)

# Qualité des ajustements ------------------------------------------------------

## MCMC ##################

# avec effort
MCMCsummary(out$avec$samples, param=c("alpha", "beta"))
MCMCtrace(out$avec$samples, pdf = FALSE, ind = TRUE, params = c("alpha", "beta"))
MCMCplot(out$avec$samples)

# sans effort
MCMCsummary(out$sans$samples, param=c("beta"))
MCMCtrace(out$sans$samples, pdf = FALSE, ind = TRUE, params = c("beta"))
MCMCplot(out$sans$samples)


## Intensité ################

l_min <- min(dat$lambda, res$lambda_avec, res$lambda_sans)
l_max <- max(dat$lambda, res$lambda_avec, res$lambda_sans)

wrap_plots(
  # Intensité fittée avec effort
  ggplot()+
    geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = res$lambda_avec))+
    # geom_point(data = dat$loc.ipp, aes(x = x, y = y), col = "white") +
    scale_fill_viridis_c(begin = 0, end = 1, limits = c(l_min, l_max)) +
    labs(fill = "l", x = "", y = "", title = "intensité fittée avec effort"),
  
  # Vraie intensité
  ggplot()+
    geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = dat$lambda))+
    # geom_point(data = dat$loc.ipp, aes(x = x, y = y), col = "white") +
    scale_fill_viridis_c(begin = 0, end = 1, limits = c(l_min, l_max)) +
    labs(fill = "l", x = "", y = "", title = "vraie intensité"),
  
  # Intensité fittée sans effort
  ggplot()+
    geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = res$lambda_sans))+
    # geom_point(data = dat$loc.ipp, aes(x = x, y = y), col = "white") +
    scale_fill_viridis_c(begin = 0, end = 1, limits = c(l_min, l_max)) +
    labs(fill = "l", x = "", y = "", title = "intensité fittée sans effort")
)

# Différence entre vraie et avec effort

dl_min <- min(dat$lambda-res$lambda_avec, dat$lambda-res$lambda_sans)
dl_max <- max(dat$lambda-res$lambda_avec, dat$lambda-res$lambda_sans)

wrap_plots(
  ggplot()+
    geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = dat$lambda-res$lambda_avec))+
    # geom_point(data = dat$loc.ipp, aes(x = x, y = y), col = "white") +
    scale_fill_viridis_c(begin = 0, end = 1, limits = c(dl_min, dl_max)) +
    labs(fill = "lv-lf", x = "", y = "", title = "différence d'intensités (vraie - estimée avec effort)"),
  
  # Différence entre vraie et sans effort
  ggplot()+
    geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = dat$lambda-res$lambda_sans))+
    # geom_point(data = dat$loc.ipp, aes(x = x, y = y), col = "white") +
    scale_fill_viridis_c(begin = 0, end = 1, limits = c(dl_min, dl_max)) +
    labs(fill = "lv-lf", x = "", y = "", title = "différence d'intensités (vraie - estimée sans effort)")
)

# Taux d'erreur entre vraie et avec effort
summary(abs(
  (dat$lambda - res$lambda_avec) / 
    dat$lambda
))

# Taux d'erreur entre vraie et sans effort
summary(abs(
  (dat$lambda - res$lambda_sans) / 
    dat$lambda
))

# Erreur quadratique
mean((dat$lambda - res$lambda_avec)^2)
mean((dat$lambda - res$lambda_sans)^2)


## Intensités normalisées ###############

sl_min <- min(scale(dat$lambda)[,1], scale(res$lambda_avec)[,1], scale(res$lambda_sans)[,1])
sl_max <- max(scale(dat$lambda)[,1], scale(res$lambda_avec)[,1], scale(res$lambda_sans)[,1])

wrap_plots(
  # Intensité normalisée fittée avec effort
  ggplot()+
    geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = scale(res$lambda_avec)[,1]))+
    # geom_point(data = dat$loc.ipp, aes(x = x, y = y), col = "white") +
    scale_fill_viridis_c(begin = 0, end = 1, limits = c(sl_min, sl_max)) +
    labs(fill = "l", x = "", y = "", title = "intensité normalisée estimée avec effort"),
  
  # Vraie intensité normalisée
  ggplot()+
    geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = scale(dat$lambda)[,1]))+
    # geom_point(data = dat$loc.ipp, aes(x = x, y = y), col = "white") +
    scale_fill_viridis_c(begin = 0, end = 1, limits = c(sl_min, sl_max)) +
    labs(fill = "l", x = "", y = "", title = "vraie intensité normalisée"),
  
  # Intensité normalisée fittée sans effort
  ggplot()+
    geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = scale(res$lambda_sans)[,1]))+
    # geom_point(data = dat$loc.ipp, aes(x = x, y = y), col = "white") +
    scale_fill_viridis_c(begin = 0, end = 1, limits = c(sl_min, sl_max)) +
    labs(fill = "l", x = "", y = "", title = "intensité normalisée estimée sans effort")
)

dsl_min <- min(scale(dat$lambda)[,1]-scale(res$lambda_avec)[,1], scale(dat$lambda)[,1]-scale(res$lambda_sans)[,1])
dsl_max <- max(scale(dat$lambda)[,1]-scale(res$lambda_avec)[,1], scale(dat$lambda)[,1]-scale(res$lambda_sans)[,1])

wrap_plots(
  # Différence entre vraie et avec effort
  ggplot()+
    geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = scale(dat$lambda)[,1]-scale(res$lambda_avec)[,1]))+
    # geom_point(data = dat$loc.ipp, aes(x = x, y = y), col = "white") +
    scale_fill_viridis_c(begin = 0, end = 1, limits = c(dsl_min, dsl_max)) +
    labs(fill = "lv-lf", x = "", y = "", title = "différence d'intensités normalisées (vraie - avec effort)"),
  
  # Différence entre vraie et sans effort
  ggplot()+
    geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = scale(dat$lambda)[,1]-scale(res$lambda_sans)[,1]))+
    # geom_point(data = dat$loc.ipp, aes(x = x, y = y), col = "white") +
    scale_fill_viridis_c(begin = 0, end = 1, limits = c(dsl_min, dsl_max)) +
    labs(fill = "lv-lf", x = "", y = "", title = "différence d'intensités normalisées (vraie - sans effort)")
)

# Taux d'erreur avec effort
summary(abs(
  (scale(dat$lambda)[, 1] - scale(res$lambda_avec)[, 1]) 
  / scale(dat$lambda)[, 1]
))

# Taux d'erreur sans effort
summary(abs(
  (scale(dat$lambda)[, 1] - scale(res$lambda_sans)[, 1]) 
  / scale(dat$lambda)[, 1]
))

# erreur quadratique
mean((scale(dat$lambda)[, 1] - scale(res$lambda_avec)[, 1])^2)
mean((scale(dat$lambda)[, 1] - scale(res$lambda_sans)[, 1])^2)

## Effort ########################

b_min <- min(dat$b, res$b)
b_max <- max(dat$b, res$b)

wrap_plots(
  # Effort fitté
  ggplot()+
    geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = res$b))+
    # geom_point(data = dat$loc.det, aes(x = x, y = y), col = "white") +
    scale_fill_viridis_c(begin = 0, end = 1, limits = c(b_min, b_max)) +
    labs(fill = "b", x = "", y = "", title = "effort fitté"),
  
  # Vrai effort
  ggplot()+
    geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = dat$b))+
    # geom_point(data = dat$loc.det, aes(x = x, y = y), col = "white") +
    scale_fill_viridis_c(begin = 0, end = 1, limits = c(b_min, b_max)) +
    labs(fill = "b", x = "", y = "", title = "vrai effort")
)

# Différence
ggplot()+
  geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = dat$b-res$b))+
  # geom_point(data = dat$loc.det, aes(x = x, y = y), col = "white") +
  scale_fill_viridis_c(begin = 0, end = 1) +
  labs(fill = "bv-bf", x = "", y = "", title = "différence d'effort")

# Taux d'erreur
summary(abs(
  (dat$b-res$b) / 
    dat$b
))



## Intensité amincie #####################

lb_min <- min(dat$lambda*dat$b, res$lambda_avec*res$b)
lb_max <- max(dat$lambda*dat$b, res$lambda_avec*res$b)

wrap_plots(
  # intensité amincie fittée
  ggplot()+
    geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = res$b*res$lambda_avec))+
    # geom_point(data = dat$loc.det, aes(x = x, y = y), col = "white") +
    scale_fill_viridis_c(begin = 0, end = 1, limits = c(lb_min, lb_max)) +
    labs(fill = "lb", x = "", y = "", title = "intensité amincie fittée"),
  
  # Vraie intensité amincie
  ggplot()+
    geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = dat$b*dat$lambda))+
    # geom_point(data = dat$loc.det, aes(x = x, y = y), col = "white") +
    scale_fill_viridis_c(begin = 0, end = 1, limits = c(lb_min, lb_max)) +
    labs(fill = "lb", x = "", y = "", title = "vraie intensité amincie")
)

# Différence
ggplot()+
  geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = dat$b*dat$lambda-res$b*res$lambda_avec))+
  # geom_point(data = dat$loc.det, aes(x = x, y = y), col = "white") +
  scale_fill_viridis_c(begin = 0, end = 1) +
  labs(fill = "lbv-lbf", x = "", y = "", title = "différence d'intensités amincies")

# Résumé du taux d'erreur
summary(abs(
  (dat$b*dat$lambda-res$b*res$lambda_avec) / 
    (dat$b*dat$lambda)
))



# Evolution de la corrélation --------------------------------------------------

## Construction des données ##################
res <- NULL
res$coeff <- c(alpha, beta)
res$gamma <- NULL
res$lambda <- NULL
res$lambda_avec <- NULL
res$lambda_sans <- NULL
res$b <- NULL
res$b_avec <- NULL
res$coeff_avec <- NULL
res$coeff_sans <- NULL

dim <- 50
sequence <- seq(-1,1,0.1)

for (g in sequence) {
  # simulaton des données
  sigma <- matrix(c(1,0,g,0,1,0,g,0,1), nrow = 3, byrow = TRUE)
  dat <- simu_data(dim, sigma)
  out <- estim(dat)
  out_a <- rbind(out$avec$samples2$chain1, out$avec$samples2$chain2)
  out_s <- rbind(out$sans$samples2$chain1, out$sans$samples2$chain2)
  
  # gamma
  res$gamma <- append(res$gamma, g)
  
  # lambda vrai
  res$lambda <- cbind(res$lambda, dat$lambda)
  
  # lambdas fittés
  mask <- str_detect(colnames(out_a), "lambda")
  res_lambda_a <- out_a[,mask]
  lambda_a <- apply(res_lambda_a, 2, median)
  mask <- str_detect(colnames(out_s), "lambda")
  res_lambda_s <- out_s[,mask]
  lambda_s <- apply(res_lambda_s, 2, median)

  res$lambda_avec <- cbind(res$lambda_avec, lambda_a)
  res$lambda_sans <- cbind(res$lambda_sans, lambda_s)
  
  # effort vrai
  res$b <- cbind(res$b, dat$b)
  
  # effort fitté
  mask <- str_detect(colnames(out_a), "b[^ed]")
  res_b_a <- out_a[,mask]
  b_a <- apply(res_b_a, 2, median)
  
  res$b_avec <- cbind(res$b_avec, b_a)
  
  # coeff fittés
  res$coeff_avec <- cbind(res$coeff_avec, 
                          MCMCsummary(out$avec$samples, probs = c(0.5), Rhat = FALSE, n.eff = FALSE)[,3])
  res$coeff_sans <- cbind(res$coeff_sans, 
                          MCMCsummary(out$sans$samples, probs = c(0.5), Rhat = FALSE, n.eff = FALSE)[,3])

  # sauvegarde de out
  save(out, file = paste0("simu_", as.character(g), ".RData"))
  
  # sauvegarde de res
  save(res, file = "resultat_sequence.RData")
}


## Coefficients fittés #############################

df_avec <- as.data.frame(t(res$coeff_avec))
colnames(df_avec) <- c("alpha_1", "alpha_2", "beta_1", "beta_2", "beta_3")
df_avec <- pivot_longer(df_avec, cols = colnames(df_avec), names_to = "coeff", values_to = "valeur")

df_sans <- as.data.frame(t(res$coeff_sans))
colnames(df_sans) <- c( "beta_1", "beta_2", "beta_3")
df_sans <- pivot_longer(df_sans, cols = colnames(df_sans), names_to = "coeff", values_to = "valeur")

ggplot() +
  geom_line(data = df_avec, aes(x =  rep(res$g, each = 5), y = valeur, color = coeff)) +
  geom_line(data = df_sans, aes(x =  rep(res$g, each = 3), y = valeur, color = coeff), linetype = "dashed") +
  geom_line(aes(x =  res$g, y = res$coeff_avec[1,]+res$coeff_avec[3,]))


## Erreur d'intensité ################

err_int_avec <- apply((res$lambda - res$lambda_avec)^2, 2, mean)
err_int_sans <- apply((res$lambda - res$lambda_sans)^2, 2, mean)

# Taux d'erreur 
ggplot() +
  geom_line(aes(x = res$g[2:21], y =  err_int_avec[2:21])) +# -1 : valeur qui explose !
  geom_line(aes(x = res$g[2:21], y =  err_int_sans[2:21]), col = "red")

# intensités réduites
lambda_red <- apply(res$lambda, 2, function(x) {
  (x-mean(x))/sd(x)
})
lambda_avec_red <- apply(res$lambda_avec, 2, function(x) {
  (x-mean(x))/sd(x)
})
lambda_sans_red <- apply(res$lambda_sans, 2, function(x) {
  (x-mean(x))/sd(x)
})

err_int_red_avec <- apply((lambda_red - lambda_avec_red)^2, 2, mean)
err_int_red_sans <- apply((lambda_red - lambda_sans_red)^2, 2, mean)


ggplot() +
  geom_line(aes(x = res$g, y =  err_int_red_avec)) +# -1 : valeur qui explose !
  geom_line(aes(x = res$g, y =  err_int_red_sans), col = "red")


## Identifiabilité #####################

invcond <- NULL
for (g in sequence) {
  sigma <- matrix(c(1,0,g,0,1,0,g,0,1), nrow = 3, byrow = TRUE)
  val <- NULL
  for (simu in 1:30) {
    dat <- simu_data(dim, sigma)
    val <- append(val, identif(dat))
  }
  invcond <- append(invcond, mean(val))
}

plot(sequence, invcond, type = "l")



