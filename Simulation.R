################################################################################
#                                Simulations                                   #
################################################################################

# Librairies utiles ------------------------------------------------------------
library(ggplot2)
library(tidyverse)
# library(sf)
library(spatstat)
library(AHMbook)
library(nimble)
library(MCMCvis)
library(plot.matrix)
library(raster)
library(patchwork)

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

# Coefficients choisis ---------------------------------------------------------

alpha <- c(-2, 1) # effort : b = plogis(...+...*w)
beta <- c(6, 1) # intensité : l = exp(...+...*x)

# Simulation des observations --------------------------------------------------

## Avec AHMbook #######################################

# On utilise ici la variable d'environnement x et la variable 
# d'effort w déterminées dans l'annexe de Koshkina

set.seed(123)

dat <- simDataDK(
  sqrt.npix = 50, 
  alpha = alpha, 
  beta = beta, 
  drop.out.prop.pb = 0, 
  quadrat.size = 2,
  # gamma = c(0, -1.5), 
  # nquadrats = 250, 
  # nsurveys = 0,
  show.plot = FALSE
) # donne x et w centrées reduites

str(dat, 1)

# Our landscape is inhabited by a total of ... individuals (dat$N.ipp)
dat$N.ipp

# of which ... are detected
dat$N.det

#Their locations represent our presence-only data: 
# loc.det contains the coordinates 
head(dat$loc.det) # Actual coordinates of points detected

# and pixel.id.det the pixel ID for each detected point.
head(dat$pixel.id.det) # Pixel ID of 302 points detected

# le thinning est donc de :
1-dat$N.det/dat$N.ipp

# surface des cellules
(logarea <- log(dat$s.area / dat$npix))

# Plot
ggplot() +
  geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = dat$xcov))+
  geom_point(data = dat$loc.ipp, aes(x = x, y = y), col = "white") +
  scale_fill_viridis_c() +
  labs(fill = "x", x = "", y = "")

# valeur de la variable d'effort h et observations
ggplot()+
  geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = dat$wcov))+
  geom_point(data = dat$loc.det, aes(x = x, y = y), col = "white") +
  scale_fill_viridis_c() +
  labs(fill = "w", x = "", y = "")

# ou...
# plot(dat$s)

## Avec le package spatstat ########################################

### Covariables aléatoires ############################
dim <- 50
X <- matrix(runif(dim*dim, -1, 1), nrow = dim)
X <- (X-mean(X))/sd(X)
X.im <- as.im(X, square(dim))
# plot(X.im)

W <- matrix(rbeta(dim*dim, 1, 3), nrow = dim)
W <- (W-mean(W))/sd(W)
W.im <- as.im(W, square(dim))
# plot(W.im)

### Covariables construites ##########################

X.im <- as.im(function(x,y){sqrt(x+y)}, owin(xrange = c(0,dim), yrange = c(0,dim)), dimyx = dim)
X <- as.matrix.im(X.im)
X <- (X-mean(X))/sd(X)
X.im <- as.im(X, square(dim))

# W.im <- as.im(function(x,y){sqrt(dim-x+dim-y)}, owin(xrange = c(0,dim), yrange = c(0,dim)), dimyx = dim)
W.im <- as.im(function(x,y){exp(-(y-dim/2)^2/50)}, owin(xrange = c(0,dim), yrange = c(0,dim)), dimyx = dim)
W <- as.matrix.im(W.im)
W <- (W-mean(W))/sd(W)
W.im <- as.im(W, square(dim))


### Construction de dat ##########################

dat <- list()
dat$npix <- dim*dim
dat$s.loc <- as.data.frame(X.im)[,1:2]
dat$xcov <- as.data.frame(X.im)$value
dat$wcov <- as.data.frame(W.im)$value

# aire des cellules
area <- 1
logarea <- log(area)

# Valeurs de lambda et b 
lambda_v <- as.im(area*exp(beta[1]+beta[2]*X), square(dim))
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


## Intensités et efforts ###############################

dat$lambda <- exp(logarea+beta[1]+beta[2]*dat$xcov)
dat$b <- plogis(alpha[1]+alpha[2]*dat$wcov)

# Illustration des données simulées --------------------------------------------

# choix de la variable : 
var = "xcov"
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

data <- list(
  log_area = rep(logarea, dat$npix),
  x = scale(dat$xcov)[, 1],
  w = scale(dat$wcov)[, 1],
  ones = rep(1, dat$N.det)
)

constants <- list(
  npixel = dat$npix,
  nobs = dat$N.det,
  CONSTANT = 50000,
  obs_pixel = dat$pixel.id.det
) 

inits <- function(){
  list(
    beta = rnorm(2, 0, 1),
    alpha = rnorm(2, 0, 1)
  )
}

## Avec effort ########################################

model_avec <- nimbleModel(
  code = code_avec,
  constants = constants,
  data = data,
  inits = inits()
)

Cmodel_avec <- compileNimble(model_avec)

conf_avec <- configureMCMC(Cmodel_avec, monitors = c("beta", "alpha"),
                        monitors2 = c("lambda", "b", "p"), thin = nt, thin2 = nt2)

Rmcmc_avec <- buildMCMC(conf_avec)
Cmcmc_avec <- compileNimble(Rmcmc_avec, project = Cmodel_avec)

out <- runMCMC(
  Cmcmc_avec,
  niter = ni,
  nburnin = nburn,
  nchains = nc
  # WAIC = TRUE
)

## Sans effort ########################################

model_s <- nimbleModel(
  code = code_sans,
  constants = constants,
  data = data,
  inits = inits()
)

Cmodel_s <- compileNimble(model_s)

conf_s <- configureMCMC(Cmodel_s, monitors = c("beta"),
                        monitors2 = c("lambda", "p"), thin = nt, thin2 = nt2)

Rmcmc_s <- buildMCMC(conf_s)
Cmcmc_s <- compileNimble(Rmcmc_s, project = Cmodel_s)

out_s <- runMCMC(
  Cmcmc_s,
  niter = ni,
  nburnin = nburn,
  nchains = nc
  # WAIC = TRUE
)


# Qualité des ajustements ------------------------------------------------------

## MCMC ##################

alpha
beta

# avec effort
MCMCsummary(out$samples, param=c("alpha", "beta"))
MCMCtrace(out$samples, pdf = FALSE, ind = TRUE, params = c("alpha", "beta"))
MCMCplot(out$samples)

# sans effort
MCMCsummary(out_s$samples, param=c("beta"))
MCMCtrace(out_s$samples, pdf = FALSE, ind = TRUE, params = c("beta"))
MCMCplot(out_s$samples)


## Intensité ################

res <- rbind(out$samples2$chain1, out$samples2$chain2)
res_s <- rbind(out_s$samples2$chain1, out_s$samples2$chain2)

# extraction des intensités (médiane)

# avec effort
mask <- str_detect(colnames(res), "lambda")
res_lambda <- res[,mask]
lambdaestim <- apply(res_lambda, 2, median)
lambdasd <- apply(res_lambda, 2, sd)

# sans effort
mask <- str_detect(colnames(res_s), "lambda")
res_s_lambda <- res_s[,mask]
lambdaestim_s <- apply(res_s_lambda, 2, median)
lambdasd_s <- apply(res_s_lambda, 2, sd)

l_min <- min(dat$lambda, lambdaestim, lambdaestim_s)
l_max <- max(dat$lambda, lambdaestim, lambdaestim_s)

wrap_plots(
  # Intensité fittée avec effort
  ggplot()+
    geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = lambdaestim))+
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
    geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = lambdaestim_s))+
    # geom_point(data = dat$loc.ipp, aes(x = x, y = y), col = "white") +
    scale_fill_viridis_c(begin = 0, end = 1, limits = c(l_min, l_max)) +
    labs(fill = "l", x = "", y = "", title = "intensité fittée sans effort")
)

# Différence entre vraie et avec effort

dl_min <- min(dat$lambda-lambdaestim, dat$lambda-lambdaestim_s)
dl_max <- max(dat$lambda-lambdaestim, dat$lambda-lambdaestim_s)

wrap_plots(
  ggplot()+
    geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = dat$lambda-lambdaestim))+
    # geom_point(data = dat$loc.ipp, aes(x = x, y = y), col = "white") +
    scale_fill_viridis_c(begin = 0, end = 1, limits = c(dl_min, dl_max)) +
    labs(fill = "lv-lf", x = "", y = "", title = "différence d'intensités (vraie - estimée avec effort)"),
  
  # Différence entre vraie et sans effort
  ggplot()+
    geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = dat$lambda-lambdaestim_s))+
    # geom_point(data = dat$loc.ipp, aes(x = x, y = y), col = "white") +
    scale_fill_viridis_c(begin = 0, end = 1, limits = c(dl_min, dl_max)) +
    labs(fill = "lv-lf", x = "", y = "", title = "différence d'intensités (vraie - estimée sans effort)")
)

# Taux d'erreur entre vraie et avec effort
summary(abs(
  (dat$lambda - lambdaestim) / 
    dat$lambda
))

# Taux d'erreur entre vraie et sans effort
summary(abs(
  (dat$lambda - lambdaestim_s) / 
    dat$lambda
))


## Intensités normalisées ###############

sl_min <- min(scale(dat$lambda)[,1], scale(lambdaestim)[,1], scale(lambdaestim_s)[,1])
sl_max <- max(scale(dat$lambda)[,1], scale(lambdaestim)[,1], scale(lambdaestim_s)[,1])

wrap_plots(
  # Intensité normalisée fittée avec effort
  ggplot()+
    geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = scale(lambdaestim)[,1]))+
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
    geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = scale(lambdaestim_s)[,1]))+
    # geom_point(data = dat$loc.ipp, aes(x = x, y = y), col = "white") +
    scale_fill_viridis_c(begin = 0, end = 1, limits = c(sl_min, sl_max)) +
    labs(fill = "l", x = "", y = "", title = "intensité normalisée estimée sans effort")
)

dsl_min <- min(scale(dat$lambda)[,1]-scale(lambdaestim)[,1], scale(dat$lambda)[,1]-scale(lambdaestim_s)[,1])
dsl_max <- max(scale(dat$lambda)[,1]-scale(lambdaestim)[,1], scale(dat$lambda)[,1]-scale(lambdaestim_s)[,1])

wrap_plots(
  # Différence entre vraie et avec effort
  ggplot()+
    geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = scale(dat$lambda)[,1]-scale(lambdaestim)[,1]))+
    # geom_point(data = dat$loc.ipp, aes(x = x, y = y), col = "white") +
    scale_fill_viridis_c(begin = 0, end = 1, limits = c(dsl_min, dsl_max)) +
    labs(fill = "lv-lf", x = "", y = "", title = "différence d'intensités normalisées (vraie - avec effort)"),
  
  # Différence entre vraie et sans effort
  ggplot()+
    geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = scale(dat$lambda)[,1]-scale(lambdaestim_s)[,1]))+
    # geom_point(data = dat$loc.ipp, aes(x = x, y = y), col = "white") +
    scale_fill_viridis_c(begin = 0, end = 1, limits = c(dsl_min, dsl_max)) +
    labs(fill = "lv-lf", x = "", y = "", title = "différence d'intensités normalisées (vraie - sans effort)")
)

# Taux d'erreur avec effort
summary(abs(
  (scale(dat$lambda)[, 1] - scale(lambdaestim)[, 1]) 
  / scale(dat$lambda)[, 1]
))

# Taux d'erreur sans effort
summary(abs(
  (scale(dat$lambda)[, 1] - scale(lambdaestim_s)[, 1]) 
  / scale(dat$lambda)[, 1]
))


## Effort ########################

# extraction des efforts (médiane)
mask <- str_detect(colnames(res), "b[^ed]")
res_b <- res[,mask]
bestim <- apply(res_b, 2, median)
bsd <- apply(res_b, 2, sd)

b_min <- min(dat$b, bestim)
b_max <- max(dat$b, bestim)

wrap_plots(
  # Effort fitté
  ggplot()+
    geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = bestim))+
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
  geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = dat$b-bestim))+
  # geom_point(data = dat$loc.det, aes(x = x, y = y), col = "white") +
  scale_fill_viridis_c(begin = 0, end = 1) +
  labs(fill = "bv-bf", x = "", y = "", title = "différence d'effort")

# Taux d'erreur
summary(abs(
  (dat$b-bestim) / 
    dat$b
))


## Effort normalisé ################

sb_min <- min(scale(dat$b)[,1], scale(bestim)[,1])
sb_max <- max(scale(dat$b)[,1], scale(bestim)[,1])

# Effort fitté
ggplot()+
  geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = scale(bestim)[,1]))+
  # geom_point(data = dat$loc.det, aes(x = x, y = y), col = "white") +
  scale_fill_viridis_c(begin = 0, end = 1, limits = c(sb_min, sb_max)) +
  labs(fill = "b", x = "", y = "", title = "effort fitté, normalisé")

# Vrai effort
ggplot()+
  geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = scale(dat$b)[,1]))+
  # geom_point(data = dat$loc.det, aes(x = x, y = y), col = "white") +
  scale_fill_viridis_c(begin = 0, end = 1, limits = c(sb_min, sb_max)) +
  labs(fill = "b", x = "", y = "", title = "vrai effort, normalisé" )

# Différence
ggplot()+
  geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = scale(dat$b)[,1]-scale(bestim)[,1]))+
  # geom_point(data = dat$loc.ipp, aes(x = x, y = y), col = "white") +
  scale_fill_viridis_c(begin = 0, end = 1) +
  labs(fill = "", x = "", y = "", title = "différence d'efforts normalisés")

# Résumé du taux d'erreur
summary(abs(
  (scale(dat$b)[,1]-scale(bestim)[,1]) 
  / scale(dat$b)[,1]
))


## Intensité amincie #####################

lb_min <- min(dat$lambda*dat$b, lambdaestim*bestim)
lb_max <- max(dat$lambda*dat$b, lambdaestim*bestim)

wrap_plots(
  # intensité amincie fittée
  ggplot()+
    geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = bestim*lambdaestim))+
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
  geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = dat$b*dat$lambda-bestim*lambdaestim))+
  # geom_point(data = dat$loc.det, aes(x = x, y = y), col = "white") +
  scale_fill_viridis_c(begin = 0, end = 1) +
  labs(fill = "lbv-lbf", x = "", y = "", title = "différence d'intensités amincies")

# Résumé du taux d'erreur
summary(abs(
  (dat$b*dat$lambda-bestim*lambdaestim) / 
    (dat$b*dat$lambda)
))


# Test matrice de Fisher -------------------------------------------------------

alpha_0 <- seq(-5, 5, by = 0.05)
alpha_1 <- -1 
beta <- c(log(8000), 0.5) 

dat <- simDataDK(
  sqrt.npix = 50, 
  alpha = c(0, alpha_1), 
  beta = beta, 
  drop.out.prop.pb = 0, 
  quadrat.size = 2,
  show.plot = FALSE
) # donne x et w centrées reduites

# surface des cellules
logarea <- log(dat$s.area / dat$npix)

rcon <- NULL
verif <- NULL
for (a in alpha_0) { 
  
  dat$lambda <- exp(logarea+beta[1]+beta[2]*dat$xcov)
  dat$b <- plogis(a+alpha_1*dat$wcov)
  
  
  partial_beta <- matrix(
    data = c(
      sum(dat$lambda*dat$b),
      sum(dat$xcov*dat$lambda*dat$b),
      sum(dat$xcov*dat$lambda*dat$b),
      sum(dat$xcov^2*dat$lambda*dat$b)
    ),
    nrow = 2, 
    byrow = TRUE
  )
  
  partial_beta_alpha <- matrix(
    data = c(
      sum(dat$lambda*dat$b*(1-dat$b)),
      sum(dat$wcov*dat$lambda*dat$b*(1-dat$b)),
      sum(dat$xcov*dat$lambda*dat$b*(1-dat$b)),
      sum(dat$xcov*dat$wcov*dat$lambda*dat$b*(1-dat$b))
    ),
    nrow = 2, 
    byrow = TRUE
  )
  
  partial_alpha <- matrix(
    data = c(
      sum(dat$lambda*dat$b*(1-dat$b)^3*(1-exp(2*logit(dat$b))) + dat$lambda*dat$b^2*(1-dat$b)),
      sum(dat$wcov*dat$lambda*dat$b*(1-dat$b)^3*(1-exp(2*logit(dat$b))) + dat$wcov*dat$lambda*dat$b^2*(1-dat$b)),
      sum(dat$wcov*dat$lambda*dat$b*(1-dat$b)^3*(1-exp(2*logit(dat$b))) + dat$wcov*dat$lambda*dat$b^2*(1-dat$b)),
      sum(dat$wcov^2*dat$lambda*dat$b*(1-dat$b)^3*(1-exp(2*logit(dat$b))) + dat$wcov^2*dat$lambda*dat$b^2*(1-dat$b))
    ),
    nrow = 2, 
    byrow = TRUE
  )
  
  I <- rbind(cbind(partial_beta, partial_beta_alpha),
             cbind(t(partial_beta_alpha), partial_alpha))
  
  verif <- cbind(verif, sum(eigen(I)$values >0) == 4)
  rcon <- cbind(rcon, eigen(I)$values[4]/ eigen(I)$values[1])
}

table(verif)
lines(alpha_0, rcon, type = "l")
