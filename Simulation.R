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




# Avec AHMbook -----------------------------------------------------------------

## Simulation d'un PPP aminci ###########################

# On utilise ici la variable d'environnement x et la variable 
# d'effort w déterminées dans l'annexe de Koshkina

set.seed(123)

alpha <- c(-2, -5) # effort : b = plogis(-2+1*w)
beta <- c(4, 1) # intensité : l = exp(6+1*x)

dat <- simDataDK(
  sqrt.npix = 100, 
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

ggplot()+
  geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = dat$wcov))+
  geom_point(data = dat$loc.det, aes(x = x, y = y), col = "white") +
  scale_fill_viridis_c() +
  labs(fill = "h", x = "", y = "")

## Ajustement du modèle avec effort ##########################

# Code (Multiples détections par cellule)

code_multi <- nimbleCode({
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

  
# Variables 
data <- list(log_area = rep(logarea, dat$npix) ,
             x = scale(dat$xcov)[,1],
             w = scale(dat$wcov)[,1],
             ones = rep(1, dat$N.det)
)
  
constants <- list(
  npixel = dat$npix,
  nobs = dat$N.det,
  CONSTANT = 50000,
  obs_pixel = dat$pixel.id.det
) 
  
# Initialisation
inits <- function(){
  list(
    beta = rnorm(2, 0, 1), 
    alpha = rnorm(2, 0, 1)
  )
}
  
# Paramètres à suivre
params <- c("alpha", "beta") 
params2 <- c("lambda", "b", "p") 
  
# MCMC settings
nc <- 2
nburn <- 10000 
ni <- nburn + 30000 
nt <- 1
  
# MCMC
model <- nimbleModel(
  code = code_multi,
  constants = constants,
  data = data,
  inits = inits()
)
  
Cmodel <- compileNimble(model)
  
conf <- configureMCMC(Cmodel, monitors = params,
                        monitors2 = params2, thin = nt, thin2 = 100)
  
Rmcmc <- buildMCMC(conf)
Cmcmc <- compileNimble(Rmcmc, project = Cmodel)
  
out <- runMCMC(
  Cmcmc,
  niter = ni,
  nburnin = nburn,
  nchains = nc
  # WAIC = TRUE
)


## Ajustement du modèle sans effort #####################

# Code (Multiples détections par cellule)

code_sans_eff <- nimbleCode({
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


# Variables 
data_sans_eff <- list(log_area = rep(logarea, dat$npix) ,
                      x = scale(dat$xcov)[,1],
                      ones = rep(1, dat$N.det)
)


# Initialisation
inits_sans_eff <- function(){
  list(
    beta = rnorm(2, 0, 1)
  )
}

# MCMC
model_s <- nimbleModel(
  code = code_sans_eff,
  constants = constants,
  data = data_sans_eff,
  inits = inits_sans_eff()
)

Cmodel_s <- compileNimble(model_s)

conf_s <- configureMCMC(Cmodel_s, monitors = c("beta"),
                        monitors2 = c("lambda", "p"), thin = nt, thin2 = 100)

Rmcmc_s <- buildMCMC(conf_s)
Cmcmc_s <- compileNimble(Rmcmc_s, project = Cmodel_s)

out_s <- runMCMC(
  Cmcmc_s,
  niter = ni,
  nburnin = nburn,
  nchains = nc
  # WAIC = TRUE
)



## Qualité des ajustements ##########

### MCMC ##################

# avec effort
MCMCsummary(out$samples, param=c("alpha", "beta"))
MCMCtrace(out$samples, pdf = FALSE, ind = TRUE, params = c("alpha", "beta"))
MCMCplot(out$samples)

# sans effort
MCMCsummary(out_s$samples, param=c("beta"))
MCMCtrace(out_s$samples, pdf = FALSE, ind = TRUE, params = c("beta"))
MCMCplot(out_s$samples)


### Intensité ################

res <- rbind(out$samples2$chain1, out$samples2$chain2)
res_s <- rbind(out_s$samples2$chain1, out_s$samples2$chain2)

# extraction des intensités (médiane)

lambdavraie <- exp(logarea+beta[1]+beta[2]*dat$xcov)

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

# Intensité fittée avec effort

l_min <- min(lambdavraie, lambdaestim, lambdaestim_s)
l_max <- max(lambdavraie, lambdaestim, lambdaestim_s)

ggplot()+
  geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = lambdaestim))+
  # geom_point(data = dat$loc.ipp, aes(x = x, y = y), col = "white") +
  scale_fill_viridis_c(begin = 0, end = 1, limits = c(l_min, l_max)) +
  labs(fill = "l", x = "", y = "", title = "intensité fittée avec effort")

# Vraie intensité
ggplot()+
  geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = lambdavraie))+
  # geom_point(data = dat$loc.ipp, aes(x = x, y = y), col = "white") +
  scale_fill_viridis_c(begin = 0, end = 1, limits = c(l_min, l_max)) +
  labs(fill = "l", x = "", y = "", title = "vraie intensité")

# Intensité fittée sans effort
ggplot()+
  geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = lambdaestim_s))+
  # geom_point(data = dat$loc.ipp, aes(x = x, y = y), col = "white") +
  scale_fill_viridis_c(begin = 0, end = 1, limits = c(l_min, l_max)) +
  labs(fill = "l", x = "", y = "", title = "intensité fittée sans effort")

# Différence entre vraie et avec effort

dl_min <- min(lambdavraie-lambdaestim, lambdavraie-lambdaestim_s)
dl_max <- max(lambdavraie-lambdaestim, lambdavraie-lambdaestim_s)

ggplot()+
  geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = lambdavraie-lambdaestim))+
  # geom_point(data = dat$loc.ipp, aes(x = x, y = y), col = "white") +
  scale_fill_viridis_c(begin = 0, end = 1, limits = c(dl_min, dl_max)) +
  labs(fill = "lv-lf", x = "", y = "", title = "différence d'intensités (vraie - estimée avec effort)")

# Différence entre vraie et sans effort
ggplot()+
  geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = lambdavraie-lambdaestim_s))+
  # geom_point(data = dat$loc.ipp, aes(x = x, y = y), col = "white") +
  scale_fill_viridis_c(begin = 0, end = 1, limits = c(dl_min, dl_max)) +
  labs(fill = "lv-lf", x = "", y = "", title = "différence d'intensités (vraie - estimée sans effort)")

# Taux d'erreur entre vraie et avec effort
summary(abs(
  (lambdavraie - lambdaestim) / 
    lambdavraie
))

# Taux d'erreur entre vraie et sans effort
summary(abs(
  (lambdavraie - lambdaestim_s) / 
    lambdavraie
))

### Intensités normalisées ###############

sl_min <- min(scale(lambdavraie)[,1], scale(lambdaestim)[,1], scale(lambdaestim_s)[,1])
sl_max <- max(scale(lambdavraie)[,1], scale(lambdaestim)[,1], scale(lambdaestim_s)[,1])

# Intensité normalisée fittée avec effort
ggplot()+
  geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = scale(lambdaestim)[,1]))+
  # geom_point(data = dat$loc.ipp, aes(x = x, y = y), col = "white") +
  scale_fill_viridis_c(begin = 0, end = 1, limits = c(sl_min, sl_max)) +
  labs(fill = "l", x = "", y = "", title = "intensité normalisée estimée avec effort")

# Vraie intensité normalisée
ggplot()+
  geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = scale(lambdavraie)[,1]))+
  # geom_point(data = dat$loc.ipp, aes(x = x, y = y), col = "white") +
  scale_fill_viridis_c(begin = 0, end = 1, limits = c(sl_min, sl_max)) +
  labs(fill = "l", x = "", y = "", title = "vraie intensité normalisée")

# Intensité normalisée fittée sans effort
ggplot()+
  geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = scale(lambdaestim_s)[,1]))+
  # geom_point(data = dat$loc.ipp, aes(x = x, y = y), col = "white") +
  scale_fill_viridis_c(begin = 0, end = 1, limits = c(sl_min, sl_max)) +
  labs(fill = "l", x = "", y = "", title = "intensité normalisée estimée sans effort")

dsl_min <- min(scale(lambdavraie)[,1]-scale(lambdaestim)[,1], scale(lambdavraie)[,1]-scale(lambdaestim_s)[,1])
dsl_max <- max(scale(lambdavraie)[,1]-scale(lambdaestim)[,1], scale(lambdavraie)[,1]-scale(lambdaestim_s)[,1])

# Différence entre vraie et avec effort
ggplot()+
  geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = scale(lambdavraie)[,1]-scale(lambdaestim)[,1]))+
  # geom_point(data = dat$loc.ipp, aes(x = x, y = y), col = "white") +
  scale_fill_viridis_c(begin = 0, end = 1, limits = c(dsl_min, dsl_max)) +
  labs(fill = "lv-lf", x = "", y = "", title = "différence d'intensités normalisées (vraie - avec effort)")

# Différence entre vraie et sans effort
ggplot()+
  geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = scale(lambdavraie)[,1]-scale(lambdaestim_s)[,1]))+
  # geom_point(data = dat$loc.ipp, aes(x = x, y = y), col = "white") +
  scale_fill_viridis_c(begin = 0, end = 1, limits = c(dsl_min, dsl_max)) +
  labs(fill = "lv-lf", x = "", y = "", title = "différence d'intensités normalisées (vraie - sans effort)")

# Taux d'erreur avec effort
summary(abs(
  (scale(lambdavraie)[, 1] - scale(lambdaestim)[, 1]) 
  / scale(lambdavraie)[, 1]
))

# Taux d'erreur sans effort
summary(abs(
  (scale(lambdavraie)[, 1] - scale(lambdaestim_s)[, 1]) 
  / scale(lambdavraie)[, 1]
))


### Effort ########################

# extraction des efforts (médiane)
mask <- str_detect(colnames(res), "b[^ed]")
res_b <- res[,mask]
bestim <- apply(res_b, 2, median)
bsd <- apply(res_b, 2, sd)

bvrai <- plogis(alpha[1]+alpha[2]*dat$wcov)

b_min <- min(bvrai, bestim)
b_max <- max(bvrai, bestim)

# Effort fitté
ggplot()+
  geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = bestim))+
  # geom_point(data = dat$loc.det, aes(x = x, y = y), col = "white") +
  scale_fill_viridis_c(begin = 0, end = 1, limits = c(b_min, b_max)) +
  labs(fill = "b", x = "", y = "", title = "effort fitté")

# Vrai effort
ggplot()+
  geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = bvrai))+
  # geom_point(data = dat$loc.det, aes(x = x, y = y), col = "white") +
  scale_fill_viridis_c(begin = 0, end = 1, limits = c(b_min, b_max)) +
  labs(fill = "b", x = "", y = "", title = "vrai effort")

# Différence
ggplot()+
  geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = bvrai-bestim))+
  # geom_point(data = dat$loc.det, aes(x = x, y = y), col = "white") +
  scale_fill_viridis_c(begin = 0, end = 1) +
  labs(fill = "bv-bf", x = "", y = "", title = "différence d'effort")

# Taux d'erreur
summary(abs(
  (bvrai-bestim) / 
    bvrai
))


### Effort normalisé ################

sb_min <- min(scale(bvrai)[,1], scale(bestim)[,1])
sb_max <- max(scale(bvrai)[,1], scale(bestim)[,1])

# Effort fitté
ggplot()+
  geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = scale(bestim)[,1]))+
  # geom_point(data = dat$loc.det, aes(x = x, y = y), col = "white") +
  scale_fill_viridis_c(begin = 0, end = 1, limits = c(sb_min, sb_max)) +
  labs(fill = "b", x = "", y = "", title = "effort fitté, normalisé")

# Vrai effort
ggplot()+
  geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = scale(bvrai)[,1]))+
  # geom_point(data = dat$loc.det, aes(x = x, y = y), col = "white") +
  scale_fill_viridis_c(begin = 0, end = 1, limits = c(sb_min, sb_max)) +
  labs(fill = "b", x = "", y = "", title = "vrai effort, normalisé" )

# Différence
ggplot()+
  geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = scale(bvrai)[,1]-scale(bestim)[,1]))+
  # geom_point(data = dat$loc.ipp, aes(x = x, y = y), col = "white") +
  scale_fill_viridis_c(begin = 0, end = 1) +
  labs(fill = "", x = "", y = "", title = "différence d'efforts normalisés")

# Résumé du taux d'erreur
summary(abs(
  (scale(bvrai)[,1]-scale(bestim)[,1]) 
  / scale(bvrai)[,1]
))


### Intensité amincie #####################

lb_min <- min(lambdavraie*bvrai, lambdaestim*bestim)
lb_max <- max(lambdavraie*bvrai, lambdaestim*bestim)

# intensité amincie fittée
ggplot()+
  geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = bestim*lambdaestim))+
  # geom_point(data = dat$loc.det, aes(x = x, y = y), col = "white") +
  scale_fill_viridis_c(begin = 0, end = 1, limits = c(lb_min, lb_max)) +
  labs(fill = "lb", x = "", y = "", title = "intensité amincie fittée")

# Vraie intensité amincie
ggplot()+
  geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = bvrai*lambdavraie))+
  # geom_point(data = dat$loc.det, aes(x = x, y = y), col = "white") +
  scale_fill_viridis_c(begin = 0, end = 1, limits = c(lb_min, lb_max)) +
  labs(fill = "lb", x = "", y = "", title = "vraie intensité amincie")

# Différence
ggplot()+
  geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = bvrai*lambdavraie-bestim*lambdaestim))+
  # geom_point(data = dat$loc.det, aes(x = x, y = y), col = "white") +
  scale_fill_viridis_c(begin = 0, end = 1) +
  labs(fill = "lbv-lbf", x = "", y = "", title = "différence d'intensités amincies")

# Résumé du taux d'erreur
summary(abs(
  (bvrai*lambdavraie-bestim*lambdaestim) / 
    (bvrai*lambdavraie)
))




# Avec le package spatstat -----------------------------------------------------

# test simu...
pp <- rpoispp(function(x,y) {100 * exp(-3*x)}, 100)
plot(pp)

data_test <- expand.grid(x = seq(0,1,0.01), y = seq(0,1,0.01))
data_test$z <- with(data_test, 100 * exp(-3*x))

ggplot(data_test, aes(x = x, y = y, fill = z)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_viridis_c()

Z <- as.im(function(x,y){sqrt(x+y)/10}, owin(xrange = c(0,100), yrange = c(0,100)), dimyx = 100)
# pp <- rpoispp(Z)
# plot(pp)
plot(Z)
# M <- as.matrix.im(Z)

W <- as.im(function(x,y){0.5*exp(-(y-50)^2/400)}, owin(xrange = c(0,100), yrange = c(0,100)), dimyx = 100)
plot(W)
plot(Z*W)
pp <- rpoispp(Z*W)
plot(pp)
