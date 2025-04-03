# Modèle 1 détection / routes_chemins ------------------------------------------

## Code et data ####

prior de beta et alpha à revoir ! double exponentielle !!!

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
      betagamma[1] * x_1[pixel] +
      betagamma[2] * x_2[pixel] +
      betagamma[3] * x_3[pixel] + 
      betagamma[4] * x_4[pixel] + 
      betagamma[5] * x_5[pixel] + 
      betagamma[6] * x_6[pixel] +
      betagamma[7] * x_7[pixel] +
      betagamma[8] * x_8[pixel] + 
      cell_area[pixel] 

    # presence only thinning prob linear predictor
    #
    # h_s = covariates for thinning probability
    # alpha  = presence-only data model regression coefficients
    #
    logit(b[pixel]) <-  alpha[1] + alphagamma[1] * h_1[pixel] + alphagamma[2] * h_2[pixel]
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
  
  for(i in 1:8){
    betagamma[i] <- beta[i+1] * gamma[i]
    gamma[i] ~ dbern(0.5)
  }
  for(j in 1:2){
    alphagamma[j] <- alpha[j+1] * gam[j]
    gam[j] ~ dbern(0.5)
  }
})

head(grid_selec$grid_id) # les ID des cellules

pixel.id.det <- grid_selec$grid_id[grid_selec$nnutria > 0] # les ID des cellules où il y a au moins une occurrence
head(pixel.id.det)

npix <- nrow(grid_selec)
# avec prise en compte de la surface intersectée
s.area <- as.numeric(units::set_units(grid_selec$area,"km^2"))
logarea <- log(s.area)

# data : choisir entre distance à l'eau ou surface...
data <- list(
  cell_area = logarea,
  # x_1 = scale(grid_selec$surface_en_eau)[,1],
  x_1 = scale(grid_selec$dist_plan_eau)[,1],
  x_2 = scale(grid_selec$logdensity)[,1],
  x_3 = scale(grid_selec$agri_cover)[,1],
  x_4 = scale(grid_selec$temp_min)[,1],
  x_5 = scale(grid_selec$temp_max)[,1],
  x_6 = scale(grid_selec$temp_mean)[,1],
  x_7 = scale(grid_selec$prec_cum)[,1],
  # x_8 = scale(grid_selec$lgr_rivieres)[,1],
  x_8 = scale(grid_selec$dist_rivieres)[,1],
  h_1 = scale(grid_selec$lgr_chemins)[,1],
  h_2 = scale(grid_selec$lgr_routes)[,1],
  ones = rep(1, length(pixel.id.det)))

constants <- list(
  npixel = npix,
  m = length(pixel.id.det), 
  CONSTANT = 50000,
  po_pixel = pixel.id.det) 

# zinit <- numeric(npix)
# zinit[pixel.id.det] <- 1
inits <- function(){
  list(
    beta = rnorm(9, 0, 1), 
    # gamma = rep(1, 8),
    alpha = rnorm(3, 0, 1)
    # gam = rep(1, 2)
    # z = zinit # pourquoi en commentaire dans le code initial ?
  )
}

params <- c("alpha", "beta", "alphagamma", "betagamma")

## MCMC ####

# MCMC settings (pour tester...)
nc <- 2
nburn <- 5000 #5000
ni <- nburn + 10000 #30000
nt <- 1

# set.seed(123) 
# iv <- inits()
# 
# # param b
# b <- plogis(iv$alpha[1] + iv$alpha[2] * data$h_1 + iv$alpha[3] * data$h_2) # pas de valeurs qui explosent
# summary(b)
# 
# # param lambda
# lambda <- exp(iv$beta[1] +
#                 iv$beta[2] * data$x_1 +
#                 iv$beta[3] * data$x_2 +
#                 iv$beta[4] * data$x_3 +
#                 iv$beta[5] * data$x_4 +
#                 iv$beta[6] * data$x_5 +
#                 iv$beta[7] * data$x_6 +
#                 iv$beta[8] * data$x_7 +
#                 iv$beta[9] * data$x_8 + data$cell_area) 
# summary(lambda)
# which.max(lambda)
# data$x_1[which.max(lambda)]
# summary(data$x_7)

set.seed(123)
start <- Sys.time()
out_env <- nimbleMCMC(
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

# save(out_env, file = "out_env_5km2_dist.RData")

MCMCsummary(out_env)

MCMCtrace(out_env, pdf = FALSE, ind = TRUE, params = "alpha")
MCMCplot(out_env, params = "beta")

res <- rbind(out_env$chain1, out_env$chain2)

# # select z
# mask <- str_detect(colnames(res), "z")
# res_z <- res[,mask]
# grid_selec$zestim <- apply(res_z, 2, median)
# grid_selec$zmoy <- apply(res_z, 2, mean)

# # viz
# ggplot() +
#   geom_sf(data = grid_selec, lwd = 0.1, aes(fill = as_factor(zestim))) +
#   labs(fill = "Présence potentielle estimée du ragondin") +
#   geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
#   # geom_sf(data = nutria) +
#   theme_void()
# 
# p <- ggplot() +
#   geom_sf(data = st_intersection(grid_selec, occitanie), lwd = 0.1, aes(fill = zmoy)) +
#   labs(fill = "Présence potentielle estimée du ragondin") +
#   scale_fill_viridis_c(begin = 0, end = 1) +
#   geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
#   # geom_sf(data = nutria) +
#   theme_light()
# p
# # ggsave(plot = p, "Images/map_env.png", dpi = 600)

# select alpha
mask <- str_detect(colnames(res), "alpha")
res_alpha <- res[,mask]
alphaestim <- apply(res_alpha, 2, median)
alphamoy <- apply(res_alpha, 2, mean)

# select beta
mask <- str_detect(colnames(res), "beta")
res_beta <- res[,mask]
betaestim <- apply(res_beta, 2, median)
betamoy <- apply(res_beta, 2, mean)

# lambda et b
grid_selec$lambda <- exp(betaestim[1] +
                           betaestim[2] * data$x_1 +
                           betaestim[3] * data$x_2 +
                           betaestim[4] * data$x_3 +
                           betaestim[5] * data$x_4 +
                           betaestim[6] * data$x_5 +
                           betaestim[7] * data$x_6 +
                           betaestim[8] * data$x_7 +
                           betaestim[9] * data$x_8 + data$cell_area)
grid_selec$b <- plogis(alphaestim[1] +
                         alphaestim[2] * data$h_1 +
                         alphaestim[3] * data$h_2)

# plot
p_lambda <- ggplot() +
  geom_sf(data = st_intersection(grid_selec, occitanie), color = NA, aes(fill = lambda)) +
  labs(fill = "Intensité") +
  scale_fill_viridis_c(begin = 0, end = 1) +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  # geom_sf(data = nutria) +
  theme_light()
p_lambda
# ggsave(plot = p_lambda, "Images/map_int_env_5km2_dist.png", dpi = 600)

p_b <- ggplot() +
  geom_sf(data = st_intersection(grid_selec, occitanie), color = NA, aes(fill = b)) +
  labs(fill = "Effort") +
  scale_fill_viridis_c(begin = 0, end = 1) +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  # geom_sf(data = nutria) +
  theme_light()
p_b
# ggsave(plot = p_b, "Images/map_eff_env_5km2_dist.png", dpi = 600)

p_p <- ggplot() +
  geom_sf(data = st_intersection(grid_selec, occitanie), color= NA, aes(fill = 1-exp(-lambda))) +
  labs(fill = "Présence potentielle estimée du ragondin") +
  scale_fill_viridis_c(begin = 0, end = 1) +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  # geom_sf(data = nutria) +
  theme_light()
p_p
# ggsave(plot = p_p, "Images/map_pres_env_5km2_dist.png", dpi = 600)