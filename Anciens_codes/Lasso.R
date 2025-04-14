# Lasso ------------------------------------------

## Code et data ####

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

    logit(b[pixel]) <-  alpha[1] + 
      alpha[2] * h_1[pixel] 
  
    }
  
  po_denominator <- inprod(lambda[1:npixel], b[1:npixel]) / m

  for(po in 1:m){
    ones[po] ~ dbern(
      exp(
        log(lambda[po_pixel[po]] * b[po_pixel[po]]) -
          po_denominator
        )
      / CONSTANT
      ) 
  }
  
  for(j in 1:2){
    alpha[j] ~ dnorm(0, sd = 2)
  }

  beta[1] ~ dnorm(0, sd = 2) # intercept
  for(i in 2:8){
    beta[i] ~ ddexp(0, tau) 
  }
  tau ~ dunif(0.001,10)

})

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


## MCMC ####

# MCMC settings (pour tester...)
nc <- 2
nburn <- 10000 #5000
ni <- nburn + 10000 #30000
nt <- 1

set.seed(123)
start <- Sys.time()
out <- nimbleMCMC(
  code = code,
  constants = constants,
  data = data,
  inits = inits(),
  monitors = c("alpha", "beta", "tau"),
  niter = ni,
  nburnin = nburn,
  nchains = nc,
  thin = nt
  # WAIC = TRUE
)
end <- Sys.time()
end - start


MCMCsummary(out)

MCMCtrace(out, pdf = FALSE, ind = TRUE, params = "all")
MCMCplot(out, params = "beta")

res <- rbind(out$chain1, out$chain2)



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
                           data$cell_area)
grid_selec$b <- plogis(alphaestim[1] +
                         alphaestim[2] * data$h_1)

# plot
p_lambda <- ggplot() +
  geom_sf(data = st_intersection(grid_selec, occitanie), color = NA, aes(fill = grid_selec$lambda)) +
  labs(fill = "Intensité") +
  scale_fill_viridis_c() +
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
s# ggsave(plot = p_p, "Images/map_pres_env_5km2_dist.png", dpi = 600)



# Lasso avec sélection de variables ? ------------------------------------------

## Code et data ####

code <- nimbleCode({
  
  for(pixel in 1:npixel){
    
    log(lambda[pixel]) <- beta[1] +
      betag[2] * x_1[pixel] +
      betag[3] * x_2[pixel] +
      betag[4] * x_3[pixel] + 
      betag[5] * x_4[pixel] + 
      betag[6] * x_5[pixel] + 
      betag[7] * x_6[pixel] +
      betag[8] * x_7[pixel] +
      cell_area[pixel] 
    
    logit(b[pixel]) <-  alpha[1] + 
      alpha[2] * h_1[pixel] 
    
  }
  
  po_denominator <- inprod(lambda[1:npixel], b[1:npixel]) / m
  
  for(po in 1:m){
    ones[po] ~ dbern(
      exp(
        log(lambda[po_pixel[po]] * b[po_pixel[po]]) -
          po_denominator
      )
      / CONSTANT
    ) 
  }
  
  for(j in 1:2){
    alpha[j] ~ dnorm(0, sd = 2)
  }
  
  beta[1] ~ dnorm(0, sd = 2) # intercept
  for(i in 2:8){
    beta[i] ~ ddexp(0, tau) 
  }
  tau ~ dunif(0.001,10)
  for (j in 1:7) {
    betag[j] <- beta[j+1] * gamma[j]
    gamma[j] ~ dbern(0.5)
  }
  
})

# Initialisation
inits <- function(){
  list(
    beta = rnorm(8, 0, 1), 
    alpha = rnorm(2, 0, 1),
    gamma = rep(1,7) # ???
  )
}

## MCMC ####

# MCMC settings (pour tester...)
nc <- 2
nburn <- 10000 #5000
ni <- nburn + 10000 #30000
nt <- 1

set.seed(123)
start <- Sys.time()
out <- nimbleMCMC(
  code = code,
  constants = constants,
  data = data,
  inits = inits(),
  monitors = c("alpha", "beta", "gamma", "tau"),
  niter = ni,
  nburnin = nburn,
  nchains = nc,
  thin = nt
  # WAIC = TRUE
)
end <- Sys.time()
end - start


MCMCsummary(out)

MCMCtrace(out, pdf = FALSE, ind = TRUE, params = "all")
MCMCplot(out, params = "beta")

res <- rbind(out$chain1, out$chain2)

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
                           data$cell_area)
grid_selec$b <- plogis(alphaestim[1] +
                         alphaestim[2] * data$h_1)

# plot
p_lambda <- ggplot() +
  geom_sf(data = st_intersection(grid_selec, occitanie), color = NA, aes(fill = grid_selec$lambda)) +
  labs(fill = "Intensité") +
  scale_fill_viridis_c() +
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