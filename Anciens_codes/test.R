library(AHMbook)
library(nimble)
library(MCMCvis)
library(raster)
library(ggplot2)

code <- nimbleCode({
  # Bayesian version of the Koshkina (2017) model.
  #
  # The latent-state model
  for(pixel in 1:npixel){
    # latent state linear predictor
    #
    # x_s  = covariates for latent state
    # beta = latent state model regression coefficients
    # cell_area = log area of grid cell 
    #
    
    log(lambda[pixel]) <- beta[1] + beta[2] * x_s[pixel] + cell_area
    # Species presence in a gridcell as a Bernoulli trial
    z[pixel] ~ dbern(1 - exp(-lambda[pixel]))
    # presence only thinning prob linear predictor
    #
    # h_s = covariates for thinning probability
    # alpha  = presence-only data model regression coefficients
    #
    logit(b[pixel]) <-  alpha[1] + alpha[2] * h_s[pixel]
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
  
  # L_i (ones trick) ? 
  for(po in 1:m){
    ones[po] ~ dbern(
      exp(
        log(lambda[po_pixel[po]] * b[po_pixel[po]]) - 
          po_denominator) 
      / CONSTANT) # attention, voir issue https://github.com/mfidino/integrated-occupancy-model/issues/1
  } 
  # Priors for latent state model
  for(i in 1:2){
    beta[i] ~ dnorm(0, sd = 2)
  }
  # Priors for presence-only data model
  for(j in 1:2){
    alpha[j] ~ dnorm(0, sd = 2)
  }
  # # Priors for det/non-det data model
  # for(pa in 1:npar_pa){
  #   a[pa] ~ dlogis(0, 1)
  # }
  # Derived parameter, the number of cells occupied
  zsum <- sum(z[1:npixel])
})


# # EXAMPLES
# # No spatial bias in point pattern data set 1: no effect of covariate W 
# str(dat <- simDataDK(alpha = c(-1, 0)), 1) 
# # Homogeneous point pattern: no effect of covariate X 
# str(dat <- simDataDK(beta = c(6, 0)), 1) # don't make intercept too big 
# # No effect of covariate W on detection in point count sampling 
# str(dat <- simDataDK(gamma = c(0, 0)), 1) 
# # No thinning of the original point pattern 
# str(dat <- simDataDK(alpha = c(20, 0), drop.out.prop.pb = 0), 1)

set.seed(123)

dat <- simDataDK(
  sqrt.npix = 100, 
  alpha = c(-2, -1), 
  beta = c(6, 0.5), 
  drop.out.prop.pb = 0, 
  quadrat.size = 4,
  gamma = c(0, -1.5), 
  nquadrats = 250, 
  nsurveys = 3, 
  show.plot = T)

str(dat, 1)

# comment apparaît la répétition des 3 observations ? 

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

data <- list(
  cell_area = logarea,
  x_s = dat$xcov,
  h_s = dat$wcov,
  ones = rep(1, length(dat$pixel.id.det))) # utiliser dat$N.det ?

constants <- list(
  npixel = dat$npix,
  m = length(dat$pixel.id.det),
  CONSTANT = 10000,
  po_pixel = dat$pixel.id.det)

# initialisation
zinit <- numeric(dat$npix)
zinit[dat$pixel.id.det] <- 1 # vecteur des détections (composé de 0 et 1)
inits <- function(){
  list(
    beta = rnorm(2, 0, 1), 
    alpha = rnorm(2, 0, 1),
    z = zinit
  )
}

params <- c("alpha", "beta", "zsum")

# MCMC settings
nc <- 2
nburn <- 2000  
ni <- nburn + 3000
nt <- 3

# on lance
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

# résultats
MCMCsummary(out)

# essai de plot de l'intensité
res <- rbind(out$chain1, out$chain2)

# select beta
beta_estim <- apply(res[,3:4], 2, median)
beta_moy <- apply(res[,3:4], 2, mean)

# Bin size control + color palette
data_p <- data.frame(x = dat$s.loc[,1], y =  dat$s.loc[,2], l = beta_estim[1]+beta_estim[2]*dat$xcov)
ggplot(data_p, aes(x=x, y=y, fill=l) ) +
  geom_raster(interpolate = TRUE) +
  scale_fill_gradientn(colors = terrain.colors(14)) +
  theme_bw()

#######################
# Ragondins
#######################
# library(AHMbook)
# library(nimble)
# library(MCMCvis)
# library(raster)
library(ggplot2)
library(tidyverse)
library(sf)

dat_poly <- st_read("../Data/Data_poly.shp") %>%
  mutate(year = year(as.Date(DateDebut)),
         sp = if_else(NomVernacu == "Ragondin", "blue", "green")) %>%
  filter(year > 2003)

str(dat_poly)

dat_lin <- st_read("../Data/Data_lin.shp")

dat_points <- st_read("../Data/Data_pts_test_infos.shp") %>%
  mutate(year = year(as.Date(jourdebut)),
         sp = if_else(nomvern == "Ragondin", "blue", "green"),
         type = if_else(is.na(obsctx), 0, 1)) %>% # 0 = pas d'info, 1 = camtrap (surtout)
  filter(year > 2003)

dat_vu <- dat_points[which(str_detect(dat_points$obsctx,"vu")),]

occitanie <- st_read("../Data/departements-d-occitanie.shp") %>%
  st_transform(crs = st_crs(dat_poly))

# plot
p <- ggplot() +
  geom_sf(data = occitanie, fill = "white", color = "black", lwd = .5) + 
  geom_sf(data = dat_poly) + 
  geom_sf(data = dat_points, aes(color = nomvern)) + 
  labs(color = "") +
  #  geom_sf(data = river_lines, color = "blue", lwd = 0.6) + 
  #  geom_sf(data = dat_lin) + 
  facet_wrap(~year, nrow = 5) + 
  theme_void()
p

#ggsave(plot = p, "fig/map.png", dpi = 600)

unique(dat_points$obsctx)
# contient des individus morts et des signes de présence uniquement... 

ggplot() +
  geom_sf(data = occitanie, fill = "white", color = "black", lwd = .5) + 
  geom_sf(data = dat_points %>% 
            mutate(year = year(as.Date(jourdebut))) %>% 
            filter(year == 2016) %>%
            filter(nomvern == "Ragondin")) + 
  #  labs(color = "") +
  theme_void()

# ragondins en 2016
nutria <- dat_points %>% 
  mutate(year = year(as.Date(jourdebut))) %>% 
  filter(year == 2016) %>%
  filter(nomvern == "Ragondin")

# départements de la région
dpts_occitanie <- st_read("../Data/departements-d-occitanie.shp") %>%
  st_transform(crs = st_crs(dat_points))

# contours de la région
occitanie <- dpts_occitanie %>% st_union()

# dataviz
ggplot() +
  geom_sf(data = occitanie, fill = "white", color = "black", lwd = .5) + 
  geom_sf(data = nutria) + 
  theme_void()



# https://urbandatapalette.com/post/2021-08-tessellation-sf/

grid <- st_make_grid(occitanie, n = 25, what = "polygons", square = FALSE)

# To sf and add grid ID
grid_sf <- st_sf(grid) %>%
  # add grid ID
  mutate(grid_id = 1:length(lengths(grid)))

ggplot() +
  geom_sf(data = grid_sf, lwd= 0.1) + 
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) + 
  theme_void()

grid_sf <- grid_sf[st_intersects(grid_sf, 
                                 occitanie, 
                                 sparse = FALSE), ]
# on renumérote
grid_sf <- grid_sf %>%
  mutate(grid_id = 1:lengths(grid_sf)[1])

ggplot() +
  geom_sf(data = grid_sf, lwd = 0.1) +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) + 
  theme_void()

# https://gis.stackexchange.com/questions/323698/counting-points-in-polygons-with-sf-package-of-r
grid_sf$nnutria <- lengths(st_intersects(grid_sf, nutria))

table(grid_sf$nnutria)

# remove grid without value of 0 (i.e. no points in side that grid)
nutria_count <- filter(grid_sf, nnutria > 0)

nutria_count %>% 
  ggplot() + 
  geom_sf(data = grid_sf, lwd = 0.1) +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) + 
  geom_sf(aes(fill=nnutria), color="#FFFFFF99") +
  scale_fill_viridis_c(direction = -1)+ 
  labs(fill = "number of removed \ncoypus") +
  geom_sf(data = nutria) + 
  theme_void()

# variables explicatives
roads <- st_read("../Data/TRONCON_ROUTE.shp")
str(roads)
unique(roads$CLASS_ADM)

routes_occ <- roads %>% 
  st_transform(crs = st_crs(grid_sf)) %>%
  st_intersection(grid_sf) %>%
  filter(CLASS_ADM %in% c("D\xe9partementale","Nationale"))

ggplot() +
  geom_sf(data = grid_sf, fill = "white", color = "black", lwd = .5) + 
  geom_sf(data = routes_occ, color = "lightblue") + 
  theme_void()

# calcule la longueur des segments
routes_occ$length <- st_length(routes_occ)

# somme des segments par cellule
longueur_par_cellule <- routes_occ %>%
  group_by(grid_id) %>%  # Remplacez id_cellule par votre nom de colonne
  summarise(lgr_totale = sum(length))  # Somme des longueurs

# jointure
grid_sf$lgr_routes <- 0
for (i in 1:nrow(grid_sf)){
  mask <- grid_sf$grid_id[i]
  if (sum(longueur_par_cellule$grid_id == mask) == 0) next
  grid_sf$lgr_routes[i] <- longueur_par_cellule$lgr_totale[longueur_par_cellule$grid_id == mask]
}

ggplot() +
  geom_sf(data = grid_sf, lwd = 0.1, aes(fill = lgr_routes/1000)) +
  labs(fill = "Longueur de routes (km)") + 
  scale_fill_viridis_c() + 
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) + 
  theme_void()

Ariege <- sf::st_read("../Data/Rivieres/Ariege/COURS_D_EAU.shp")
Aude <- sf::st_read("../Data/Rivieres/Aude/COURS_D_EAU.shp")
Aveyron <- sf::st_read("../Data/Rivieres/Aveyron/COURS_D_EAU.shp")
Gard <- sf::st_read("../Data/Rivieres/Gard/COURS_D_EAU.shp")
HauteGaronne <- sf::st_read("../Data/Rivieres/HauteGaronne/COURS_D_EAU.shp")
Gers <- sf::st_read("../Data/Rivieres/Gers/COURS_D_EAU.shp")
Herault <- sf::st_read("../Data/Rivieres/Herault/COURS_D_EAU.shp")
Lot <- sf::st_read("../Data/Rivieres/Lot/COURS_D_EAU.shp")
Lozere <- sf::st_read("../Data/Rivieres/Lozere/COURS_D_EAU.shp")
HautesPyrenees <- sf::st_read("../Data/Rivieres/HautesPyrenees/COURS_D_EAU.shp")
PyreneesOrientales <- sf::st_read("../Data/Rivieres/PO/COURS_D_EAU.shp")
Tarn <- sf::st_read("../Data/Rivieres/Tarn/COURS_D_EAU.shp")
TarnetGaronne <- sf::st_read("../Data/Rivieres/TarnEtGaronne/COURS_D_EAU.shp")

# bind river
river_lines <- Ariege %>%
  rbind(Aude) %>%
  rbind(Aveyron) %>%
  rbind(Gard) %>%
  rbind(Gers) %>%
  rbind(HauteGaronne) %>%
  rbind(HautesPyrenees) %>%
  rbind(Herault) %>%
  rbind(Lot) %>%
  rbind(Lozere) %>%
  rbind(PyreneesOrientales) %>%
  rbind(Tarn) %>%
  rbind(TarnetGaronne) %>%
  sf::st_transform(crs = st_crs(occitanie)) %>%
  sf::st_simplify()

# checks
nrow(river_lines %>%
       filter(IMPORTANCE %in% c(3,4)))

# rivières d'Occitanie
rivers_occ <- river_lines %>% 
  st_transform(crs = st_crs(grid_sf)) %>%
  st_intersection(grid_sf)

# longueur de chaque segment
rivers_occ$length <- st_length(rivers_occ)  

# somme des segments par cellule
longueur_par_cellule <- rivers_occ %>%
  group_by(grid_id) %>%  # Remplacez id_cellule par votre nom de colonne
  summarise(lgr_totale = sum(length))  # Somme des longueurs

# jointure
grid_sf$lgr_rivieres <- 0
for (i in 1:nrow(grid_sf)){
  mask <- grid_sf$grid_id[i]
  if (sum(longueur_par_cellule$grid_id == mask) == 0) next
  grid_sf$lgr_rivieres[i] <- longueur_par_cellule$lgr_totale[longueur_par_cellule$grid_id == mask]
}

# viz
ggplot() +
  geom_sf(data = grid_sf, lwd = 0.1, aes(fill = lgr_rivieres/1000)) +
  labs(fill = "Longueur de cours d'eau (km)") + 
  scale_fill_viridis_c() + 
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) + 
  theme_void()

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
    
    log(lambda[pixel]) <- beta[1] + beta[2] * x_s[pixel] + cell_area
    # Species presence in a gridcell as a Bernoulli trial
    z[pixel] ~ dbern(1 - exp(-lambda[pixel]))
    # presence only thinning prob linear predictor
    #
    # h_s = covariates for thinning probability
    # alpha  = presence-only data model regression coefficients
    #
    logit(b[pixel]) <-  alpha[1] + alpha[2] * h_s[pixel]
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
  for(i in 1:2){
    beta[i] ~ dnorm(0, sd = 2)
  }
  # Priors for presence-only data model
  for(j in 1:2){
    alpha[j] ~ dnorm(0, sd = 2)
  }
  # Derived parameter, the number of cells occupied
  zsum <- sum(z[1:npixel])
})

head(grid_sf$grid_id) # les ID de toutes les cellules

pixel.id.det <- grid_sf$grid_id[grid_sf$nnutria > 0] # les ID des cellules où il y a au moins une occurrence
head(pixel.id.det)

npix <- nrow(grid_sf)
s.area <- as.numeric(units::set_units(st_area(grid_sf)[1],"km^2"))
(logarea <- log(s.area / npix))

data <- list(
  cell_area = logarea,
  x_s = (grid_sf$lgr_rivieres - mean(grid_sf$lgr_rivieres))/sd(grid_sf$lgr_rivieres), # distribution
  h_s = (grid_sf$lgr_routes - mean(grid_sf$lgr_routes))/sd(grid_sf$lgr_routes), # detection
  ones = rep(1, length(pixel.id.det)))

constants <- list(
  npixel = npix,
  m = length(pixel.id.det), # à modifier pour prendre en compte toutes les détections
  CONSTANT = 10000,
  po_pixel = pixel.id.det) # à modifier pour prendre en compte toutes les détections

zinit <- numeric(npix)
zinit[pixel.id.det] <- 1
inits <- function(){
  list(
    beta = rnorm(2, 0, 1), 
    alpha = rnorm(2, 0, 1)
    #z = zinit
  )
}

params <- c("alpha", "beta", "z")

# MCMC settings (pour tester...)
nc <- 2
nburn <- 1000 #5000  
ni <- nburn + 1000 #30000
nt <- 1

start <- Sys.time()
# run the model!
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

res <- out
#res <- rbind(out$chain1, out$chain2)


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
  geom_sf(data = nutria) + 
  theme_void()

ggplot() +
  geom_sf(data = grid_sf, lwd = 0.1, aes(fill = zmoy)) +
  labs(fill = "Présence potentielle estimée du ragondin") + 
  scale_fill_viridis_c() + 
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) + 
  geom_sf(data = nutria) + 
  theme_void()


###################################
# Test avec multiples détections par cellules...
###################################

# problèmes : messages de présence de NA, estimation des paramètres bizarre (NA),
# très variable d'une exécution à l'autre... 

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
    
    log(lambda[pixel]) <- beta[1] + beta[2] * x_s[pixel] + cell_area
    # Species presence in a gridcell as a Bernoulli trial
    # z[pixel] ~ dbern(1 - exp(-lambda[pixel]))
    # presence only thinning prob linear predictor
    #
    # h_s = covariates for thinning probability
    # alpha  = presence-only data model regression coefficients
    #
    logit(b[pixel]) <-  alpha[1] + alpha[2] * h_s[pixel]
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
  for(i in 1:2){
    beta[i] ~ dnorm(0, sd = 2)
  }
  # Priors for presence-only data model
  for(j in 1:2){
    alpha[j] ~ dnorm(0, sd = 2)
  }
  # Derived parameter, the number of cells occupied
  # zsum <- sum(z[1:npixel])
})

head(grid_sf$grid_id) # les ID de toutes les cellules

pixel.id.det <- grid_sf$grid_id[grid_sf$nnutria > 0] # les ID des cellules où il y a au moins une occurrence
head(pixel.id.det)

npix <- nrow(grid_sf)
s.area <- as.numeric(units::set_units(st_area(grid_sf)[1],"km^2"))
(logarea <- log(s.area / npix))

data <- list(
  cell_area = logarea,
  x_s = (grid_sf$lgr_rivieres - mean(grid_sf$lgr_rivieres))/sd(grid_sf$lgr_rivieres), # distribution
  h_s = (grid_sf$lgr_routes - mean(grid_sf$lgr_routes))/sd(grid_sf$lgr_routes), # detection
  ones = rep(1, length(pixel.id.det)))

constants <- list(
  npixel = npix,
  m = length(pixel.id.det), 
  CONSTANT = 10000,
  po_pixel = pixel.id.det,
  nobs_pixel = grid_sf$nnutria[grid_sf$nnutria > 0]) # modifié pour prendre en compte toutes les détections

# zinit <- numeric(npix)
# zinit[pixel.id.det] <- 1
inits <- function(){
  list(
    beta = rnorm(2, 0, 1), 
    alpha = rnorm(2, 0, 1)
    # z = zinit # pourquoi mis en commentaire dans code initial ?
  )
}

params <- c("alpha", "beta") #, "z")

# MCMC settings (pour tester...)
nc <- 2
nburn <- 500 # 5000
ni <- nburn + 1000 #10000
nt <- 1

start <- Sys.time()
# run the model!
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
  geom_sf(data = nutria) + 
  theme_void()

ggplot() +
  geom_sf(data = grid_sf, lwd = 0.1, aes(fill = zmoy)) +
  labs(fill = "Présence potentielle estimée du ragondin") + 
  scale_fill_viridis_c() + 
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) + 
  geom_sf(data = nutria) + 
  theme_void()
