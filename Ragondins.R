#######################
# Ragondins
#######################

# libraries utiles

# library(nimble)
# library(MCMCvis)
# library(raster)
library(ggplot2)
library(tidyverse)
library(sf)
library(tidyterra)
library(lubridate)
library(KrigR)

#######################
# Occitanie

# contours des départements d'Occitanie
dpts_occitanie <- st_read("Data/departements-d-occitanie.shp") 

# contours de la région
occitanie <- dpts_occitanie %>% st_union()

#######################
# Import des données de ragondins

# import des données (points uniquement)
dat_points <- st_read("Data/Data_pts_test_infos.shp") %>%
  mutate(year = year(as.Date(jourdebut)),
         sp = if_else(nomvern == "Ragondin", "blue", "green"),
         type = if_else(is.na(obsctx), 0, 1)) %>% # 0 = pas d'info, 1 = camtrap (surtout)
  filter(year > 2003)

# plot des données annuelles
p <- ggplot() +
  geom_sf(data = occitanie, fill = "white", color = "black", lwd = .5) + 
  geom_sf(data = dat_points, aes(color = nomvern)) + 
  labs(color = "") +
  #  geom_sf(data = river_lines, color = "blue", lwd = 0.6) + 
  #  geom_sf(data = dat_lin) + 
  facet_wrap(~year, nrow = 5) + 
  theme_void()
p
# ggsave(plot = p, "fig/map.png", dpi = 600)

table(dat_points$year, dat_points$nomvern)

# changement de format
dpts_occitanie <- st_read("Data/departements-d-occitanie.shp") %>%
  st_transform(crs = st_crs(dat_points))
occitanie <- dpts_occitanie %>% st_union()

#####################
# Ragondins - année 2019 (avec le plus d'observations)

# Sélection des données sur les ragondins seuls en 2019
nutria <- dat_points %>% 
  mutate(year = year(as.Date(jourdebut))) %>% 
  filter(year == 2019) %>%
  filter(nomvern == "Ragondin")

# visualisation
ggplot() +
  geom_sf(data = dpts_occitanie, fill = "white", color = "black", lwd = .5) + 
  geom_sf(data = nutria) + 
  theme_void()

#####################
# Création de la grille
# https://urbandatapalette.com/post/2021-08-tessellation-sf/
grid <- st_make_grid(occitanie, n = 25, what = "polygons", square = FALSE)

# sf
grid_sf <- st_sf(grid)

ggplot() +
  geom_sf(data = grid_sf, lwd= 0.1) + 
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) + 
  theme_void()

# Cellules utiles uniquement
grid_sf <- grid_sf[st_intersects(grid_sf, 
                                 occitanie, 
                                 sparse = FALSE), ]

ggplot() +
  geom_sf(data = grid_sf, lwd = 0.1) +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) + 
  theme_void()

# Ajout de l'identifiant des cellules
grid_sf <- grid_sf %>%
  mutate(grid_id = 1:lengths(grid_sf)[1])

# Nombre d'observations par cellule
# https://gis.stackexchange.com/questions/323698/counting-points-in-polygons-with-sf-package-of-r
grid_sf$nnutria <- lengths(st_intersects(grid_sf, nutria))

table(grid_sf$nnutria)

# cellules ayant au moins une observation
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

######################
# Variables explicatives
######################

Dir.Base <- getwd() # identifying the current directory
Dir.Data <- file.path(Dir.Base, "Data") # folder path for data

# contours de la région
loc_site <- occitanie %>% 
  st_transform(crs = 4326)

st_bbox(loc_site)

occitanie_raster <- st_bbox(c(xmin = -0.4, 
                              xmax = 5, 
                              ymax = 45.1, 
                              ymin = 42), 
                            crs = st_crs(4326)) %>%
  st_as_sfc() %>%
  terra::vect() %>%
  terra::rast()

######################
# Températures
toccitanie <- CDownloadS(
  Variable = "2m_temperature",
  DataSet = "reanalysis-era5-land-monthly-means",
  Type = "monthly_averaged_reanalysis",
  DateStart = "2003-01-01 00:00",
  DateStop = "2024-12-31 23:00",
  TZone = "Europe/Paris",
  TResolution = "month",
  TStep = 1,
  Extent = occitanie_raster, # our data.frame with Lat and Lon columns
  Dir = Dir.Data,
  FileName = "Toccitanie",
  API_User = "camille.mottier@cefe.cnrs.fr",
  API_Key = "4700337b-89cd-42e8-98b6-a450f702d332"
)

toccitanie <- weathermetrics::kelvin.to.celsius(toccitanie)

# save(toccitanie, file = "Data/toccitanie.RData")
# load("Data/toccitanie.RData") marche pas...

# variables annuelles de 2003 (1) à 2024 (22)
# mois min
toccitanie_min <- terra::aggregate(toccitanie, fact=c(1,1,12), fun=min)
# mois max 
toccitanie_max <- terra::aggregate(toccitanie, fact=c(1,1,12), fun=max)
# mois mean
toccitanie_mean <- terra::aggregate(toccitanie, fact=c(1,1,12), fun=mean)

# moyenne des voisins pour cellules manquantes (à voir si nécessaire)
weight_matrix <- matrix(1, nrow=3, ncol=3)
weight_matrix[2,2] <- NA
toccitanie_min <- terra::focal(toccitanie_min, w=weight_matrix, fun = mean, na.policy = "only", na.rm = T)
toccitanie_max <- terra::focal(toccitanie_max, w=weight_matrix, fun = mean, na.policy = "only", na.rm = T)
toccitanie_mean <- terra::focal(toccitanie_mean, w=weight_matrix, fun = mean, na.policy = "only", na.rm = T)

# plot de 2019 (17)
ggplot() +
  geom_spatraster(data=toccitanie_min$lyr.17) +
  scale_fill_viridis_c() +
  geom_sf(data = loc_site, 
          fill = NA)

ggplot() +
  geom_spatraster(data=toccitanie_max$lyr.17) +
  scale_fill_viridis_c() +
  geom_sf(data = loc_site, 
          fill = NA)

ggplot() +
  geom_spatraster(data=toccitanie_mean$lyr.17) +
  scale_fill_viridis_c() +
  geom_sf(data = loc_site, 
          fill = NA)

# avec la grille (2019)
temp_min <- terra::extract(toccitanie_min, grid_sf, 
                           # exact = T, # pour toutes les cellules... 
                           ) %>%
  group_by(ID) %>%
  summarise(tmin2019 = mean(lyr.17, na.rm = T))

# pour toutes les années d'un coup
temp_min_all <- terra::extract(toccitanie_min, grid_sf, fun = mean, na.rm = T,
                               # exact = T # pour utiliser toutes les cases qui intersectent les cellules
                               # voir https://tmieno2.github.io/R-as-GIS-for-Economists/extracting-values-from-raster-layers-for-vector-data.html
                               ) 

ggplot() +
  geom_sf(data = grid_sf, lwd = 0.1, aes(fill = temp_min$tmin2019)) +
  scale_fill_viridis_c() +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()
# problème : des cellules manquantes... 

temp_max <- terra::extract(toccitanie_max, grid_sf, 
                           # exact = T, # pour toutes les cellules... 
) %>%
  group_by(ID) %>%
  summarise(tmax2019 = mean(lyr.17, na.rm = T))

temp_mean <- terra::extract(toccitanie_mean, grid_sf, 
                           # exact = T, # pour toutes les cellules... 
) %>%
  group_by(ID) %>%
  summarise(tmean2019 = mean(lyr.17, na.rm = T))

# ajout des covariables dans grif_sf
grid_sf$temp_min <- temp_min$tmin2019
grid_sf$temp_max <- temp_max$tmax2019
grid_sf$temp_mean <- temp_mean$tmean2019

######################
# Précipitations
poccitanie <- CDownloadS(
  Variable = "total_precipitation",
  DataSet = "reanalysis-era5-land-monthly-means",
  CumulVar = TRUE,
  Type = "monthly_averaged_reanalysis",
  DateStart = "2003-01-01 00:00",
  DateStop = "2024-12-31 23:00",
  TZone = "Europe/Paris",
  TResolution = "month",
  TStep = 1,
  Extent = occitanie_raster, # our data.frame with Lat and Lon columns
  Dir = Dir.Data,
  FileName = "Poccitanie",
  API_User = "camille.mottier@cefe.cnrs.fr",
  API_Key = "4700337b-89cd-42e8-98b6-a450f702d332"
)

# save(poccitanie, file = "Data/poccitanie.RData")
# load("Data/poccitanie.RData")

# Cumuls annuels
poccitanie_cum <- terra::aggregate(poccitanie, fact=c(1,1,12), fun=sum)
# Moyenne pour les cellules manquantes
poccitanie_cum <- terra::focal(poccitanie_cum, w=weight_matrix, fun = mean, na.policy = "only", na.rm = T)

# plot de 2019 (17)
ggplot() +
  geom_spatraster(data=poccitanie_cum$lyr.17) +
  scale_fill_viridis_c() +
  geom_sf(data = loc_site, 
          fill = NA)

# avec la grille (2019)
prec_cum <- terra::extract(poccitanie_cum, grid_sf) %>%
  group_by(ID) %>%
  summarise(pcum2019 = mean(lyr.17, na.rm = T)
  )

# pour toutes les années d'un coup
prec_cum_all <- terra::extract(poccitanie_cum, grid_sf, fun = mean) 

ggplot() +
  geom_sf(data = grid_sf, lwd = 0.1, aes(fill = prec_cum$pcum2019)) +
  scale_fill_viridis_c() +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()
# problème : des cellules manquantes... 

# Ajout de la covariable à grid_sf
grid_sf$prec_cum <- prec_cum$pcum2019

######################
# Variables de sol
grid_sf <- grid_sf %>%
  mutate(area = st_area(.), .before = 1)

# Données de CORINE révisées 2012 (faut-il tenir compte de l'évolution ?)
Cori12 <- st_read("Data/CLC12_RLRMP_RGF.shp")
Cori12 <- Cori12 %>%
  st_transform(crs = st_crs(grid_sf)) %>%
  st_crop(occitanie) 

str(Cori12)
unique(Cori12$CODE_12)

# Proportion de surface agricole 
agri <- Cori12 %>%
  filter(str_sub(CODE_12,1,1)=="2")

grid_agri <- agri %>%
  st_intersection(grid_sf) %>%
  mutate(area = st_area(.)) %>%
  group_by(grid_id) %>%
  summarise(aera_agri = sum(area)) %>%
  as_tibble() %>%
  select(-geometry)

grid_sf <- grid_sf %>% 
  full_join(grid_agri, by = "grid_id") %>%
  mutate(.before = 1, agri_cover = aera_agri/area) %>%
  select(-aera_agri)

ggplot() +
  geom_sf(data = grid_sf, lwd = 0.1, aes(fill = as.numeric(agri_cover))) +
  scale_fill_viridis_c(
    labels = scales::percent_format()
  ) +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()

# Cellules manquantes (différent de 0) ! 

# # proportion des zones humides -> trop incomplet
# # faire le point avec Olivier sur les codes à utiliser
# humide <- Cori12 %>%
#   filter(CODE_12 %in% c("411", "412", "511", "512"))
# 
# grid_hum <- humide %>%
#   st_intersection(grid_sf) %>%
#   mutate(area = st_area(.)) %>%
#   group_by(grid_id) %>%
#   summarise(aera_hum = sum(area)) %>%
#   as_tibble() %>%
#   select(-geometry)
# 
# grid_sf <- grid_sf %>% 
#   full_join(grid_hum, by = "grid_id") %>%
#   mutate(.before = 1, hum_cover = aera_hum/area)
# 
# rm(list = c("agri", "Cori12", "humide"))
# gc()
# 
# ggplot() +
#   geom_sf(data = grid_sf, lwd = 0.1, aes(fill = as.numeric(hum_cover))) +
#   scale_fill_viridis_c(
#     labels = scales::percent_format()
#   ) +
#   geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
#   theme_void()

# cellules manquantes !


######################
# Densité de population 

pop <- st_read("Data/pop2021.gpkg")
pop <- pop %>% 
  st_transform(crs = st_crs(grid_sf)) %>%
  st_crop(occitanie)

grid_pop <- pop %>%
  st_intersection(grid_sf) %>%
  group_by(grid_id) %>%
  summarise(hab = sum(TOT_P_2021)) %>%
  as_tibble() %>%
  select(-geom)

grid_sf <- grid_sf %>% 
  full_join(grid_pop, by = "grid_id") %>%
  mutate(.before = 1, logdensity = log(hab/area)) %>%
  select(-hab)

# problème des cellules -inf (vides) et NA (cellules manquantes sur la frontière)... 

ggplot() +
  geom_sf(data = grid_sf, lwd = 0.1, aes(fill = as.numeric(logdensity))) +
  scale_fill_viridis_c() +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()

######################
# Longueur de chemins

# OpenStreetMap https://www.data.gouv.fr/fr/datasets/itineraires-de-randonnee-dans-openstreetmap/
chemins_osm <- st_read("Data/hiking_foot_routes_lineLine.shp")
chemins_osm <- chemins_osm %>% 
  st_transform(crs = st_crs(grid_sf)) 

# une même ligne peut apparaître plusieurs fois ? 

chemins_occ <- chemins_osm %>%
  st_intersection(grid_sf)

ggplot() +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  geom_sf(data = chemins_occ, color = "grey") +
  theme_void()

# calcule la longueur des segments
chemins_occ$length <- st_length(chemins_occ)

# somme des segments par cellule
longueur_par_cellule <- chemins_occ %>%
  group_by(grid_id) %>% 
  summarise(lgr_totale = sum(length))  

# jointure
grid_sf$lgr_chemins <- 0
for (i in 1:nrow(grid_sf)){
  mask <- grid_sf$grid_id[i]
  if (sum(longueur_par_cellule$grid_id == mask) == 0) next
  grid_sf$lgr_chemins[i] <- longueur_par_cellule$lgr_totale[longueur_par_cellule$grid_id == mask]
}

ggplot() +
  geom_sf(data = grid_sf, lwd = 0.1, aes(fill = lgr_chemins/1000)) +
  labs(fill = "Longueur de chemins (km)") +
  scale_fill_viridis_c() +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()

######################
# Longueur de routes

roads <- st_read("Data/TRONCON_ROUTE.shp")
str(roads)
unique(roads$CLASS_ADM)

routes_occ <- roads %>%
  st_transform(crs = st_crs(grid_sf)) %>%
  st_intersection(grid_sf) %>%
  filter(CLASS_ADM %in% c("D\xe9partementale","Nationale"))

ggplot() +
  geom_sf(data = grid_sf, fill = "white", color = "black", lwd = .5) +
  geom_sf(data = routes_occ, color = "gray") +
  theme_void()

# calcule la longueur des segments
routes_occ$length <- st_length(routes_occ)

# somme des segments par cellule
longueur_par_cellule <- routes_occ %>%
  group_by(grid_id) %>%
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

######################
# Longueur de cours d'eau
# https://bdtopoexplorer.ign.fr/detail_hydrographique

Ariege <- sf::st_read("Data/Rivieres/Ariege/COURS_D_EAU.shp")
Aude <- sf::st_read("Data/Rivieres/Aude/COURS_D_EAU.shp")
Aveyron <- sf::st_read("Data/Rivieres/Aveyron/COURS_D_EAU.shp")
Gard <- sf::st_read("Data/Rivieres/Gard/COURS_D_EAU.shp")
HauteGaronne <- sf::st_read("Data/Rivieres/HauteGaronne/COURS_D_EAU.shp")
Gers <- sf::st_read("Data/Rivieres/Gers/COURS_D_EAU.shp")
Herault <- sf::st_read("Data/Rivieres/Herault/COURS_D_EAU.shp")
Lot <- sf::st_read("Data/Rivieres/Lot/COURS_D_EAU.shp")
Lozere <- sf::st_read("Data/Rivieres/Lozere/COURS_D_EAU.shp")
HautesPyrenees <- sf::st_read("Data/Rivieres/HautesPyrenees/COURS_D_EAU.shp")
PyreneesOrientales <- sf::st_read("Data/Rivieres/PO/COURS_D_EAU.shp")
Tarn <- sf::st_read("Data/Rivieres/Tarn/COURS_D_EAU.shp")
TarnetGaronne <- sf::st_read("Data/Rivieres/TarnEtGaronne/COURS_D_EAU.shp")

# on regroupe
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

rm(list=c("Ariege", "Aude","Aveyron","Gard","Gers","HauteGaronne",
          "HautesPyrenees","Herault","Lot","Lozere","PyreneesOrientales","Tarn","TarnetGaronne"))

# rivières d'Occitanie
rivers_occ <- river_lines %>%
  st_transform(crs = st_crs(grid_sf)) %>%
  st_intersection(grid_sf)

# longueur de chaque segment
rivers_occ$length <- st_length(rivers_occ)

# somme des segments par cellule
longueur_par_cellule <- rivers_occ %>%
  group_by(grid_id) %>%  
  summarise(lgr_totale = sum(length))  # Somme des longueurs

# jointure
grid_sf$lgr_rivieres <- 0
for (i in 1:nrow(grid_sf)){
  mask <- grid_sf$grid_id[i]
  if (sum(longueur_par_cellule$grid_id == mask) == 0) next
  grid_sf$lgr_rivieres[i] <- longueur_par_cellule$lgr_totale[longueur_par_cellule$grid_id == mask]
}

# visualisation
ggplot() +
  geom_sf(data = grid_sf, lwd = 0.1, aes(fill = lgr_rivieres/1000)) +
  labs(fill = "Longueur de cours d'eau (km)") +
  scale_fill_viridis_c() +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()

######################
# Proportion de plan d'eau
# https://bdtopoexplorer.ign.fr/detail_hydrographique

# Plan d'eau
Ariege <- sf::st_read("Data/Rivieres/Ariege/PLAN_D_EAU.shp")
Aude <- sf::st_read("Data/Rivieres/Aude/PLAN_D_EAU.shp")
Aveyron <- sf::st_read("Data/Rivieres/Aveyron/PLAN_D_EAU.shp")
Gard <- sf::st_read("Data/Rivieres/Gard/PLAN_D_EAU.shp")
HauteGaronne <- sf::st_read("Data/Rivieres/HauteGaronne/PLAN_D_EAU.shp")
Gers <- sf::st_read("Data/Rivieres/Gers/PLAN_D_EAU.shp")
Herault <- sf::st_read("Data/Rivieres/Herault/PLAN_D_EAU.shp")
Lot <- sf::st_read("Data/Rivieres/Lot/PLAN_D_EAU.shp")
Lozere <- sf::st_read("Data/Rivieres/Lozere/PLAN_D_EAU.shp")
HautesPyrenees <- sf::st_read("Data/Rivieres/HautesPyrenees/PLAN_D_EAU.shp")
PyreneesOrientales <- sf::st_read("Data/Rivieres/PO/PLAN_D_EAU.shp")
Tarn <- sf::st_read("Data/Rivieres/Tarn/PLAN_D_EAU.shp")
TarnetGaronne <- sf::st_read("Data/Rivieres/TarnEtGaronne/PLAN_D_EAU.shp")

# On regroupe
plan_eau <- Ariege %>%
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

rm(list=c("Ariege", "Aude","Aveyron","Gard","Gers","HauteGaronne",
          "HautesPyrenees","Herault","Lot","Lozere","PyreneesOrientales","Tarn","TarnetGaronne"))

# à voir s'il faut faire un tri
unique(plan_eau$NATURE)

# avec la grille
plan_eau_occ <- plan_eau %>%
  st_transform(crs = st_crs(grid_sf)) %>%
  st_intersection(grid_sf) %>%
  mutate(area = st_area(.)) %>%
  group_by(grid_id) %>%
  summarise(aera_eau = sum(area)) %>%
  as_tibble() %>%
  select(-geometry)

grid_sf <- grid_sf %>% 
  full_join(plan_eau_occ, by = "grid_id") %>%
  mutate(.before = 1, surface_en_eau = aera_eau/area) %>%
  select(-aera_eau)

ggplot() +
  geom_sf(data = grid_sf, lwd = 0.1, aes(fill = as.numeric(surface_en_eau))) +
  scale_fill_viridis_c(
    labels = scales::percent_format()
  ) +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()

# save(grid_sf, file = "grid_sf.RData")

# # Bayesian version of the Koshkina (2017) model.
# code <- nimbleCode({
#   # latent-state model
#   for(pixel in 1:npixel){
#     # latent state linear predictor
#     #
#     # x_s  = covariates for latent state
#     # beta = latent state model regression coefficients
#     # cell_area = log area of grid cell 
#     #
#     
#     log(lambda[pixel]) <- beta[1] + beta[2] * x_s[pixel] + cell_area
#     # Species presence in a gridcell as a Bernoulli trial
#     z[pixel] ~ dbern(1 - exp(-lambda[pixel]))
#     # presence only thinning prob linear predictor
#     #
#     # h_s = covariates for thinning probability
#     # alpha  = presence-only data model regression coefficients
#     #
#     logit(b[pixel]) <-  alpha[1] + alpha[2] * h_s[pixel]
#   }
#   # The presence only data model.
#   #
#   # This part of the model just uses the
#   #  what we have calculated above (lambda
#   #  and b). The denominator of this likelihood
#   #  is actually a scalar so we can calculate it
#   #  outside of a for loop. Let's do that first.
#   #
#   # The presence_only data model denominator, which
#   #  is the thinned poisson process across the
#   #  whole region (divided by the total number of 
#   #  data points because it has to be 
#   #  evaluated for each data point).
#   # m is the number of presence-only data points
#   po_denominator <- inprod(lambda[1:npixel], b[1:npixel]) / m
#   #
#   # Loop through each presence-only data point
#   #  using Bernoulli one's trick. The numerator
#   #  is just the thinned poisson process for
#   #  the ith data point.
#   #  po_pixel denotes the grid cell of the ith presence only data point
#   for(po in 1:m){
#     ones[po] ~ dbern(
#       exp(
#         log(lambda[po_pixel[po]] * b[po_pixel[po]]) - 
#           po_denominator) 
#       / CONSTANT) # attention, voir issue https://github.com/mfidino/integrated-occupancy-model/issues/1
#   } 
#   # Priors for latent state model
#   for(i in 1:2){
#     beta[i] ~ dnorm(0, sd = 2)
#   }
#   # Priors for presence-only data model
#   for(j in 1:2){
#     alpha[j] ~ dnorm(0, sd = 2)
#   }
#   # Derived parameter, the number of cells occupied
#   zsum <- sum(z[1:npixel])
# })
# 
# head(grid_sf$grid_id) # les ID de toutes les cellules
# 
# pixel.id.det <- grid_sf$grid_id[grid_sf$nnutria > 0] # les ID des cellules où il y a au moins une occurrence
# head(pixel.id.det)
# 
# npix <- nrow(grid_sf)
# s.area <- as.numeric(units::set_units(st_area(grid_sf)[1],"km^2"))
# (logarea <- log(s.area / npix))
# 
# data <- list(
#   cell_area = logarea,
#   x_s = (grid_sf$lgr_rivieres - mean(grid_sf$lgr_rivieres))/sd(grid_sf$lgr_rivieres), # distribution
#   h_s = (grid_sf$lgr_routes - mean(grid_sf$lgr_routes))/sd(grid_sf$lgr_routes), # detection
#   ones = rep(1, length(pixel.id.det)))
# 
# constants <- list(
#   npixel = npix,
#   m = length(pixel.id.det), # à modifier pour prendre en compte toutes les détections
#   CONSTANT = 10000,
#   po_pixel = pixel.id.det) # à modifier pour prendre en compte toutes les détections
# 
# zinit <- numeric(npix)
# zinit[pixel.id.det] <- 1
# inits <- function(){
#   list(
#     beta = rnorm(2, 0, 1), 
#     alpha = rnorm(2, 0, 1)
#     #z = zinit
#   )
# }
# 
# params <- c("alpha", "beta", "z")
# 
# # MCMC settings (pour tester...)
# nc <- 1 #2
# nburn <- 1000 #5000  
# ni <- nburn + 1000 #30000
# nt <- 1
# 
# start <- Sys.time()
# # run the model!
# out <- nimbleMCMC(
#   code = code, 
#   constants = constants, 
#   data = data, 
#   inits = inits(),
#   monitors = params, 
#   niter = ni, 
#   nburnin = nburn, 
#   nchains = nc, 
#   thin = nt)
# end <- Sys.time()
# end - start
# 
# MCMCsummary(out)
# 
# res <- out
# #res <- rbind(out$chain1, out$chain2)
# 
# 
# # select z
# mask <- str_detect(colnames(res), "z")
# res_z <- res[,mask]
# grid_sf$zestim <- apply(res_z, 2, median)
# grid_sf$zmoy <- apply(res_z, 2, mean)
# 
# # viz
# ggplot() +
#   geom_sf(data = grid_sf, lwd = 0.1, aes(fill = as_factor(zestim))) +
#   labs(fill = "Présence potentielle estimée du ragondin") + 
#   geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) + 
#   geom_sf(data = nutria) + 
#   theme_void()
# 
# ggplot() +
#   geom_sf(data = grid_sf, lwd = 0.1, aes(fill = zmoy)) +
#   labs(fill = "Présence potentielle estimée du ragondin") + 
#   scale_fill_viridis_c() + 
#   geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) + 
#   geom_sf(data = nutria) + 
#   theme_void()
# 
# 
# ###################################
# # Test avec multiples détections par cellules...
# ###################################
# 
# # problèmes : messages de présence de NA, estimation des paramètres bizarre (NA),
# # très variable d'une exécution à l'autre... 
# 
# # Bayesian version of the Koshkina (2017) model.
# code <- nimbleCode({
#   # latent-state model
#   for(pixel in 1:npixel){
#     # latent state linear predictor
#     #
#     # x_s  = covariates for latent state
#     # beta = latent state model regression coefficients
#     # cell_area = log area of grid cell 
#     #
#     
#     log(lambda[pixel]) <- beta[1] + beta[2] * x_s[pixel] + cell_area
#     # Species presence in a gridcell as a Bernoulli trial
#     z[pixel] ~ dbern(1 - exp(-lambda[pixel]))
#     # presence only thinning prob linear predictor
#     #
#     # h_s = covariates for thinning probability
#     # alpha  = presence-only data model regression coefficients
#     #
#     logit(b[pixel]) <-  alpha[1] + alpha[2] * h_s[pixel]
#   }
#   # The presence only data model.
#   #
#   # This part of the model just uses the
#   #  what we have calculated above (lambda
#   #  and b). The denominator of this likelihood
#   #  is actually a scalar so we can calculate it
#   #  outside of a for loop. Let's do that first.
#   #
#   # The presence_only data model denominator, which
#   #  is the thinned poisson process across the
#   #  whole region (divided by the total number of 
#   #  data points because it has to be 
#   #  evaluated for each data point).
#   # m is the number of presence-only data points
#   po_denominator <- inprod(lambda[1:npixel], b[1:npixel]) / m
#   #
#   # Loop through each presence-only data point
#   #  using Bernoulli one's trick. The numerator
#   #  is just the thinned poisson process for
#   #  the ith data point.
#   #  po_pixel denotes the grid cell of the ith presence only data point
#   for(po in 1:m){
#     ones[po] ~ dbern(
#       exp(
#         nobs_pixel[po]*log(lambda[po_pixel[po]] * b[po_pixel[po]]) - 
#           po_denominator)   # modif ici pour prendre en compte le nombre d'observations par pixel : nobs_pixel[po]
#       / CONSTANT) # attention, voir issue https://github.com/mfidino/integrated-occupancy-model/issues/1
#   } 
#   # Priors for latent state model
#   for(i in 1:2){
#     beta[i] ~ dnorm(0, sd = 2)
#   }
#   # Priors for presence-only data model
#   for(j in 1:2){
#     alpha[j] ~ dnorm(0, sd = 2)
#   }
#   # Derived parameter, the number of cells occupied
#   zsum <- sum(z[1:npixel])
# })
# 
# head(grid_sf$grid_id) # les ID de toutes les cellules
# 
# pixel.id.det <- grid_sf$grid_id[grid_sf$nnutria > 0] # les ID des cellules où il y a au moins une occurrence
# head(pixel.id.det)
# 
# npix <- nrow(grid_sf)
# s.area <- as.numeric(units::set_units(st_area(grid_sf)[1],"km^2"))
# (logarea <- log(s.area / npix))
# 
# data <- list(
#   cell_area = logarea,
#   x_s = (grid_sf$lgr_rivieres - mean(grid_sf$lgr_rivieres))/sd(grid_sf$lgr_rivieres), # distribution
#   h_s = (grid_sf$lgr_routes - mean(grid_sf$lgr_routes))/sd(grid_sf$lgr_routes), # detection
#   ones = rep(1, length(pixel.id.det)))
# 
# constants <- list(
#   npixel = npix,
#   m = length(pixel.id.det), 
#   CONSTANT = 50000,
#   po_pixel = pixel.id.det,
#   nobs_pixel = grid_sf$nnutria[grid_sf$nnutria > 0]) # modifié pour prendre en compte toutes les détections
# 
# zinit <- numeric(npix)
# zinit[pixel.id.det] <- 1
# inits <- function(){
#   list(
#     beta = rnorm(2, 0, 1), 
#     alpha = rnorm(2, 0, 1),
#     z = zinit # pourquoi mis en commentaire dans code initial ?
#   )
# }
# 
# params <- c("alpha", "beta", "z")
# 
# # MCMC settings (pour tester...)
# nc <- 2
# nburn <- 5000 
# ni <- nburn + 10000 #30000
# nt <- 1
# 
# start <- Sys.time()
# # run the model!
# out <- nimbleMCMC(
#   code = code, 
#   constants = constants, 
#   data = data, 
#   inits = inits(),
#   monitors = params, 
#   niter = ni, 
#   nburnin = nburn, 
#   nchains = nc, 
#   thin = nt)
# end <- Sys.time()
# end - start
# 
# MCMCsummary(out)
# 
# res <- rbind(out$chain1, out$chain2)
# 
# 
# # select z
# mask <- str_detect(colnames(res), "z")
# res_z <- res[,mask]
# grid_sf$zestim <- apply(res_z, 2, median)
# grid_sf$zmoy <- apply(res_z, 2, mean)
# 
# # viz
# ggplot() +
#   geom_sf(data = grid_sf, lwd = 0.1, aes(fill = as_factor(zestim))) +
#   labs(fill = "Présence potentielle estimée du ragondin") + 
#   geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) + 
#   geom_sf(data = nutria) + 
#   theme_void()
# 
# ggplot() +
#   geom_sf(data = grid_sf, lwd = 0.1, aes(fill = zmoy)) +
#   labs(fill = "Présence potentielle estimée du ragondin") + 
#   scale_fill_viridis_c() + 
#   geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) + 
#   geom_sf(data = nutria) + 
#   theme_void()
