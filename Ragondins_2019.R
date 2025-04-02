#######################
# Ragondins 2019 uniquement
#######################

# librairies utiles

library(nimble)
library(MCMCvis)
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

rm(dpts_occitanie)

#######################
# Ragondins Occitanie 2019

# Import des données
nutria <- st_read("Data/Data_pts_test_infos.shp") %>%
  mutate(year = year(as.Date(jourdebut))) %>% 
  filter(year == 2019) %>%
  filter(nomvern == "Ragondin") 

# changement de format
occitanie <- occitanie %>%
  st_transform(crs = st_crs(nutria))

# intersection avec occitanie
nutria <- st_intersection(nutria, occitanie)

# plot des données 
p <- ggplot() +
  geom_sf(data = occitanie, fill = "white", color = "black", lwd = .5) + 
  geom_sf(data = nutria) + 
  labs(color = "") +
  theme_void()
p
# ggsave(plot = p, "fig/map.png", dpi = 600)

#####################
# Création de la grille
# n = 25
# grid <- st_make_grid(occitanie, n = n, what = "polygons", square = FALSE)

grid_area <- units::set_units(5000000,"m^2") #5km2
grid <- st_make_grid(occitanie, cellsize = grid_area, what = "polygons", square = FALSE)

# sf
grid_sf <- st_sf(grid)
rm(grid)

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

# Surface exacte des cellules en occitanie
grid_sf$area <- st_area(st_intersection(grid_sf, occitanie))

# # intersections nulles, faut-il les enlever ?
# ggplot() +
#   geom_sf(data = grid_sf, lwd = 0.1, aes(fill = as.numeric((as.numeric(area)==0)))) +
#   geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
#   theme_void()

######################
# Ragondins par cellules

# Nombre d'observations par cellule
# https://gis.stackexchange.com/questions/323698/counting-points-in-polygons-with-sf-package-of-r
grid_sf$nnutria <- lengths(st_intersects(grid_sf, nutria))

table(grid_sf$nnutria)

# cellules ayant au moins une observation
nutria_count <- filter(grid_sf, nnutria > 0)

# plot
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
# elles doivent toutes être intersectées avec Occitanie
# et ramenées à des variables par unité d'aire

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

# avec la grille (2019 uniquement)
temp_min <- terra::extract(toccitanie_min, grid_sf # exact = T, # pour toutes les cellules...
    ) %>%
  group_by(ID) %>%
  summarise(tmin2019 = mean(lyr.17, na.rm = T))

ggplot() +
  geom_sf(data = grid_sf, lwd = 0.1, color = NA, aes(fill = temp_min$tmin2019)) +
  scale_fill_viridis_c() +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()

temp_max <- terra::extract(toccitanie_max, grid_sf, # exact = T, # pour toutes les cellules... 
    ) %>%
  group_by(ID) %>%
  summarise(tmax2019 = mean(lyr.17, na.rm = T))

temp_mean <- terra::extract(toccitanie_mean, grid_sf,  # exact = T, # pour toutes les cellules... 
    ) %>%
  group_by(ID) %>%
  summarise(tmean2019 = mean(lyr.17, na.rm = T))

# ajout des covariables dans grif_sf
grid_sf$temp_min <- temp_min$tmin2019
grid_sf$temp_max <- temp_max$tmax2019
grid_sf$temp_mean <- temp_mean$tmean2019

# save(grid_sf, file = "grid_sf_5km2.RData")

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

# avec la grille (2019 uniquement)
prec_cum <- terra::extract(poccitanie_cum, grid_sf) %>%
  group_by(ID) %>%
  summarise(pcum2019 = mean(lyr.17, na.rm = T)
  )

ggplot() +
  geom_sf(data = grid_sf, lwd = 0.1, color = NA, aes(fill = prec_cum$pcum2019)) +
  scale_fill_viridis_c() +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()

# Ajout de la covariable à grid_sf
grid_sf$prec_cum <- prec_cum$pcum2019

# save(grid_sf, file = "grid_sf_5km2.RData")

######################
# Variables de sol

# Données de CORINE révisées 2012 (faut-il tenir compte de l'évolution ?)
Cori12 <- st_read("Data/CLC12_RLRMP_RGF.shp")
Cori12 <- Cori12 %>%
  st_transform(crs = st_crs(grid_sf)) %>%
  st_intersection(occitanie) 

str(Cori12)
unique(Cori12$CODE_12)

# Surface agricole 
agri <- Cori12 %>%
  filter(str_sub(CODE_12,1,1)=="2")

# Avec la grille
grid_agri <- agri %>%
  st_intersection(grid_sf) %>%
  mutate(area = st_area(.)) %>%
  group_by(grid_id) %>%
  summarise(aera_agri = sum(area)) %>%
  as_tibble() %>%
  select(-geometry)

# Proportion de terres agricoles
grid_sf <- grid_sf %>% 
  full_join(grid_agri, by = "grid_id") %>%
  mutate(.before = 1, agri_cover = aera_agri/area) %>%
  select(-aera_agri)

ggplot() +
  geom_sf(data = grid_sf, lwd = 0.1, color = NA, aes(fill = as.numeric(agri_cover))) +
  scale_fill_viridis_c(
    labels = scales::percent_format()
  ) +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()

# Cellules manquantes ! 
# save(grid_sf, file = "grid_sf_5km2.RData")


######################
# Densité de population 

pop <- st_read("Data/pop2021.gpkg")
pop <- pop %>% 
  st_transform(crs = st_crs(grid_sf)) %>%
  st_intersection(occitanie)

ggplot() +
  geom_sf(data = pop, lwd = 0.1, color = NA, aes(fill = TOT_P_2021)) +
  scale_fill_viridis_c() +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()

# la frontière n'est pas donnée !

# Avec la grille
grid_pop <- pop %>%
  st_intersection(grid_sf) %>%
  group_by(grid_id) %>%
  summarise(hab = sum(TOT_P_2021)) %>%
  as_tibble() %>%
  select(-geom)

# Log densité
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

# save(grid_sf, file = "Data/grid_sf_5km2.RData")


######################
# Longueur de chemins

# OpenStreetMap https://www.data.gouv.fr/fr/datasets/itineraires-de-randonnee-dans-openstreetmap/
chemins_osm <- st_read("Data/hiking_foot_routes_lineLine.shp")
chemins_osm <- chemins_osm %>% 
  st_transform(crs = st_crs(grid_sf)) %>%
  st_intersection(occitanie)

# Avec la grille
grid_chemin <- chemins_osm  %>%
  st_intersection(grid_sf) %>%
  mutate(len = st_length(.)) %>%
  group_by(grid_id) %>%
  summarise(lgr = sum(len)) %>%
  as_tibble() %>%
  select(-geometry)

# longueur/surface 
grid_sf <- grid_sf %>% 
  full_join(grid_chemin, by = "grid_id") %>%
  mutate(.before = 1, lgr_chemins = lgr/area) %>%
  select(-lgr)


ggplot() +
  geom_sf(data = grid_sf, lwd = 0.1, aes(fill = as.numeric(lgr_chemins))) +
  labs(fill = "Longueur de chemins (km/m^2)") +
  scale_fill_viridis_c() +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()

# save(grid_sf, file = "Data/grid_sf_5km2.RData")

######################
# Longueur de routes

roads <- st_read("Data/TRONCON_ROUTE.shp")
str(roads)
unique(roads$CLASS_ADM)

routes_occ <- roads %>%
  st_transform(crs = st_crs(grid_sf)) %>%
  st_intersection(occitanie) %>%
  filter(CLASS_ADM != "Autoroute" & CLASS_ADM !="Sans objet") # D\xe9partementale ne marche pas...

# Avec la grille
grid_route <- routes_occ  %>%
  st_intersection(grid_sf) %>%
  mutate(len = st_length(.)) %>%
  group_by(grid_id) %>%
  summarise(lgr = sum(len)) %>%
  as_tibble() %>%
  select(-geometry)

# longueur/surface 
grid_sf <- grid_sf %>% 
  full_join(grid_route, by = "grid_id") %>%
  mutate(.before = 1, lgr_routes = lgr/area) %>%
  select(-lgr)

ggplot() +
  geom_sf(data = grid_sf, lwd = 0.1, aes(fill = as.numeric(lgr_routes))) +
  labs(fill = "Longueur de routes (km/m^2)") +
  scale_fill_viridis_c() +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()

# save(grid_sf, file = "Data/grid_sf_5km2.RData")

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
  st_intersection(occitanie)

# Avec la grille
grid_rivieres <- rivers_occ  %>%
  st_intersection(grid_sf) %>%
  mutate(len = st_length(.)) %>%
  group_by(grid_id) %>%
  summarise(lgr = sum(len)) %>%
  as_tibble() %>%
  select(-geometry)

# longueur/surface 
grid_sf <- grid_sf %>% 
  full_join(grid_rivieres, by = "grid_id") %>%
  mutate(.before = 1, lgr_rivieres = lgr/area) %>%
  select(-lgr)

# visualisation
ggplot() +
  geom_sf(data = grid_sf, lwd = 0.1, aes(fill = as.numeric(lgr_rivieres))) +
  labs(fill = "Longueur de cours d'eau par unité d'aire (km/m^2)") +
  scale_fill_viridis_c() +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()

# save(grid_sf, file = "Data/grid_sf_5km2.RData")

######################
# Distance aux rivières

dist = 10000 # rayon de recherche autour de la cellule
dist_buff <- function (i) {
  cell = grid_sf[i,]
  buffer <- st_buffer(cell, dist) 
  plan_filtre <- river_lines[st_intersects(river_lines, buffer, sparse = F),]
  return(min(st_distance(cell, plan_filtre)))
}
# renvoit des Inf si pas d'intersection

for (i in 1:nrow(grid_sf)) {
  grid_sf$dist_rivieres[i] <- dist_buff(i)
}

ggplot() +
  geom_sf(data = grid_sf, color = NA, aes(fill = dist_rivieres)) +
  scale_fill_viridis_c(
  ) +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()

# save(grid_sf, file = "RData/grid_sf_5km2.RData")

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

# plans d'eau d'Occitanie
plan_eau_occ <- plan_eau %>%
  st_transform(crs = st_crs(grid_sf)) %>%
  st_intersection(occitanie)

# Avec la grille
grid_plan <-  plan_eau_occ  %>%
  st_intersection(grid_sf) 

grid_plan_union <- grid_plan %>% 
  group_by(grid_id) %>% 
  summarise(geometry = st_union(geometry)) %>% # pour enlever les chevauchements étranges...
  mutate(area_eau_tot = st_area(.)) %>%
  as_tibble() %>%
  select(-geometry)

grid_sf <- grid_sf %>% 
  full_join(grid_plan_union, by = "grid_id") %>%
  mutate(.before = 1, surface_en_eau = area_eau_tot/area) %>%
  select(-area_eau_tot)

ggplot() +
  geom_sf(data = grid_sf, lwd = 0.1, aes(fill = as.numeric(surface_en_eau))) +
  scale_fill_viridis_c(
    labels = scales::percent_format()
  ) +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()

# save(grid_sf, file = "grid_sf_5km2.RData")

######################
# Distance aux plans d'eau

# # distance aux plans d'eau
# dist_eau <- st_distance(grid_sf, st_transform(plan_eau, crs = st_crs(grid_sf)))
# # save(dist_eau, file = "RData/dist_eau.RData")
# 
# # On regarde la distance minimale
# grid_sf$dist_plan_eau <- apply(dist_eau, 1, min)

# Pour essayer d'être plus efficace en utilisant un buffer
# # format -> utile ?
# plan_eau <- plan_eau %>%
#   st_transform(crs = st_crs(grid_sf))

dist = 30000 # rayon de recherche autour de la cellule
dist_buff <- function (i) {
  cell = grid_sf[i,]
  buffer <- st_buffer(cell, dist) 
  plan_filtre <- plan_eau[st_intersects(plan_eau, buffer, sparse = F),]
  return(min(st_distance(cell, plan_filtre)))
}
# renvoit des Inf si pas d'intersection

for (i in 1:nrow(grid_sf)) {
  grid_sf$dist_plan_eau[i] <- dist_buff(i)
}

ggplot() +
  geom_sf(data = grid_sf, color = NA, aes(fill = dist_plan_eau)) +
  scale_fill_viridis_c(
  ) +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()

# save(grid_sf, file = "RData/grid_sf_5km2.RData")

######################
# Observations toutes espèces (GBIF)

library(rgbif)

# Récupération des clés taxonomiques des groupes d'intérêt
oiseaux_key <- name_backbone(name = "Aves")$usageKey
mammiferes_key <- name_backbone(name = "Mammalia")$usageKey
amphibiens_key <- name_backbone(name = "Amphibia", rank = "class")$usageKey
reptiles_key <- name_backbone(name = "Reptilia")$usageKey

# Affichage des clés
print(c(oiseaux_key, mammiferes_key, amphibiens_key, reptiles_key))

# Définition de la zone géographique
bbox_wkt <- "POLYGON((-0.4 42, 5 42, 5 45.1, -0.4 45.1, -0.4 42))"

# Fonction pour récupérer les données d'un taxon donné
get_gbif_data <- function(taxon_key, month) {
  occ_search(
    classKey = taxon_key, 
    year = 2019,
    month = month, 
    hasCoordinate = TRUE,  # Seulement les données géolocalisées
    limit = 500,  # Maximum d'enregistrements récupérés (à adapter)
    geometry = paste(bbox_wkt, collapse = ",")  # Zone Occitanie
  )$data
}

# on télécharge 500 données mensuelles pour chaque taxon
# on fait varier les mois pour éviter la saisonnalité... 

month = 1
data_oiseaux <- get_gbif_data(oiseaux_key, month)
data_mammiferes <- get_gbif_data(mammiferes_key, month)
data_amphibiens <- get_gbif_data(amphibiens_key, month)
data_reptiles <- get_gbif_data(reptiles_key, month)
gbif_data <- bind_rows(data_oiseaux, data_mammiferes, data_amphibiens, data_reptiles)

for (month in 2:12) {
  data_oiseaux <- get_gbif_data(oiseaux_key, month)
  data_mammiferes <- get_gbif_data(mammiferes_key, month)
  data_amphibiens <- get_gbif_data(amphibiens_key, month)
  data_reptiles <- get_gbif_data(reptiles_key, month)
  gbif_data <- bind_rows(gbif_data, data_oiseaux, data_mammiferes, data_amphibiens, data_reptiles)
}
# Pas de réptiles...
# save(gbif_data, file = "gbif_data.RData")

table(gbif_data$month)

dim(gbif_data)  # Nombre de lignes et colonnes
glimpse(gbif_data)  # Aperçu des colonnes disponibles

# On garde uniquement les colonnes utiles pour l'analyse :
gbif_clean <- gbif_data %>%
  select(species, class, order, family, decimalLatitude, decimalLongitude, year, basisOfRecord)

head(gbif_clean)

# Transformer en objet spatial
gbif_sf <- st_as_sf(gbif_clean, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) 

gbif_sf <- gbif_sf %>%
  st_transform(crs = st_crs(occitanie))

# Afficher sur une carte
ggplot() +
  geom_sf(data = occitanie, fill="white", color = "black", lwd = .5) + 
  geom_sf(data = gbif_sf, aes(color = class), alpha = 0.6) +
  theme_minimal() +
  labs(title = "Observations GBIF en Occitanie (2019)", color = "Classe") +
  theme_void()

# En occitanie
gbif_occ <-  gbif_sf %>% st_intersection(occitanie)

# Plot
ggplot() +
  geom_sf(data = occitanie, fill = "white", color = "black", lwd = .5) + 
  geom_sf(data = gbif_occ, aes(color = class), alpha = 0.6) +
  theme_minimal() +
  labs(title = "Observations GBIF en Occitanie (2019)", color = "Classe") +
  theme_void()

# Nombre d'observations par cellule
nobs_gbif <- lengths(st_intersects(grid_sf, gbif_occ))

# sur une carte
ggplot() +
  geom_sf(data = grid_sf, lwd = 0.1, aes(fill = as.numeric(nobs_gbif))) +
  scale_fill_viridis_c() +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()

hist(nobs_gbif)

# avec une troncature à 50 (voir Twinings)
# Troncature
nobs_gbif[nobs_gbif>=50] <- 50

# sur une carte
ggplot() +
  geom_sf(data = grid_sf, lwd = 0.1, aes(fill = as.numeric(nobs_gbif))) +
  scale_fill_viridis_c() +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()

# Ajout de la covariable (variable par unité d'aire)
grid_sf$GBIF <- nobs_gbif/grid_sf$area 

# sur une carte
ggplot() +
  geom_sf(data = grid_sf, lwd = 0.1, aes(fill = as.numeric(GBIF))) +
  scale_fill_viridis_c() +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()

# save(grid_sf, file = "grid_sf_5km2.RData")


