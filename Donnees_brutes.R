################################################################################
#            Chargement des données en Occitanie et à proximité                #
################################################################################

# librairies utiles ------------------------------------------------------------

library(ggplot2)
library(tidyverse)
library(sf)
library(tidyterra)
library(lubridate)
library(KrigR)
library(rgbif)
library(osmdata)
library(terra)


# Occitanie et données de ragondins --------------------------------------------

## Occitanie ###############

# contours des départements d'Occitanie
dpts_occitanie <- st_read("Data/departements-d-occitanie.shp") 

# contours de la région
occitanie <- dpts_occitanie %>% st_union()

rm(dpts_occitanie)


## Ragondins #############

periode = (2010:2024)

# Import des données
nutria_periode <- st_read("Data/CEN_2025/Ragondin_rat_musque_pts_2025.shp") %>%
  mutate(year = year(as.Date(DateDebut))) %>% 
  filter(year %in% periode) %>%
  filter(NomVernacu == "Ragondin") 

table(nutria_periode$year)

# changement de format
occitanie <- occitanie %>%
  st_transform(crs = st_crs(nutria_periode))

# Les ragondins d'Occitanie
nutria_periode <- st_intersection(nutria_periode, occitanie)

table(nutria_periode$year)

# plot des données 
ggplot() +
  geom_sf(data = occitanie, fill = "white", color = "black", lwd = .5) + 
  geom_sf(data = nutria_periode, color = "chartreuse4", cex = 0.5) + 
  facet_wrap(~ year, nrow = 4) +
  labs(color = "") +
  theme_void()

# ggsave(plot = p, "fig/map.png", dpi = 600)


## Buffers ###################

# buffer autour de l'Occitanie (10km)
occitanie_buff <- st_buffer(occitanie, dist = 10000, endCapStyle = "ROUND")

# buffer côté espagnol (10km)
fr <- st_read("Data/regions_2015_metropole_region.shp", options = "ENCODING=ISO-8859-15") %>% 
  st_union() %>% 
  st_transform(crs = st_crs(occitanie))

espagne_buff <- st_difference(occitanie_buff, fr) 

box_espagne <- st_bbox(c(xmin = -0.49, ymin = 42.19, xmax = 3.33, ymax = 42.87), crs = 4326) %>%
  st_as_sfc() %>%
  st_transform(crs = st_crs(occitanie))

espagne_buff <- st_crop(espagne_buff, box_espagne) # on recadre un peu

ggplot() +
  geom_sf(data = espagne_buff, fill = "grey") +
  geom_sf(data = occitanie, fill = NA)

rm(fr)


# Variables explicatives -------------------------------------------------------


## Températures ####################

Dir.Base <- getwd() # identifying the current directory
Dir.Data <- file.path(Dir.Base, "Data") # folder path for data

# contours de la région
loc_site <- st_transform(occitanie, crs = 4326)

st_bbox(loc_site)

occitanie_raster <- st_bbox(c(xmin = -0.4, 
                              xmax = 5, 
                              ymax = 45.1, 
                              ymin = 42), 
                            crs = st_crs(4326)) %>%
  st_as_sfc() %>%
  terra::vect() %>%
  terra::rast()

toccitanie <- CDownloadS(
  Variable = "2m_temperature",
  DataSet = "reanalysis-era5-land-monthly-means",
  Type = "monthly_averaged_reanalysis",
  DateStart = "2010-01-01 00:00",
  DateStop = "2024-12-31 23:00",
  TZone = "Europe/Paris",
  TResolution = "month",
  TStep = 1,
  Extent = occitanie_raster, 
  Dir = Dir.Data,
  FileName = "Toccitanie",
  API_User = "camille.mottier@cefe.cnrs.fr",
  API_Key = "4700337b-89cd-42e8-98b6-a450f702d332"
)

toccitanie <- weathermetrics::kelvin.to.celsius(toccitanie)

# variables annuelles sur la période
# mois min
toccitanie_min <- terra::aggregate(toccitanie, fact=c(1,1,12), fun=min)
# mois max 
toccitanie_max <- terra::aggregate(toccitanie, fact=c(1,1,12), fun=max)
# mois mean
toccitanie_mean <- terra::aggregate(toccitanie, fact=c(1,1,12), fun=mean)

# moyenne des voisins pour cellules manquantes 
weight_matrix <- matrix(1, nrow=3, ncol=3)
weight_matrix[2,2] <- NA
toccitanie_min <- terra::focal(toccitanie_min, w=weight_matrix, fun = mean, na.policy = "only", na.rm = T)
toccitanie_max <- terra::focal(toccitanie_max, w=weight_matrix, fun = mean, na.policy = "only", na.rm = T)
toccitanie_mean <- terra::focal(toccitanie_mean, w=weight_matrix, fun = mean, na.policy = "only", na.rm = T)

# plot 
ggplot() +
  geom_spatraster(data=toccitanie_min) +
  # geom_spatraster(data=toccitanie_max) +
  # geom_spatraster(data=toccitanie_mean) +
  facet_wrap(~ lyr, nrow = 4) +
  scale_fill_viridis_c() +
  geom_sf(data = loc_site, 
          fill = NA)

rm(toccitanie)

## Précipitations ####################

poccitanie <- CDownloadS(
  Variable = "total_precipitation",
  DataSet = "reanalysis-era5-land-monthly-means",
  CumulVar = TRUE,
  Type = "monthly_averaged_reanalysis",
  DateStart = "2010-01-01 00:00",
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

# Cumuls annuels
poccitanie_cum <- terra::aggregate(poccitanie, fact=c(1,1,12), fun=sum)
# Moyenne pour les cellules manquantes
poccitanie_cum <- terra::focal(
  poccitanie_cum,
  w = weight_matrix,
  fun = mean,
  na.policy = "only",
  na.rm = T
)

# plot 
ggplot() +
  geom_spatraster(data=poccitanie_cum) +
  facet_wrap(~ lyr, nrow = 4) +
  scale_fill_viridis_c() +
  geom_sf(data = loc_site, 
          fill = NA)

rm(poccitanie, Dir.Base, weight_matrix, Dir.Data, loc_site, occitanie_raster)


## Surface agricole ####################

# Données de CORINE révisées 2012 (faut-il tenir compte de l'évolution ?)
Cori12 <- st_read("Data/CLC12_RLRMP_RGF.shp")
Cori12 <- Cori12 %>%
  st_transform(crs = st_crs(occitanie)) %>%
  st_intersection(occitanie) 

str(Cori12)
unique(Cori12$CODE_12)

# Surface agricole 
agri <- Cori12 %>%
  filter(str_sub(CODE_12,1,1)=="2")

rm(Cori12)


## Densité de population ####################

pop <- st_read("Data/Densite/popoccitanie.shp") # Occitanie uniquement
pop <- pop %>% 
  st_transform(crs = st_crs(occitanie)) 

ggplot() +
  geom_sf(data = pop, lwd = 0.1, color = NA, aes(fill = TOT_P_2021)) +
  scale_fill_viridis_c() +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()



## Distances aux chemins ####################

# Open Street Map
# Pour les labels à utiliser, voir : 
# https://wiki.openstreetmap.org/wiki/Map_features#Highway

# Côté français
# https://www.data.gouv.fr/fr/datasets/itineraires-de-randonnee-dans-openstreetmap/
chemins_osm <- st_read("Data/Chemins/hiking_foot_routes_lineLine.shp")
chemins_osm <- chemins_osm %>% 
  st_transform(crs = st_crs(occitanie)) %>%
  st_intersection(occitanie_buff)

ggplot() +
  geom_sf(data = chemins_osm, color = "royalblue") +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5)

# Côté espagnol
# Créer la requête OSM pour extraire les chemins de randonnée (frontière uniquement !)
query <- opq(bbox = st_transform(box_espagne, crs = 4326)) %>%
  add_osm_feature(key = "route", value = c("foot", "hiking"))

# Exécuter la requête et récupérer les résultats
data_osm <- osmdata_sf(query)

# Extraire les chemins de randonnée (les autres formats donnent la même chose)
lines_osm <- data_osm$osm_lines
# points_osm <- data_osm$osm_points
# poly_osm <- data_osm$osm_polygons
# multilines_osm <- data_osm$osm_multilines

# On intersecte avec le buffer
lines_osm <- lines_osm %>% 
  st_transform(crs = st_crs(occitanie)) %>%
  st_intersection(espagne_buff)

# save(lines_osm, file = "RData/chemins_espagne.RData")

ggplot() +
  geom_sf(data = lines_osm, color = "royalblue") +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5)

# On réunit Espagne et France 
chemins <- rbind(
  chemins_osm %>% select(osm_id, geometry), 
  lines_osm %>% select(osm_id, geometry)
)

ggplot() +
  geom_sf(data = chemins, color = "royalblue") +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5)

rm(chemins_osm, lines_osm, query, data_osm)



## Distances aux routes ####################

# Côté français
roads <- st_read("Data/TRONCON_ROUTE.shp") # routes de la France entière
str(roads)
unique(roads$CLASS_ADM)

# Dans le buffer uniquement
routes_occ <- roads %>%
  st_transform(crs = st_crs(occitanie)) %>%
  st_intersection(occitanie_buff) %>%
  filter(CLASS_ADM != "Autoroute" & CLASS_ADM !="Sans objet") # D\xe9partementale ne marche pas...

ggplot() +
  geom_sf(data = routes_occ, color = "royalblue") +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) 

# Côté espagnol
# Créer la requête OSM pour extraire les routes (frontière uniquement !)
query <- opq(bbox = st_transform(box_espagne, crs = 4326)) %>%
  add_osm_feature(key = "route", value = c("road"))

# Exécuter la requête et récupérer les résultats
road_osm <- osmdata_sf(query)

# Les routes du buffer (sauf autoroute)
roads_lines <- road_osm$osm_lines %>%
  st_transform(crs = st_crs(occitanie)) %>%
  filter(highway %in% c("trunk", "primary", "secondary", "tertiary", "unclassified")) %>%
  st_intersection(espagne_buff)

ggplot() +
  geom_sf(data = roads_lines, color = "royalblue") +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) # pas très convaincant

# On réunit espagne et france
routes <- rbind(
  roads_lines %>% select(geometry), 
  routes_occ %>% select(geometry)
)

ggplot() +
  geom_sf(data = routes, color = "royalblue") +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5)


rm(query, road_osm, roads, roads_lines, routes_occ)



## Distances aux cours d'eau ####################

# Côté français
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
  rbind(TarnetGaronne) 

rm(list=c("Ariege", "Aude","Aveyron","Gard","Gers","HauteGaronne",
          "HautesPyrenees","Herault","Lot","Lozere","PyreneesOrientales","Tarn","TarnetGaronne"))

# on ajoute les départements du buffer
for (i in 1:10) {
  river_lines <- river_lines %>%
    rbind(st_read(paste("Data/Rivieres/Buffer/Buffer", i,  "/COURS_D_EAU.shp", sep = "")))
}

# on garde Occitanie + buffer
river_lines <- river_lines %>%
  st_transform(crs = st_crs(occitanie)) %>%
  st_intersection(occitanie_buff) %>%
  st_simplify()

ggplot() +
  geom_sf(data = river_lines, lwd = 0.1) +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()

# Côté espagnol

river_esp <- st_read("Data/Rivieres/Espagne/HydroRIVERS_v10_eu.shp")

# On garde uniquement le buffer autour de la frontiere
river_esp <- river_esp %>%
  st_transform(crs = st_crs(occitanie)) %>%
  st_crop(box_espagne) %>% # on garde que la box
  st_intersection(espagne_buff) # Espagne uniquement

ggplot() +
  geom_sf(data = rivieres_espagne, color = "royalblue") +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5)

# On réunit espagne et france (doublons possibles)
rivieres <- rbind(
  river_lines %>% select(geometry), 
  river_esp %>% select(geometry)
)


rm(river_esp, river_lines)


## Distance aux plans d'eau ####################

# Côté fançais

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
  rbind(TarnetGaronne) 

rm(list=c("Ariege", "Aude","Aveyron","Gard","Gers","HauteGaronne",
          "HautesPyrenees","Herault","Lot","Lozere","PyreneesOrientales","Tarn","TarnetGaronne"))

# on ajoute les départements du buffer
for (i in 1:10) {
  plan_eau <- plan_eau %>%
    rbind(st_read(paste("Data/Rivieres/Buffer/Buffer", i,  "/PLAN_D_EAU.shp", sep = "")))
}

# on garde le buffer
plan_eau <- plan_eau %>%
  st_transform(crs = st_crs(occitanie)) %>%
  st_intersection(occitanie_buff) %>%
  st_simplify()

# Côté espagnol

plan_eau_esp <- st_read("Data/Rivieres/Espagne/HydroLAKES_polys_v10.shp")

# On garde uniquement le buffer autour de la frontiere
plan_eau_esp <- plan_eau_esp %>%
  st_transform(crs = st_crs(occitanie)) %>%
  st_crop(box_espagne) %>% # on garde que la box
  st_intersection(espagne_buff) # Espagne uniquement

ggplot() +
  geom_sf(data = plan_eau_esp, color = "royalblue") +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5)

# On réunit espagne et france (doublons possibles)
plan_eau_tot <- rbind(
  plan_eau %>% select(geometry), 
  plan_eau_esp %>% select(geometry)
)

rm(plan_eau, plan_eau_esp)


## Observations toutes espèces (GBIF) ####################

# /!\ corrigé ! Plus de données chargées, dans une région plus adaptée

# Récupération des clés taxonomiques des groupes d'intérêt
oiseaux_key <- name_backbone(name = "Aves")$usageKey
mammiferes_key <- name_backbone(name = "Mammalia")$usageKey
amphibiens_key <- name_backbone(name = "Amphibia", rank = "class")$usageKey
reptiles_key <- name_backbone(name = "Reptilia")$usageKey

# Affichage des clés
print(c(oiseaux_key, mammiferes_key, amphibiens_key, reptiles_key))

# # Définition de la zone géographique
# bbox_wkt <- "POLYGON((-0.4 42, 5 42, 5 45.1, -0.4 45.1, -0.4 42))"

# Définition de la zone géographique (enveloppe convexe d'occitanie)
bbox_wkt <- occitanie %>% 
  st_transform(crs = "WGS84") %>%
  st_convex_hull() %>%
  st_as_text()

# Fonction pour récupérer les données d'un taxon donné
get_gbif_data <- function(taxon_key, month, year) {
  occ_search(
    classKey = taxon_key, 
    year = year,
    month = month, 
    hasCoordinate = TRUE,  # Seulement les données géolocalisées
    limit = 2000,  # Maximum d'enregistrements récupérés 
    geometry = bbox_wkt  # Zone Occitanie
  )$data
}

# on télécharge 2000 données mensuelles pour chaque taxon et chaque année
# on fait varier les mois pour éviter la saisonnalité... 

gbif_data <- NULL

for (annee in periode) {
  print(annee)
  for (mois in 1:12) {
    data_oiseaux <- get_gbif_data(oiseaux_key, mois, annee)
    data_mammiferes <- get_gbif_data(mammiferes_key, mois, annee)
    data_amphibiens <- get_gbif_data(amphibiens_key, mois, annee)
    # data_reptiles <- get_gbif_data(reptiles_key, mois, annee) # pas de reptiles sur la periode
    gbif_data <- bind_rows(gbif_data,
                           data_oiseaux,
                           data_mammiferes,
                           data_amphibiens)
    # data_reptiles
  }
}

# save(gbif_data, file = "RData/gbif_data_periode.RData")

table(gbif_data$year)
table(gbif_occ$year, gbif_occ$month)
table(gbif_occ$year, gbif_occ$class)

glimpse(gbif_data)  # Aperçu des colonnes disponibles

# On garde uniquement les colonnes utiles pour l'analyse :
gbif_data <- gbif_data %>%
  select(species, class, order, family, decimalLatitude, decimalLongitude, year, basisOfRecord)

head(gbif_data)

# Transformer en objet spatial
gbif_data <- st_as_sf(gbif_data, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) 

gbif_data <- gbif_data %>%
  st_transform(crs = st_crs(occitanie))

# En occitanie
gbif_occ <- gbif_data[st_intersects(gbif_data, occitanie, sparse = FALSE), ]

# Plot
ggplot() +
  geom_sf(data = occitanie, fill = "white", color = "black", lwd = .5) + 
  geom_sf(data = gbif_occ, aes(color = class), alpha = 0.6) +
  theme_minimal() +
  facet_wrap(~year, nrow = 4) +
  labs(color = "Classe") +
  theme_void()

rm(gbif_data)

# Nettoyage des données inutiles et sauvegarde ---------------------------------

rm(box_espagne, espagne_buff, occitanie_buff, i)
writeRaster(toccitanie_min, filename = "Data/A_charger/toccitanie_min.tif")
writeRaster(toccitanie_max, filename = "Data/A_charger/toccitanie_max.tif")
writeRaster(toccitanie_mean, filename = "Data/A_charger/toccitanie_mean.tif")
writeRaster(poccitanie_cum, filename = "Data/A_charger/poccitanie_cum.tif")

save(agri, chemins, gbif_occ, nutria_periode, occitanie, plan_eau_tot, pop, rivieres, routes, periode, file = "Data/A_charger/donnees_a_charger.RData")
