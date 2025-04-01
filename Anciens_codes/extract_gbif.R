library(rgbif)
library(tidyverse)
library(sf)

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
    # taxonKey = taxon_key,
    classKey = taxon_key, # proposition de remplacement ? 
    year = 2019,
    month = month, 
    hasCoordinate = TRUE,  # Seulement les données géolocalisées
    limit = 500,  # !!!!!!!!!!!!!!! Maximum d'enregistrements récupérés (à adapter)
    geometry = paste(bbox_wkt, collapse = ",")  # Zone Occitanie
  )$data
}

# on fait varier les mois pour éviter la saisonnalité... 

month = 1
data_oiseaux <- get_gbif_data(oiseaux_key, month)
data_mammiferes <- get_gbif_data(mammiferes_key, month)
data_amphibiens <- get_gbif_data(amphibiens_key, month)
data_reptiles <- get_gbif_data(reptiles_key, month)
gbif_data <- bind_rows(data_oiseaux, data_mammiferes, data_amphibiens, data_reptiles)


# Télécharger les données pour chaque groupe
for (month in 2:12) {
  data_oiseaux <- get_gbif_data(oiseaux_key, month)
  data_mammiferes <- get_gbif_data(mammiferes_key, month)
  data_amphibiens <- get_gbif_data(amphibiens_key, month)
  data_reptiles <- get_gbif_data(reptiles_key, month)
  gbif_data <- bind_rows(gbif_data, data_oiseaux, data_mammiferes, data_amphibiens, data_reptiles)
}

table(gbif_data$month)

# # ragondins ?
# data_ragondins <- occ_search(
#     taxonKey = 4264680,
#     year = 2019,
#     hasCoordinate = TRUE,  # Seulement les données géolocalisées
#     limit = 5000,  # !!!!!!!!!!!!!!! Maximum d'enregistrements récupérés (à adapter)
#     geometry = paste(bbox_wkt, collapse = ",")  # Zone Occitanie
#   )$data

dim(gbif_data)  # Nombre de lignes et colonnes
glimpse(gbif_data)  # Aperçu des colonnes disponibles

# On garde uniquement les colonnes utiles pour l'analyse :
gbif_clean <- gbif_data %>%
  select(species, class, order, family, decimalLatitude, decimalLongitude, year, basisOfRecord)
head(gbif_clean)

# Transformer en objet spatial
gbif_sf <- st_as_sf(gbif_clean, coords = c("decimalLongitude", "decimalLatitude"), crs = st_crs(occitanie)) 
save(gbif_sf, file = "gbif_sf.RData")

# Afficher sur une carte
ggplot() +
  geom_sf(data = occitanie, fill="white", color = "black", lwd = .5) + 
  geom_sf(data = gbif_sf, aes(color = class), alpha = 0.6) +
  theme_minimal() +
  labs(title = "Observations GBIF en Occitanie (2019)", color = "Classe") +
  theme_void()

# Avec la grille
# Question : faut-il intersecter avec grid ou occitanie ??? (cellules de la frontière et en bord de mer !)
grid_gbif <- st_transform(gbif_sf, crs = st_crs(grid_sf)) %>%
  st_intersection(st_union(grid_sf))

# dans l'Occitanie
ggplot() +
  geom_sf(data = grid_sf, fill = "white", color = "black", lwd = .5) + 
  geom_sf(data = grid_gbif, aes(color = class), alpha = 0.6) +
  theme_minimal() +
  labs(title = "Observations GBIF en Occitanie (2019)", color = "Classe") +
  theme_void()

# Nombre d'observations par cellule
nobs_gbif <- lengths(st_intersects(grid_sf, grid_gbif))
grid_sf$nobs_gbif <- nobs_gbif/st_area(grid_sf) 

# sur une carte
ggplot() +
  geom_sf(data = grid_sf, lwd = 0.1, aes(fill = as.numeric(nobs_gbif))) +
  scale_fill_viridis_c() +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()

# avec une troncature à 150 (voir Twinings)
# Troncature
nobs_gbif[nobs_gbif>=150] <- 150
grid_sf$nobs_gbif <- nobs_gbif/st_area(grid_sf) 

# sur une carte
ggplot() +
  geom_sf(data = grid_sf, lwd = 0.1, aes(fill = as.numeric(nobs_gbif))) +
  scale_fill_viridis_c() +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()

# save(grid_sf, file = "grid_sf.RData")

# # ragondins dans l'Occitanie
# ggplot() +
#   geom_sf(data = occitanie, fill = "white", color = "black", lwd = .5) + 
#   geom_sf(data = st_as_sf(data_ragondins, coords = c("decimalLongitude", "decimalLatitude"), crs = st_crs(occitanie)) %>% st_intersection(occitanie), color="red") +
#   geom_sf(data = nutria, color="royalblue") + 
#   theme_minimal() +
#   labs(title = "Ragondins : GBIF (rouge), INPN (bleu) en Occitanie (2019)") +
#   theme_void()

