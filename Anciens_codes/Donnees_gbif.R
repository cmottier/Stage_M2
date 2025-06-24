################################################################################
#                        Chargement des données gbif                           #
################################################################################

# librairies utiles ------------------------------------------------------------

library(tidyverse)
library(sf)
library(rgbif)

# Occitanie et données de ragondins --------------------------------------------

## Occitanie ###############

# contours des départements d'Occitanie
dpts_occitanie <- st_read("departements-d-occitanie.shp") 

# contours de la région
occitanie <- dpts_occitanie %>% st_union()

rm(dpts_occitanie)


## Ragondins #############

periode = (2010:2024)

# Import des données
nutria_periode <- st_read("Ragondin_rat_musque_pts_2025.shp") %>%
  mutate(year = year(as.Date(DateDebut))) %>% 
  filter(year %in% periode) %>%
  filter(NomVernacu == "Ragondin") 

table(nutria_periode$year)

# changement de format
occitanie <- occitanie %>%
  st_transform(crs = st_crs(nutria_periode))


## Observations toutes espèces (GBIF) ####################

# Récupération des clés taxonomiques des groupes d'intérêt
oiseaux_key <- name_backbone(name = "Aves")$usageKey
mammiferes_key <- name_backbone(name = "Mammalia")$usageKey
amphibiens_key <- name_backbone(name = "Amphibia", rank = "class")$usageKey
reptiles_key <- name_backbone(name = "Reptilia")$usageKey

# Affichage des clés
print(c(oiseaux_key, mammiferes_key, amphibiens_key, reptiles_key))

# Définition de la zone géographique
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
    limit = 2000,  # Maximum d'enregistrements récupérés (à adapter)
    geometry = bbox_wkt  # Zone Occitanie
  )$data
}

# on télécharge 500 données mensuelles pour chaque taxon et chaque année
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
table(gbif_occ$year, gbif_occ$class)

# glimpse(gbif_data)  # Aperçu des colonnes disponibles

# On garde uniquement les colonnes utiles pour l'analyse :
gbif_data <- gbif_data %>%
  select(species, class, order, family, decimalLatitude, decimalLongitude, year, basisOfRecord)

# head(gbif_data)

# Transformer en objet spatial
gbif_data <- st_as_sf(gbif_data, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) 

gbif_data <- gbif_data %>%
  st_transform(crs = st_crs(occitanie))

# En occitanie
gbif_occ <- gbif_data[st_intersects(gbif_data, occitanie, sparse = FALSE), ]

# # Plot
# ggplot() +
#   geom_sf(data = occitanie, fill = "white", color = "black", lwd = .5) + 
#   geom_sf(data = gbif_occ, aes(color = class), alpha = 0.6) +
#   theme_minimal() +
#   facet_wrap(~year, nrow = 4) +
#   labs(color = "Classe") +
#   theme_void()

rm(gbif_data)

save(gbif_occ, file = "donnees_gbif.RData")
