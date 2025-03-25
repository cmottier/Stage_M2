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
get_gbif_data <- function(taxon_key) {
  occ_search(
    # taxonKey = taxon_key,
    classKey = taxon_key, # proposition de remplacement ? 
    year = 2019,
    hasCoordinate = TRUE,  # Seulement les données géolocalisées
    limit = 5000,  # !!!!!!!!!!!!!!! Maximum d'enregistrements récupérés (à adapter)
    geometry = paste(bbox_wkt, collapse = ",")  # Zone Occitanie
  )$data
}
# problème : données classées par mois -> on conserve ici des données de janvier à avril uniquement...

# Télécharger les données pour chaque groupe
data_oiseaux <- get_gbif_data(oiseaux_key)
data_mammiferes <- get_gbif_data(mammiferes_key)
data_amphibiens <- get_gbif_data(amphibiens_key)
data_reptiles <- get_gbif_data(reptiles_key)

# ragondins ?
data_ragondins <- occ_search(
    taxonKey = 4264680,
    year = 2019,
    hasCoordinate = TRUE,  # Seulement les données géolocalisées
    limit = 5000,  # !!!!!!!!!!!!!!! Maximum d'enregistrements récupérés (à adapter)
    geometry = paste(bbox_wkt, collapse = ",")  # Zone Occitanie
  )$data




# Fusionner les jeux de données
gbif_data <- bind_rows(data_oiseaux, data_mammiferes, data_amphibiens, data_reptiles)


dim(gbif_data)  # Nombre de lignes et colonnes
glimpse(gbif_data)  # Aperçu des colonnes disponibles

# On garde uniquement les colonnes utiles pour l'analyse :
gbif_clean <- gbif_data %>%
  select(species, class, order, family, decimalLatitude, decimalLongitude, year, basisOfRecord)
head(gbif_clean)

# Transformer en objet spatial
gbif_sf <- st_as_sf(gbif_clean, coords = c("decimalLongitude", "decimalLatitude"), crs = st_crs(occitanie))

# Afficher sur une carte
ggplot() +
  geom_sf(data = occitanie, fill="white", color = "black", lwd = .5) + 
  geom_sf(data = gbif_sf, aes(color = class), alpha = 0.6) +
  theme_minimal() +
  labs(title = "Observations GBIF en Occitanie (2019)", color = "Classe") +
  theme_void()

# classes étranges...
table(gbif_clean$class)

# dans l'Occitanie
ggplot() +
  geom_sf(data = occitanie, fill = "white", color = "black", lwd = .5) + 
  geom_sf(data = gbif_sf %>% st_intersection(occitanie), aes(color = class), alpha = 0.6) +
  theme_minimal() +
  labs(title = "Observations GBIF en Occitanie (2019)", color = "Classe") +
  theme_void()

# ragondins dans l'Occitanie
ggplot() +
  geom_sf(data = occitanie, fill = "white", color = "black", lwd = .5) + 
  geom_sf(data = st_as_sf(data_ragondins, coords = c("decimalLongitude", "decimalLatitude"), crs = st_crs(occitanie)) %>% st_intersection(occitanie), color="red") +
  geom_sf(data = nutria, color="royalblue") + 
  theme_minimal() +
  labs(title = "Ragondins : GBIF (rouge), INPN (bleu) en Occitanie (2019)") +
  theme_void()

