################################################################################
#                              Validation                                      #
################################################################################


# Librairies utiles ------------------------------------------------------------
library(ggplot2)
library(tidyverse)
library(sf)
library(spatstat)

# Chargement -------------------------------------------------------------------

load("RData/5km2/grid_sf_5km2.RData")

periode = (2010:2024)

dpts_occitanie <- st_read("Data/departements-d-occitanie.shp") 
occitanie <- dpts_occitanie %>% st_union()
rm(dpts_occitanie)

nutria_periode <- st_read("Data/CEN_2025/Ragondin_rat_musque_pts_2025.shp") %>%
  mutate(year = year(as.Date(DateDebut))) %>% 
  filter(year %in% periode) %>%
  filter(NomVernacu == "Ragondin") %>%
  st_transform(crs = st_crs(occitanie)) %>%
  st_intersection(occitanie)



# Essai quadrat method (pas adapté à notre modèle) -----------------------------
# https://www.paulamoraga.com/book-spatial/complete-spatial-randomness.html

## à la main #################

# On garde les cellules pleines uniquement
grille <- grid_sf %>%
  filter(as.numeric(area) > max(as.numeric(area))-1)

# ggplot() +
#   geom_sf(data = grid_sf, fill = "red") +
#   geom_sf

m <- nrow(grille)

for (annee in 2010:2024) {
  n_star <- sum(grille[[paste0("nnutria", annee)]]) / m
  chi <- sum((grille[[paste0("nnutria", annee)]] - n_star)^2 / n_star)
  print(2*min(pchisq(chi, m - 1, lower.tail = FALSE), pchisq(chi, m - 1)))
}
# On rejette l'hypothèse nulle d'une CSR 

## avec le package spatstat #######################


annee <- 2013

nutria <- nutria_periode %>%
  filter(year == annee)

X <- as.ppp(st_coordinates(nutria), st_bbox(nutria))
Q <- quadratcount(X, nx = 100, ny = 50)
# Q
# plot(X)
# axis(1)
# axis(2)
# plot(Q, add = TRUE, cex = 2)

quadrat.test(Q) # pas significatif : pas CSR
quadrat.test(Q, "regular") # ce test est significatif : clustered pas rejeté
quadrat.test(Q, "clustered")
