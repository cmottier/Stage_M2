################################################################################
#              Ragondins 2010 - 2024 : Construction de la grille               #
################################################################################

# librairies utiles ------------------------------------------------------------

library(ggplot2)
library(tidyverse)
library(sf)
# library(tidyterra)
library(terra)

load("Data/A_charger/donnees_a_charger.RData")
load("RData/5km2/grid_sf_5km2.RData")
load("donnees_gbif.RData")

## Densité de population ####################

# ajout de la densité et de log(densité + 1) à grid_sf

# Avec la grille
grid_pop <- pop %>%
  st_intersection(grid_sf) %>%
  group_by(grid_id) %>%
  summarise(hab = sum(TOT_P_2021)) %>%
  as_tibble() %>%
  select(-geometry)

# Densité nulle ...
sum(grid_pop$hab == 0)
summary(grid_pop$hab/grid_sf$area)
min(grid_pop$hab[grid_pop$hab != 0]/grid_sf$area[grid_pop$hab != 0])

# Densité 
grid_sf <- grid_sf %>% 
  full_join(grid_pop, by = "grid_id") %>%
  mutate(density = as.numeric(hab/area)) %>% 
  select(-hab)

grid_sf <- grid_sf %>%
  mutate(logdensity1 = log(density + 1)) # valeur arbitraire pour éviter le 0

rm(grid_pop, pop)


## Observations toutes espèces (GBIF) ####################

periode <- 2010:2024

# recalcul après avoir chargé plus de données

table(gbif_occ$year)
table(gbif_occ$year, gbif_occ$class)

# Nombre d'observations par cellule et par an
nobs_gbif <- NULL
for (annee in periode) {
  obs <- gbif_occ %>%
    st_transform(crs = st_crs(grid_sf)) %>%
    filter(year == annee)
  nobs_gbif[[paste0("gbif_",annee)]] <- lengths(st_intersects(grid_sf, obs))
}
nobs_gbif <- as_tibble(nobs_gbif)

a <- 2019
table(nobs_gbif[[paste0("gbif_",a)]])

# avec une troncature à 100 (voir Twinings)
# Troncature
nobs_gbif[nobs_gbif>=100] <- 100

# Ajout des covariables dans grid_sf (variables par unité d'aire)
grid_sf <- cbind(grid_sf,nobs_gbif) %>%
  mutate(across(starts_with("gbif"), ~.x/area, .names = "d{.col}")) %>%
  select(-starts_with("gbif"))

grid_sf <- grid_sf %>%
  mutate(across(starts_with("dgbif"), ~log(as.numeric(.x) + 10^(-12)), .names = "log_{.col}"))  
# 10^(-12) valeur artificielle à déterminer...

# nombre de cellules avec observation 
for (a in periode) {
  print(length(unique(grid_sf$grid_id[as.numeric(grid_sf[[paste0("dgbif_",a)]])>0])))
  }

# # plot
# ggplot() +
#   geom_sf(data = grid_sf, lwd = 0.1, aes(fill = as.numeric(log_dgbif_2019))) +
#   scale_fill_viridis_c() +
#   geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
#   theme_void()

rm(nobs_gbif)

# Enregistrement de la grille --------------------------------------------------

save(grid_sf, file = "RData/5km2/grid_sf_5km2_bis.RData")

