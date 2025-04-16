################################################################################
#              Ragondins 2010 - 2024 : Construction de la grille               #
################################################################################

# librairies utiles ------------------------------------------------------------

library(ggplot2)
library(tidyverse)
library(sf)
library(tidyterra)
library(corrplot)
library(terra)


# Chargement des données (indépendantes de la grille) --------------------------

load("Data/A_charger/donnees_a_charger.RData")
toccitanie_min <- rast("Data/A_charger/toccitanie_min.tif")
toccitanie_max <- rast("Data/A_charger/toccitanie_max.tif")
toccitanie_mean <- rast("Data/A_charger/toccitanie_mean.tif")
poccitanie_cum <- rast("Data/A_charger/poccitanie_cum.tif")



# Grille -----------------------------------------------------------------------

## Construction de la grille ###########

# Choix de la surface des cellules
grid_area <- units::set_units(50000000,"m^2") #50km2
grid <- st_make_grid(occitanie, cellsize = grid_area, what = "polygons", square = FALSE)

# sf
grid_sf <- st_sf(grid)

rm(grid, grid_area)

ggplot() +
  geom_sf(data = grid_sf, lwd= 0.1) + 
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) + 
  theme_void()

# Cellules d'Occitanie uniquement
grid_sf <- grid_sf[st_intersects(grid_sf, 
                                 occitanie, 
                                 sparse = FALSE), ] %>%
  st_intersection(occitanie) # on rogne les cellules de la frontière

ggplot() +
  geom_sf(data = grid_sf$grid, lwd = 0.1) +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) + 
  theme_void()

# Surface exacte des cellules en occitanie
grid_sf$area <- st_area(grid_sf)

# On supprime les cellules de surface nulle
grid_sf <- grid_sf %>%
  filter(as.numeric(area)!=0)

# Ajout de l'identifiant des cellules
grid_sf <- grid_sf %>%
  mutate(.before = 1, grid_id = 1:lengths(grid_sf)[1])

# Centroides des cellules
centroide <- st_centroid(grid_sf$grid)

# save(grid_sf, file = "RData/50km2/grid_sf_50km2.RData")


## Ragondins par cellule et par an ####################

# Nombre d'observations de ragondins par cellules
for (annee in periode) {
  nutria_annee <- nutria_periode %>%
    filter(year == annee)
  grid_sf[,paste0("nnutria", annee)] <- lengths(st_intersects(grid_sf, nutria_annee))
}

rm(nutria_annee, nutria_periode)

# # Résumé
# nnutria_names <- paste0(rep("nnutria", length(periode)), periode)
# 
# for (name in nnutria_names) {
#   print(name)
#   print(table(grid_sf[[name]]))
#   }
# 
# for (name in nnutria_names) {
#   print(name)
#   print(which.max(grid_sf[[name]]))
# }
# 
# rm(name)

# save(grid_sf, file = "RData/50km2/grid_sf_50km2.RData")




# Variables explicatives -------------------------------------------------------

## Températures ####################

# avec la grille 
temp_min <- terra::extract(toccitanie_min, grid_sf) %>%
  group_by(ID) %>%
  summarise(across(starts_with("lyr"), mean)) %>%
  select(-ID) %>%
  setNames(paste0("tmin_", periode))

temp_max <- terra::extract(toccitanie_max, grid_sf) %>%
  group_by(ID) %>%
  summarise(across(starts_with("lyr"), mean)) %>%
  select(-ID) %>%
  setNames(paste0("tmax_", periode))

temp_mean <- terra::extract(toccitanie_mean, grid_sf) %>%
  group_by(ID) %>%
  summarise(across(starts_with("lyr"), mean)) %>%
  select(-ID) %>%
  setNames(paste0("tmean_", periode))

# Corrélations
for (annee in periode) {
  print(cor(cbind(temp_min[[paste0("tmin_", annee)]], 
            temp_max[[paste0("tmax_", annee)]], 
            temp_mean[[paste0("tmean_", annee)]])))
}
# on garde uniquement tmin...

# ajout des covariables dans grif_sf
grid_sf <- cbind(grid_sf, temp_min)

# plot
ggplot() +
  geom_sf(data = grid_sf, lwd = 0.1, color = NA, aes(fill = tmin_2019)) +
  scale_fill_viridis_c() +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()

# save(grid_sf, file = "RData/50km2/grid_sf_50km2.RData")

rm(annee, temp_min, temp_max, temp_mean, toccitanie_min, toccitanie_max, toccitanie_mean)



## Précipitations ####################

# avec la grille 
prec_cum <- terra::extract(poccitanie_cum, grid_sf) %>%
  group_by(ID) %>%
  summarise(across(starts_with("lyr"), mean)) %>%
  select(-ID) %>%
  setNames(paste0("pcum_", periode))

# ajout des covariables dans grif_sf
grid_sf <- cbind(grid_sf, prec_cum)

# plot
ggplot() +
  geom_sf(data = grid_sf, lwd = 0.1, color = NA, aes(fill = pcum_2019)) +
  scale_fill_viridis_c() +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()

# save(grid_sf, file = "RData/50km2/grid_sf_50km2.RData")

rm(poccitanie_cum, prec_cum)



## Surface agricole ####################

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
  mutate(agri_cover = aera_agri/area) %>%
  select(-aera_agri)

# On remplace les NA par 0 (pas de jointure = pas de terres agricoles)
grid_sf$agri_cover[is.na(grid_sf$agri_cover)] <- 0

# plot
ggplot() +
  geom_sf(data = grid_sf, lwd = 0.1, color = NA, aes(fill = as.numeric(agri_cover))) +
  scale_fill_viridis_c(
    labels = scales::percent_format()
  ) +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()

# save(grid_sf, file = "RData/50km2/grid_sf_50km2.RData")

rm(agri, grid_agri)
   
   
## Densité de population ####################

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
  mutate(logdensity = log(as.numeric(hab/area) + 10^(-20))) %>% # valeur arbitraire pour éviter le 0
  select(-hab)

# Plot
ggplot() +
  geom_sf(data = grid_sf, lwd = 0.1, aes(fill = logdensity)) +
  scale_fill_viridis_c() +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()

# save(grid_sf, file = "RData/50km2/grid_sf_50km2.RData")

rm(grid_pop, pop)


## Distances aux chemins ####################

# On cherche l'indice du chemin le plus proche
index <- st_nearest_feature(centroide, chemins)

# On garde la distance associée
grid_sf$dist_chemins <- st_distance(centroide, chemins[index,], by_element = TRUE)

# plot
ggplot() +
  geom_sf(data = grid_sf, color = NA, aes(fill = as.numeric(dist_chemins))) +
  scale_fill_viridis_c(
  ) +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()

# save(grid_sf, file = "RData/50km2/grid_sf_50km2.RData")

rm(chemins, index)


## Distances aux routes ####################

# On cherche l'indice de la route la plus proche
index <- st_nearest_feature(centroide, routes)

# On garde la distance associée
grid_sf$dist_routes <- st_distance(centroide, routes[index,], by_element = TRUE)

# plot
ggplot() +
  geom_sf(data = grid_sf, color = NA, aes(fill = as.numeric(dist_routes))) +
  scale_fill_viridis_c(
  ) +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()

# save(grid_sf, file = "RData/50km2/grid_sf_50km2.RData")

rm(routes)


## Routes et chemins ensemble #################

grid_sf$dist_acces <- pmin(grid_sf$dist_chemins, grid_sf$dist_routes)

# plot
ggplot() +
  geom_sf(data = grid_sf, color = NA, aes(fill = as.numeric(dist_acces))) +
  scale_fill_viridis_c(
  ) +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()

# save(grid_sf, file = "RData/50km2/grid_sf_50km2.RData")


## Distances aux cours d'eau ####################

# On cherche l'indice de la rivière la plus proche
index <- st_nearest_feature(centroide, rivieres)

# On garde la distance associée
grid_sf$dist_rivieres <- st_distance(centroide, rivieres[index,], by_element = TRUE)

# plot 
ggplot() +
  geom_sf(data = grid_sf[as.numeric(grid_sf$dist_rivieres)!=0,], color = NA, aes(fill = as.numeric(dist_rivieres))) +
  scale_fill_viridis_c(
  ) +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()

# save(grid_sf, file = "RData/50km2/grid_sf_50km2.RData")

rm(rivieres, index)


## Distance aux plans d'eau ####################

# On cherche l'indice du plan d'eau le plus proche
index <- st_nearest_feature(centroide, plan_eau_tot)

# On garde la distance associée
grid_sf$dist_plan_eau <- st_distance(centroide, plan_eau_tot[index,], by_element = TRUE)

# plot
ggplot() +
  geom_sf(data = grid_sf, color = NA, aes(fill = as.numeric(dist_plan_eau))) +
  scale_fill_viridis_c(
  ) +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()

# save(grid_sf, file = "RData/50km2/grid_sf_50km2.RData")

rm(plan_eau_tot, index)

## Rivières et plans d'eau ensemble ######################

grid_sf$dist_eau <- pmin(grid_sf$dist_rivieres, grid_sf$dist_plan_eau)

# plot
ggplot() +
  geom_sf(data = grid_sf, color = NA, aes(fill = as.numeric(dist_eau))) +
  scale_fill_viridis_c(
  ) +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()

# save(grid_sf, file = "RData/50km2/grid_sf_50km2.RData")



## Observations toutes espèces (GBIF) ####################

# Nombre d'observations par cellule et par an
nobs_gbif <- NULL
for (annee in periode) {
  obs <- gbif_occ %>%
    filter(year == annee)
  nobs_gbif[[paste0("gbif_",annee)]] <- lengths(st_intersects(grid_sf, obs))
}
nobs_gbif <- as_tibble(nobs_gbif)

table(nobs_gbif$gbif_2010)

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

# plot
ggplot() +
  geom_sf(data = grid_sf, lwd = 0.1, aes(fill = as.numeric(log_dgbif_2019))) +
  scale_fill_viridis_c() +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()


# Enregistrement de la grille --------------------------------------------------

save(grid_sf, file = "RData/50km2/grid_sf_50km2.RData")

