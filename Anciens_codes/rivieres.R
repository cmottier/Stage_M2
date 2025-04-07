## Longueur de cours d'eau ####################

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

ggplot() +
  geom_sf(data = river_lines, lwd = 0.1) +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()

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

# save(grid_sf, file = "RData/grid_sf_5km2.RData")

# Essai des distances à la rivière la plus proche

index <- st_nearest_feature(grid_sf[8578:nrow(grid_sf),], river_lines)
grid_sf$dist_rivieres[8578:nrow(grid_sf)] <- st_distance(grid_sf[8578:nrow(grid_sf),], element[index,])


distance_min <- function(cell, grid, element) {
  index <- st_nearest_feature(grid[cell,], element)
  retrun(st_distance(grid[cell,], element[index,]))
}

# on applique aux rivières
for (i in 7600:nrow(grid_sf)) {
  grid_sf$dist_rivieres[i] <- distance_min(i, grid_sf, river_lines)
}

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

##################
# Espagne


rivieres <- st_read("Data/Rivieres/Espagne/HydroRIVERS_v10_eu.shp")

# On garde uniquement le buffer autour de la frontiere
rivieres_espagne <- rivieres %>%
  st_transform(crs = st_crs(grid_sf)) %>%
  st_crop(box) %>% # on garde que la box
  st_intersection(espagne_buff) # Espagne uniquement

ggplot() +
  geom_sf(data = rivieres_espagne, color = "royalblue") +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5)
