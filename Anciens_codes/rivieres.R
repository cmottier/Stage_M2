# Rivières 
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

# bind river
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

# format
plan_eau <- plan_eau %>%
  st_transform(crs = st_crs(grid_sf))

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

# distance aux plans d'eau

dist_eau <- st_distance(grid_sf, st_transform(plan_eau, crs = st_crs(grid_sf)))
# save(dist_eau, file = "Data/dist_eau.RData")

# On regarde la distance minimale
dist_eau_min <- apply(dist_eau, 1, min)

ggplot() +
  geom_sf(data = grid_sf, lwd = 0.1, aes(fill = dist_eau_min)) +
  scale_fill_viridis_c(
  ) +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()

# Pour essayer d'être plus efficace...
dist = 15000 # rayon de recherche autour de la cellule
f_buff <- function (i) {
  cell = grid_sf[i,]
  buffer <- st_buffer(cell, dist) 
  plan_filtre <- plan_eau[st_intersects(plan_eau, buffer, sparse = F),]
  return(min(st_distance(cell, plan_filtre)))
}

for (i in 1:nrow(grid_sf)) {
  grid_sf$dist_plan_eau[i] <- f_buff(i)
}
# apply ne marche pas...

ggplot() +
  geom_sf(data = grid_sf, lwd = 0.1, aes(fill = as.numeric(dist_plan_eau))) +
  scale_fill_viridis_c(
  ) +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()

ggplot() +
  geom_sf(data = grid_sf[14,]) +
  geom_sf(data = st_buffer(grid_sf[14,], dist), fill = NA) +
  geom_sf(data=plan_eau_occ[st_intersects(plan_eau_occ, st_buffer(grid_sf[14,], dist), sparse = F),])

min(st_distance(grid_sf[14,], plan_eau_occ[st_intersects(plan_eau_occ, st_buffer(grid_sf[14,], dist), sparse = F),]))
