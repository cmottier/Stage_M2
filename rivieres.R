# Rivières 
Ariege <- sf::st_read("../Data/Rivieres/Ariege/COURS_D_EAU.shp")
Aude <- sf::st_read("../Data/Rivieres/Aude/COURS_D_EAU.shp")
Aveyron <- sf::st_read("../Data/Rivieres/Aveyron/COURS_D_EAU.shp")
Gard <- sf::st_read("../Data/Rivieres/Gard/COURS_D_EAU.shp")
HauteGaronne <- sf::st_read("../Data/Rivieres/HauteGaronne/COURS_D_EAU.shp")
Gers <- sf::st_read("../Data/Rivieres/Gers/COURS_D_EAU.shp")
Herault <- sf::st_read("../Data/Rivieres/Herault/COURS_D_EAU.shp")
Lot <- sf::st_read("../Data/Rivieres/Lot/COURS_D_EAU.shp")
Lozere <- sf::st_read("../Data/Rivieres/Lozere/COURS_D_EAU.shp")
HautesPyrenees <- sf::st_read("../Data/Rivieres/HautesPyrenees/COURS_D_EAU.shp")
PyreneesOrientales <- sf::st_read("../Data/Rivieres/PO/COURS_D_EAU.shp")
Tarn <- sf::st_read("../Data/Rivieres/Tarn/COURS_D_EAU.shp")
TarnetGaronne <- sf::st_read("../Data/Rivieres/TarnEtGaronne/COURS_D_EAU.shp")

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

# Plan d'eau
Ariege <- sf::st_read("../Data/Rivieres/Ariege/PLAN_D_EAU.shp")
Aude <- sf::st_read("../Data/Rivieres/Aude/PLAN_D_EAU.shp")
Aveyron <- sf::st_read("../Data/Rivieres/Aveyron/PLAN_D_EAU.shp")
Gard <- sf::st_read("../Data/Rivieres/Gard/PLAN_D_EAU.shp")
HauteGaronne <- sf::st_read("../Data/Rivieres/HauteGaronne/PLAN_D_EAU.shp")
Gers <- sf::st_read("../Data/Rivieres/Gers/PLAN_D_EAU.shp")
Herault <- sf::st_read("../Data/Rivieres/Herault/PLAN_D_EAU.shp")
Lot <- sf::st_read("../Data/Rivieres/Lot/PLAN_D_EAU.shp")
Lozere <- sf::st_read("../Data/Rivieres/Lozere/PLAN_D_EAU.shp")
HautesPyrenees <- sf::st_read("../Data/Rivieres/HautesPyrenees/PLAN_D_EAU.shp")
PyreneesOrientales <- sf::st_read("../Data/Rivieres/PO/PLAN_D_EAU.shp")
Tarn <- sf::st_read("../Data/Rivieres/Tarn/PLAN_D_EAU.shp")
TarnetGaronne <- sf::st_read("../Data/Rivieres/TarnEtGaronne/PLAN_D_EAU.shp")

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

# avec la grille
plan_eau_occ <- plan_eau %>%
  st_transform(crs = st_crs(grid_sf)) %>%
  st_intersection(grid_sf) %>%
  mutate(area = st_area(.)) %>%
  group_by(grid_id) %>%
  summarise(aera_eau = sum(area)) %>%
  as_tibble() %>%
  select(-geometry)

grid_sf <- grid_sf %>% 
  full_join(plan_eau_occ, by = "grid_id") %>%
  mutate(.before = 1, surface_en_eau = aera_eau/area) %>%
  select(-aera_eau)

ggplot() +
  geom_sf(data = grid_sf, lwd = 0.1, aes(fill = as.numeric(surface_en_eau))) +
  scale_fill_viridis_c(
    labels = scales::percent_format()
  ) +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()




