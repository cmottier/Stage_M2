# # Chemins de randonnées https://geoservices.ign.fr/sites/default/files/2023-10/DC_BDTOPO_Transport.pdf
# chemins_ariege <- st_read("../Data/Rando/Ariege/ITINERAIRE_AUTRE.shp")
# chemins_aude <- st_read("../Data/Rando/Aude/ITINERAIRE_AUTRE.shp")
# chemins_aveyron <- st_read("../Data/Rando/Aveyron/ITINERAIRE_AUTRE.shp")
# chemins_gard <- st_read("../Data/Rando/Gard/ITINERAIRE_AUTRE.shp")
# chemins_gers <- st_read("../Data/Rando/Gers/ITINERAIRE_AUTRE.shp")
# chemins_hautegaronne <- st_read("../Data/Rando/HauteGaronne/ITINERAIRE_AUTRE.shp")
# chemins_hautespyrenees <- st_read("../Data/Rando/HautesPyrenees/ITINERAIRE_AUTRE.shp")
# chemins_herault <- st_read("../Data/Rando/Herault/ITINERAIRE_AUTRE.shp")
# chemins_lot <- st_read("../Data/Rando/Lot/ITINERAIRE_AUTRE.shp")
# chemins_lozere <- st_read("../Data/Rando/Lozere/ITINERAIRE_AUTRE.shp")
# chemins_po <- st_read("../Data/Rando/PO/ITINERAIRE_AUTRE.shp")
# chemins_tarn <- st_read("../Data/Rando/Tarn/ITINERAIRE_AUTRE.shp")
# chemins_tarnetgaronne <- st_read("../Data/Rando/TarnEtGaronne/ITINERAIRE_AUTRE.shp")
# 
# # on rassemble
# chemins <- chemins_ariege %>%
#   rbind(chemins_aude) %>%
#   rbind(chemins_aveyron) %>%
#   rbind(chemins_gard) %>%
#   rbind(chemins_gers) %>%
#   rbind(chemins_hautegaronne) %>%
#   rbind(chemins_hautespyrenees) %>%
#   rbind(chemins_herault) %>%
#   rbind(chemins_lot) %>%
#   rbind(chemins_lozere) %>%
#   rbind(chemins_po) %>%
#   rbind(chemins_tarn) %>%
#   rbind(chemins_tarnetgaronne) %>%
#   sf::st_transform(crs = st_crs(grid_sf)) %>%
#   sf::st_simplify()
# 
# chemins_occ <- chemins %>%
#   st_intersection(grid_sf)
# 
# ggplot() +
#   geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
#   geom_sf(data = chemins_occ, color = "grey") +
#   theme_void()
# 
# # calcule la longueur des segments
# chemins_occ$length <- st_length(chemins_occ)
# 
# # somme des segments par cellule
# longueur_par_cellule <- chemins_occ %>%
#   group_by(grid_id) %>% 
#   summarise(lgr_totale = sum(length))  
# 
# # jointure
# grid_sf$lgr_chemins <- 0
# for (i in 1:nrow(grid_sf)){
#   mask <- grid_sf$grid_id[i]
#   if (sum(longueur_par_cellule$grid_id == mask) == 0) next
#   grid_sf$lgr_chemins[i] <- longueur_par_cellule$lgr_totale[longueur_par_cellule$grid_id == mask]
# }
# 
# ggplot() +
#   geom_sf(data = grid_sf, lwd = 0.1, aes(fill = lgr_chemins/1000)) +
#   labs(fill = "Longueur de chemins (km)") +
#   scale_fill_viridis_c() +
#   geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
#   theme_void()

# Chemins (OpenStreetMap) https://www.data.gouv.fr/fr/datasets/itineraires-de-randonnee-dans-openstreetmap/
chemins_osm <- st_read("../Data/hiking_foot_routes_lineLine.shp")
chemins_osm <- chemins_osm %>% 
  st_transform(crs = st_crs(grid_sf)) 

# une même ligne peut apparaître plusieurs fois ? 

chemins_occ <- chemins_osm %>%
  st_intersection(grid_sf)

ggplot() +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  geom_sf(data = chemins_occ, color = "grey") +
  theme_void()

# calcule la longueur des segments
chemins_occ$length <- st_length(chemins_occ)

# somme des segments par cellule
longueur_par_cellule <- chemins_occ %>%
  group_by(grid_id) %>% 
  summarise(lgr_totale = sum(length))  

# jointure
grid_sf$lgr_chemins <- 0
for (i in 1:nrow(grid_sf)){
  mask <- grid_sf$grid_id[i]
  if (sum(longueur_par_cellule$grid_id == mask) == 0) next
  grid_sf$lgr_chemins[i] <- longueur_par_cellule$lgr_totale[longueur_par_cellule$grid_id == mask]
}

ggplot() +
  geom_sf(data = grid_sf, lwd = 0.1, aes(fill = lgr_chemins/1000)) +
  labs(fill = "Longueur de chemins (km)") +
  scale_fill_viridis_c() +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()


# Routes
roads <- st_read("../Data/TRONCON_ROUTE.shp")
str(roads)
unique(roads$CLASS_ADM)

routes_occ <- roads %>%
  st_transform(crs = st_crs(grid_sf)) %>%
  st_intersection(grid_sf) %>%
  filter(CLASS_ADM %in% c("D\xe9partementale","Nationale"))

ggplot() +
  geom_sf(data = grid_sf, fill = "white", color = "black", lwd = .5) +
  geom_sf(data = routes_occ, color = "gray") +
  theme_void()

# calcule la longueur des segments
routes_occ$length <- st_length(routes_occ)

# somme des segments par cellule
longueur_par_cellule <- routes_occ %>%
  group_by(grid_id) %>%
  summarise(lgr_totale = sum(length))  # Somme des longueurs

# jointure
grid_sf$lgr_routes <- 0
for (i in 1:nrow(grid_sf)){
  mask <- grid_sf$grid_id[i]
  if (sum(longueur_par_cellule$grid_id == mask) == 0) next
  grid_sf$lgr_routes[i] <- longueur_par_cellule$lgr_totale[longueur_par_cellule$grid_id == mask]
}

ggplot() +
  geom_sf(data = grid_sf, lwd = 0.1, aes(fill = lgr_routes/1000)) +
  labs(fill = "Longueur de routes (km)") +
  scale_fill_viridis_c() +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()
