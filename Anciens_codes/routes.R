# Calcul des longueurs par cellules -> idée abandonnée

###########
# Chemins
# OpenStreetMap https://www.data.gouv.fr/fr/datasets/itineraires-de-randonnee-dans-openstreetmap/
chemins_osm <- st_read("Data/hiking_foot_routes_lineLine.shp")
chemins_osm <- chemins_osm %>% 
  st_transform(crs = st_crs(grid_sf)) %>%
  st_intersection(occitanie_buff)

ggplot() +
  geom_sf(data = chemins_osm, color = "royalblue") +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5)

# Avec la grille
grid_chemin <- chemins_osm  %>%
  st_intersection(grid_sf) %>%
  mutate(len = st_length(.)) %>%
  group_by(grid_id) %>%
  summarise(lgr = sum(len)) %>%
  as_tibble() %>%
  select(-geometry)

# longueur/surface
grid_sf <- grid_sf %>%
  full_join(grid_chemin, by = "grid_id") %>%
  mutate(.before = 1, lgr_chemins = lgr/area) %>%
  select(-lgr)


ggplot() +
  geom_sf(data = grid_sf, lwd = 0.1, aes(fill = as.numeric(lgr_chemins))) +
  labs(fill = "Longueur de chemins (km/m^2)") +
  scale_fill_viridis_c() +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()


###########
# Routes
roads <- st_read("Data/TRONCON_ROUTE.shp")
str(roads)
unique(roads$CLASS_ADM)

routes_occ <- roads %>%
  st_transform(crs = st_crs(grid_sf)) %>%
  st_intersection(occitanie) %>%
  filter(CLASS_ADM != "Autoroute" & CLASS_ADM !="Sans objet") # D\xe9partementale ne marche pas...

# Avec la grille
grid_route <- routes_occ  %>%
  st_intersection(grid_sf) %>%
  mutate(len = st_length(.)) %>%
  group_by(grid_id) %>%
  summarise(lgr = sum(len)) %>%
  as_tibble() %>%
  select(-geometry)

# longueur/surface 
grid_sf <- grid_sf %>% 
  full_join(grid_route, by = "grid_id") %>%
  mutate(.before = 1, lgr_routes = lgr/area) %>%
  select(-lgr)

ggplot() +
  geom_sf(data = grid_sf, lwd = 0.1, aes(fill = as.numeric(lgr_routes))) +
  labs(fill = "Longueur de routes (km/m^2)") +
  scale_fill_viridis_c() +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()
