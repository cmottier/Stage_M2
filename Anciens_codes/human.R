library(sf)

pop <- st_read("../Data/pop2021.gpkg")
pop <- pop %>% 
  st_transform(crs = st_crs(grid_sf)) %>%
  st_crop(occitanie)

ggplot() +
  geom_sf(data = pop, lwd = 0.1, aes(fill = log(TOT_P_2021))) +
  scale_fill_viridis_c() +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()

ggplot() +
  geom_sf(data = grid_sf[grid_sf$grid_id %in% c(70, 94, 123, 196, 235),], lwd = 0.1, fill = "black") +
  scale_fill_viridis_c() +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()


grid_pop <- pop %>%
  st_intersection(grid_sf) %>%
  group_by(grid_id) %>%
  summarise(hab = sum(TOT_P_2021)) %>%
  as_tibble() %>%
  select(-geom)

grid_pop <- grid_sf %>% 
  full_join(grid_pop, by = "grid_id") %>%
  mutate(.before = 1, logdensity = log(hab/area))

# problème des cellules -inf et NA (cellules manquantes sur la frontière)... 

ggplot() +
  geom_sf(data = grid_pop, lwd = 0.1, aes(fill = as.numeric(logdensity))) +
  scale_fill_viridis_c() +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()
