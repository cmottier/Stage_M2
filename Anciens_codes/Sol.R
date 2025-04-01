
grid_sf <- grid_sf %>%
  mutate(area = st_area(.), .before = 1)

# Données de CORINE révisées en 2006 et 2012...
Cori12 <- st_read("../Data/CLC12_RLRMP_RGF.shp")
str(Cori12)
unique(Cori12$CODE_12)

# Proportion de la surface agricole par maille
agri <- Cori12 %>%
  st_transform(crs = st_crs(grid_sf)) %>%
  filter(str_sub(CODE_12,1,1)=="2")

grid_agri <- agri %>%
  st_crop(occitanie) %>%
  st_intersection(grid_sf) %>%
  mutate(area = st_area(.)) %>%
  group_by(grid_id) %>%
  summarise(aera_agri = sum(area)) %>%
  as_tibble() %>%
  select(-geometry)

grid_agri <- grid_sf %>% 
  inner_join(grid_agri, by = "grid_id") %>%
  mutate(.before = 1, agri_cover = aera_agri/area)

ggplot() +
  geom_sf(data = grid_agri, lwd = 0.1, aes(fill = as.numeric(agri_cover))) +
  scale_fill_viridis_c(
    labels = scales::percent_format()
  ) +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()

# proportion des zones humides 
# faire le point avec Olivier sur les codes à utiliser
humide <- Cori12 %>%
  st_transform(crs = st_crs(grid_sf)) %>%
  filter(CODE_12 %in% c("411", "412", "511", "512"))

grid_hum <- humide %>%
  st_crop(occitanie) %>%
  st_intersection(grid_sf) %>%
  mutate(area = st_area(.)) %>%
  group_by(grid_id) %>%
  summarise(aera_hum = sum(area)) %>%
  as_tibble() %>%
  select(-geometry)

grid_hum <- grid_sf %>% 
  inner_join(grid_hum, by = "grid_id") %>%
  mutate(.before = 1, hum_cover = aera_hum/area)

ggplot() +
  geom_sf(data = grid_hum, lwd = 0.1, aes(fill = as.numeric(hum_cover))) +
  scale_fill_viridis_c(
    labels = scales::percent_format()
  ) +
  geom_sf(data = occitanie, fill = NA, color = "black", lwd = .5) +
  theme_void()
