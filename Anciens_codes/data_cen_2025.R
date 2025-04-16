
library(ggplot2)
library(tidyverse)
library(sf)

# Occitanie et données de ragondins --------------------------------------------

## Occitanie ###############

# contours des départements d'Occitanie
dpts_occitanie <- st_read("Data/departements-d-occitanie.shp") 

# contours de la région
occitanie <- dpts_occitanie %>% st_union()

rm(dpts_occitanie)


## Ragondins 2019 #############

# Import des données
data_CEN_lin <- st_read("Data/CEN_2025/Ragondin_rat_musque_lin_2025.shp") %>%
  mutate(year = year(as.Date(DateDebut))) 
data_CEN_poly <- st_read("Data/CEN_2025/Ragondin_rat_musque_poly_2025.shp")  %>%
  mutate(year = year(as.Date(DateDebut))) 
data_CEN_pts <- st_read("Data/CEN_2025/Ragondin_rat_musque_pts_2025.shp")  %>%
  mutate(year = year(as.Date(DateDebut))) 

names(data_CEN_lin)
names(data_CEN_poly)
names(data_CEN_pts)

table(data_CEN_lin$year)
table(data_CEN_poly$year)
table(data_CEN_pts$year)

# changement de format
occitanie <- occitanie %>%
  st_transform(crs = st_crs(data_CEN_lin))

# plot des données 
p <- ggplot() +
  geom_sf(data = occitanie, fill = "white", color = "black", lwd = .5) + 
  # geom_sf(data = data_CEN_lin, color = "blue") +
  geom_sf(data = data_CEN_poly, color = "red") +
  geom_sf(data = data_CEN_pts, color = "chartreuse4") +
  labs(color = "") +
  theme_void()
p

# intersection avec occitanie
nutria <- st_intersection(nutria, occitanie)

# plot des données 
p <- ggplot() +
  geom_sf(data = occitanie, fill = "white", color = "black", lwd = .5) + 
  # geom_sf(data = occitanie_buff, fill = NA, color = "red", lwd = .5) + 
  geom_sf(data = nutria) + 
  labs(color = "") +
  theme_void()
p