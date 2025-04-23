################################################################################
#              Evolution de la couverture agricole 2012-2018                   #
################################################################################

# librairies utiles ------------------------------------------------------------

library(ggplot2)
library(tidyverse)
library(sf)



## Occitanie ###############

# contours des départements d'Occitanie
dpts_occitanie <- st_read("Data/departements-d-occitanie.shp") 

# contours de la région
occitanie <- dpts_occitanie %>% st_union()

rm(dpts_occitanie)

# Box
occitanie_box <- st_bbox(c(
  xmin = -0.4,
  xmax = 5,
  ymax = 45.1,
  ymin = 42
), crs = st_crs(4326)) %>%
  st_as_sfc() %>%
  st_transform(crs = st_crs(occitanie))


## Surface agricole ####################

# Données de CORINE révisées 2012 
Cori12 <- st_read("Data/Corine_Land_Cover/CLC12_RLRMP_RGF.shp")

unique(Cori12$CODE_12) # agricole en 2..

# Surface agricole 
agri12 <- Cori12 %>%
  filter(str_sub(CODE_12,1,1)=="2")

agri12 <- agri12 %>%
  st_transform(crs = st_crs(occitanie)) %>%
  st_crop(occitanie_box) %>%
  st_intersection(occitanie) 

# save(agri12, file = "Data/Corine_Land_Cover/agri12.RData")

# Données de CORINE révisées 2018
Cori18 <- st_read("Data/Corine_Land_Cover/CLC2018.gpkg")

str(Cori18)
unique(Cori18$code_18)

# Surface agricole 
agri18 <- Cori18 %>%
  filter(str_sub(code_18,1,1)=="2") 

agri18 <- agri18 %>%
  select(geom)

rm(Cori18)

agri18 <- agri18 %>%
  st_transform(crs = st_crs(occitanie)) %>%
  st_crop(occitanie_box) %>%
  st_intersection(occitanie) 

# save(agri18, file = "Data/Corine_Land_Cover/agri18.RData")


# Comparaison
# Diff <- st_sym_difference(agri12, agri18)

changements <- st_read("Data/Corine_Land_Cover/U2018_CHA1218_V2020_20u1.shp")
changements <- changements %>%
  st_transform(crs = st_crs(occitanie)) %>%
  st_intersection(occitanie) %>%
  filter(xor(str_sub(Code_18,1,1)=="2", str_sub(Code_12,1,1)=="2" ))

augmentations <- changements %>%
  filter(str_sub(Code_18,1,1)=="2" & str_sub(Code_12,1,1)!="2")

diminutions <- changements %>%
  filter(str_sub(Code_18,1,1)!="2" & str_sub(Code_12,1,1)=="2")

ggplot() +
  geom_sf(data = occitanie) +
  geom_sf(data = augmentations, fill = "royalblue") +
  geom_sf(data = diminutions, fill = "coral")

