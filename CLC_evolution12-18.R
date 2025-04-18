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


## Surface agricole ####################

# Données de CORINE révisées 2012 
Cori12 <- st_read("Data/CLC12_RLRMP_RGF.shp")
Cori12 <- Cori12 %>%
  st_transform(crs = st_crs(occitanie)) %>%
  st_intersection(occitanie) 

unique(Cori12$CODE_12) # agricole en 2..

# Surface agricole 
agri12 <- Cori12 %>%
  filter(str_sub(CODE_12,1,1)=="2")


# Données de CORINE révisées 2012 
Cori18 <- st_read("Data/Corine_Land_Cover/CLC2018.gpkg")
Cori18 <- Cori18 %>%
  st_transform(crs = st_crs(occitanie)) %>%
  st_intersection(occitanie) 

str(Cori18)
unique(Cori18$code_18)

# Surface agricole 
agri18 <- Cori18 %>%
  filter(str_sub(code_18,1,1)=="2")

# Comparaison

Diff <- st_sym_difference(agri12, agri18)

ggplot() +
  geom_sf(data = Diff)
