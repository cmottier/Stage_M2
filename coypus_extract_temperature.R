# get temperature through kriging

library(tidyverse)
library(tidyterra)
library(lubridate)
library(sf)
library(KrigR) # https://www.erikkusch.com/courses/krigr/

Dir.Base <- getwd() # identifying the current directory
Dir.Data <- file.path(Dir.Base, "data") # folder path for data

# départements de la région
dpts_occitanie <- st_read("../Data/departements-d-occitanie.shp") 

# contours de la région
loc_site <- dpts_occitanie %>% 
  st_union() %>%
  st_transform(crs = 4326)

st_bbox(loc_site)

occitanie_raster <- st_bbox(c(xmin = -0.35, 
                              xmax = 5, 
                              ymax = 45.1, 
                              ymin = 42), 
                            crs = st_crs(4326)) %>%
  st_as_sfc() %>%
  terra::vect() %>%
  terra::rast()

toccitanie <- CDownloadS(
  Variable = "2m_temperature",
  DataSet = "reanalysis-era5-land-monthly-means",
  Type = "monthly_averaged_reanalysis",
  DateStart = "2003-01-01 00:00",
  DateStop = "2024-12-31 23:00",
  TZone = "Europe/Paris",
  TResolution = "month",
  TStep = 1,
  Extent = occitanie_raster, # our data.frame with Lat and Lon columns
  Dir = Dir.Data,
  FileName = "Toccitanie",
  API_User = "camille.mottier@cefe.cnrs.fr",
  API_Key = "4700337b-89cd-42e8-98b6-a450f702d332"
)

toccitanie <- weathermetrics::kelvin.to.celsius(toccitanie)

# mois le plus froid par an, de 2003 (1) à 2024 (22)
toccitanie_min <- terra::aggregate(toccitanie, fact=c(1,1,12), fun=min)
ggplot() +
  geom_spatraster(data=toccitanie_min$lyr.22) +
  scale_fill_viridis_c() +
  geom_sf(data = dpts_occitanie %>% 
            st_union() %>%
            st_transform(crs = 4326), 
          fill = NA)



# mois le plus chaud par an, de 2003 (1) à 2024 (22)
toccitanie_max <- terra::aggregate(toccitanie, fact=c(1,1,12), fun=max)
ggplot() +
  geom_spatraster(data=toccitanie_max$lyr.22) +
  scale_fill_viridis_c() +
  geom_sf(data = dpts_occitanie %>% 
            st_union() %>%
            st_transform(crs = 4326), 
          fill = NA)

# moyenne annuelle, de 2003 (1) à 2024 (22)
toccitanie_mean <- terra::aggregate(toccitanie, fact=c(1,1,12), fun=mean)
ggplot() +
  geom_spatraster(data=toccitanie_mean$lyr.22) +
  scale_fill_viridis_c() +
  geom_sf(data = dpts_occitanie %>% 
            st_union() %>%
            st_transform(crs = 4326), 
          fill = NA)

# pas compris temp.tif...
terra::writeRaster(x = toccitanie,
                  filename = "temp.tif",
                  overwrite = TRUE)

wtoccitanie <- terra::rast("temp.tif")

# plot d'un certain nombre de mois
Plot.SpatRast(toccitanie$Toccitanie_11) + 
  geom_sf(data = dpts_occitanie %>% 
            st_union() %>%
            st_transform(crs = 4326), 
             fill = NA)

Plot.SpatRast(toccitanie$Toccitanie_12) + 
  geom_sf(data = dpts_occitanie %>% 
            st_union() %>%
            st_transform(crs = 4326), 
          fill = NA)

Plot.SpatRast(toccitanie$Toccitanie_13) + 
  geom_sf(data = dpts_occitanie %>% 
            st_union() %>%
            st_transform(crs = 4326), 
          fill = NA)

Plot.SpatRast(toccitanie$Toccitanie_14) + 
  geom_sf(data = dpts_occitanie %>% 
            st_union() %>%
            st_transform(crs = 4326), 
          fill = NA)



# extract temperature per commune
departements <- st_read("shp/departements-d-occitanie.shp")
mask <- departements$nom_officie.2 == "HERAULT" 
herault <- departements[mask, ]
communes <- st_read("shp/georef-france-commune-millesime.shp")
communes_herault <- communes %>% st_intersection(herault)
noms_communes_piegeage <- c("BAILLARGUES",
                            "CANDILLARGUES",
                            "LA GRANDE-MOTTE",
                            "LANSARGUES",
                            "MARSILLARGUES",
                            "MAUGUIO",
                            "SAINT-NAZAIRE-DE-PÉZAN",
                            "TEYRAN",
                            "ENTRE-VIGNES", #"SAINT-CHRISTOL",
                            "SAINT-AUNÈS",
                            "VALERGUES",
                            "SAINT-JUST",
                            "SAINT-GENIÈS-DES-MOURGUES",
                            "PÉROLS",
                            "LUNEL-VIEL",
                            "LUNEL",
                            "SAINT-VINCENT-DE-BARBEYRARGUES")
communes_piegeage <- communes_herault[communes_herault$com_name_up %in% noms_communes_piegeage,]

temp_communes <- terra::extract(toccitanie, communes_piegeage) %>%
  group_by(ID) %>%
  summarise(tdec2021 = mean(Toccitanie_11, na.rm = T),
            tjan2022 = mean(Toccitanie_12, na.rm = T),
            tfeb2022 = mean(Toccitanie_13, na.rm = T),
            tmar2022 = mean(Toccitanie_14, na.rm = T),
  )



temp_communes <- temp_communes %>% 
  mutate(tdec2021 = weathermetrics::kelvin.to.celsius(tdec2021),
         tjan2022 = weathermetrics::kelvin.to.celsius(tjan2022),
         tfeb2022 = weathermetrics::kelvin.to.celsius(tfeb2022),
         tmar2022 = weathermetrics::kelvin.to.celsius(tmar2022))

saveRDS(temp_communes, "data/temperature.rds")
