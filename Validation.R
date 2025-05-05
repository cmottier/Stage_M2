################################################################################
#                              Validation                                      #
################################################################################


# Librairies utiles ------------------------------------------------------------
library(ggplot2)
library(tidyverse)
library(sf)
library(spatstat)

# Chargement des données utiles ------------------------------------------------

# Occitanie
dpts_occitanie <- st_read("Data/departements-d-occitanie.shp") 
occitanie <- dpts_occitanie %>% st_union()
rm(dpts_occitanie)

# observations
periode = (2010:2024)

nutria_periode <- st_read("Data/CEN_2025/Ragondin_rat_musque_pts_2025.shp") %>%
  mutate(year = year(as.Date(DateDebut))) %>% 
  filter(year %in% periode) %>%
  filter(NomVernacu == "Ragondin") 

crs_commun <- st_crs(nutria_periode)

occitanie <- st_transform(occitanie, crs = crs_commun)
  
nutria_periode <- nutria_periode %>%
  st_intersection(occitanie)

# grille à utiliser
load("RData/5km2/grid_sf_5km2.RData")

# année étudiée et coefficients estimés associés
annee <- 2016

load("Resultats_MCMC/5km2_avec_lasso/IEP_2016.RData")

IEP_2016 <- IEP_2016 %>%
  st_transform(crs = crs_commun)

estim <- IEP_2016 %>%
  select(grid, lambda_med, b_med)




# Quadrat method (met en évidence le clustering ?) -----------------------------

# https://www.paulamoraga.com/book-spatial/complete-spatial-randomness.html

## à la main #################

# On garde les cellules pleines uniquement
grille <- grid_sf %>%
  filter(as.numeric(area) > max(as.numeric(area))-1)

# ggplot() +
#   geom_sf(data = grid_sf, fill = "red") +
#   geom_sf

m <- nrow(grille)

for (annee in 2010:2024) {
  n_star <- sum(grille[[paste0("nnutria", annee)]]) / m
  chi <- sum((grille[[paste0("nnutria", annee)]] - n_star)^2 / n_star)
  print(2*min(pchisq(chi, m - 1, lower.tail = FALSE), pchisq(chi, m - 1)))
}
# On rejette l'hypothèse nulle d'une CSR 


## avec le package spatstat #######################

nutria <- nutria_periode %>%
  filter(year == annee)

X <- as.ppp(st_coordinates(nutria), st_bbox(nutria))
Q <- quadratcount(X, nx = 100, ny = 50)
# Q
# plot(X)
# axis(1)
# axis(2)
# plot(Q, add = TRUE, cex = 2)

quadrat.test(Q) # pas significatif : pas CSR
quadrat.test(Q, "regular") # ce test est significatif : clustered pas rejeté
quadrat.test(Q, "clustered")



# Fonction K de Ripley ---------------------------------------------------------

nutria <- nutria_periode %>%
  filter(year == annee)

X <- as.ppp(st_coordinates(nutria), st_bbox(nutria))

# test CSR (PPP homogène)
envelope_K_csr <- envelope(X, Kest, nsim = 99)
plot(envelope_K_csr, main = "CSR Test with Envelopes")
# évidemment pas adapté !

# estimate intensity * thinning at each point location 
# doit être donné en unité d'aire ! donc à diviser par l'aire de la cellule...
lambda_b_hat <- function(x,y) { 
  # doit marcher en vectoriel
  pts <- st_as_sf(data.frame(x = x, y = y), coords = c("x", "y")) %>%
    st_set_crs(crs_commun)
  cells <- st_nearest_feature(pts, estim)
  return(estim$lambda_med[cells]*estim$b_med[cells]/as.numeric(grid_sf$area[cells]))
}

# quelques simulation
simu1 <- rpoispp(lambda_b_hat,  win = as.owin(occitanie))
simu2 <- rpoispp(lambda_b_hat,  win = as.owin(occitanie))

pts_simu1 <- data.frame(x1 = simu1$x, y1 = simu1$y)  %>%
  st_as_sf(coords = c("x1", "y1"), crs = crs_commun)
pts_simu2 <- data.frame(x2 = simu2$x, y2 = simu2$y) %>%
  st_as_sf(coords = c("x2", "y2"), crs = crs_commun)


ggplot() +
  geom_sf(data = pts_simu1, col = "chartreuse4") +
  geom_sf(data = pts_simu2, col = "coral") +
  geom_sf(data = nutria, col = "royalblue") + # en bleu les vraies observations
  geom_sf(data = occitanie, fill = NA)

# Fonction K et enveloppe
# faut-il mettre une correction ? (edge)
env_inhom <- envelope(
  X,
  Kinhom,
  nsim = 99,
  correction = "best",
  simulate = expression(rpoispp(lambda_b_hat, win = as.owin(occitanie)))
)
plot(env_inhom, main = "Inhomogeneous K-function with Envelopes")


# Now we're asking: 
# "Given the inhomogeneous intensity, are the points still more/less clustered than expected?"

