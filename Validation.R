################################################################################
#                              Validation                                      #
################################################################################


# Librairies utiles ------------------------------------------------------------
library(ggplot2)
library(tidyverse)
library(sf)
library(spatstat)

# Chargement -------------------------------------------------------------------

load("RData/5km2/grid_sf_5km2.RData")

load("Resultats_MCMC/5km2_sans_lasso/IEP_2021.RData")

periode = (2010:2024)

dpts_occitanie <- st_read("Data/departements-d-occitanie.shp") 
occitanie <- dpts_occitanie %>% st_union()
rm(dpts_occitanie)

nutria_periode <- st_read("Data/CEN_2025/Ragondin_rat_musque_pts_2025.shp") %>%
  mutate(year = year(as.Date(DateDebut))) %>% 
  filter(year %in% periode) %>%
  filter(NomVernacu == "Ragondin") %>%
  st_transform(crs = st_crs(occitanie)) %>%
  st_intersection(occitanie)



# Essai quadrat method (pas adapté à notre modèle) -----------------------------

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


annee <- 2021

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


# K fonction de Ripley ---------------------------------------------------------

# test complete spatial randomness (PPP homogène)
envelope_K_csr <- envelope(X, Kest, nsim = 99)
plot(envelope_K_csr, main = "CSR Test with Envelopes")
# évidemment pas adapté !


# estimate intensity at each point location
intensite <- function(x,y) {
  pt <- st_point(c(x,y))
  cell <- which(lengths(st_intersects(IEP_2021, pt))>0)
  return(IEP_2021$lambda_med[cell])
}

lambda_hat <- function(x_vec,y_vec) { 
  # doit marcher en vectoriel
  apply(matrix(c(x_vec, y_vec), ncol = 2, byrow = FALSE), 1, intensite)
}

env_inhom <- envelope(pinhom, Kinhom, nsim = 99,
                      simulate = expression(rpoispp(lambda_hat, win = as.owin(occitanie))))
plot(env_inhom, main = "Inhomogeneous K-function with Envelopes")

# Now we're asking: 
# "Given the inhomogeneous intensity, are the points still more/less clustered than expected?"

# if from a fitted model, just replace lambda_hat in rpoispp(lambda_hat)
# by the estimated lambda
