################################################################################
#                           Modèle multi-année                                 #
################################################################################

# Librairies utiles ------------------------------------------------------------

library(nimble)
library(sf)
library(tidyverse)

# fonctions utiles -------------------------------------------------------------

source("code_multi_annees.R")

# Chargement  de la grille -----------------------------------------------------
# /!\ à modifier

load("grid_sf_5km2.RData")

# Chargement des paramètres initiaux -------------------------------------------
# /!\ à modifier

load("inits_modele0_5km2_2010-2024.RData")

# Lancement et sauvegarde ------------------------------------------------------
# /!\ à modifier

periode = 2010:2024
nb_annees <- length(periode)
modele <- 0

out <- estim_param(
  grid = grid_sf,
  modele = modele,
  periode = periode,
  inits = inits
  # inits = inits_list
)

save(out, file = "modele2_50km2_2010-2024.RData")


# Étude du résultat obtenu -----------------------------------------------------

# library(MCMCvis)
# library(patchwork)
# 
# out$WAIC
# 
# (resume <- MCMCsummary(out$samples))
# MCMCtrace(out$samples, pdf = FALSE, ind = TRUE)
# MCMCplot(out$samples)
# 
# 
# int_eff <- grid_sf %>%
#   select(grid)
# 
# res <- rbind(out$samples2$chain1, out$samples2$chain2)
# mask <- str_detect(colnames(res), "lambda")
# res_lambda <- res[,mask]
# res_lambda <- apply(res_lambda, 2, median)
# max_lambda <- max(res_lambda)
# res_lambda <- matrix(res_lambda, ncol = length(periode), byrow = FALSE)
# for (a in 1:length(periode)) {
#   int_eff[[paste0("lambda_med", periode[a]) ]] <- res_lambda[,a]
# }
# mask <- str_detect(colnames(res), "b[^d]")
# res_b <- res[,mask]
# res_b <- apply(res_b, 2, median)
# res_b <- matrix(res_b, ncol = length(periode), byrow = FALSE)
# for (a in 1:length(periode)) {
#   int_eff[[paste0("b_med", periode[a]) ]] <- res_b[,a]
# }
# 
# # Intensité
# plot_l <- do.call(
#   wrap_plots,
#   lapply(periode,
#          function(x) {
#            ggplot() +
#              geom_sf(data = int_eff, color = NA, aes(fill = get(paste0("lambda_med", x)))) +
#              labs(fill = "Intensité") +
#              scale_fill_viridis_c(begin = 0, end = 1, limits = c(0, max_lambda)) + # pour fixer couleurs
#            labs(title = x) +
#              theme_light()
#          })
# )
# 
# plot_l
# 
# # Intensité normalisée
# plot_l_scaled <- do.call(
#   wrap_plots,
#   lapply(periode,
#          function(x) {
#            ggplot() +
#              geom_sf(data = int_eff, color = NA, aes(fill = scale(get(paste0("lambda_med", x)))[,1])) +
#              labs(fill = "Intensité") +
#              scale_fill_viridis_c(begin = 0, end = 1, limits = c(-1.5,5)) + # pour fixer couleurs
#              labs(title = x) +
#              theme_light()
#          })
# )
# 
# plot_l_scaled
# 
# # effort
# plot_b <- do.call(
#   wrap_plots,
#   lapply(periode,
#          function(x) {
#            ggplot() +
#              geom_sf(data = int_eff, color = NA, aes(fill = get(paste0("b_med", x)))) +
#              labs(fill = "Effort") +
#              scale_fill_viridis_c(begin = 0, end = 1, limits = c(0, 1)) + # pour fixer couleurs
#              labs(title = x) +
#              theme_light()
#          })
# )
# 
# plot_b
# 
