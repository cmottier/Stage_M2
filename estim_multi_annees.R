################################################################################
#                           Modèle multi-année                                 #
################################################################################

# Librairies utiles ------------------------------------------------------------

library(nimble)
library(sf)
library(tidyverse)

# fonctions utiles -------------------------------------------------------------

source("code_multi_annees.R")
source("recherche_inits.R")


# Chargement et travail préalable sur grid -------------------------------------

# Chargement de la grille complète
load("grid_sf_5km2_bis.RData")

# on garde les cellules qui intersectent à au moins 5% (densité...)
grid_sf <- grid_sf %>%
  filter(as.numeric(area) > 0.05*max(as.numeric(area)))

# on renumérote les cellules de la sélection
grid_sf$grid_id <- 1:nrow(grid_sf)


# Choix de la période et du modèle ---------------------------------------------
# /!\ à modifier 

periode = 2010:2024
nb_annees <- length(periode)

donnees_utiles <- extract_data(grid_sf, periode) 

modele <- 0
# choix du modèle - type d'intercepts : 
# 0 : constant, 1 : indépendants, 2 : marche aléatoire, 3 : linéaire, 5 : effet aléatoire



# Obtention/chargement des paramètres initiaux ---------------------------------

# load("inits_modele0_50km2_2010-2024.RData")

inits <- valeurs_inits(
  code = choix_code(modele),
  constants = donnees_utiles$constants,
  data = donnees_utiles$data,
  params = params(modele),
  inits = get(paste0("inits", modele))
)

save(inits, file = paste0("inits_modele", modele, "_5km2_2010-2024.RData"))


# Lancement et sauvegarde ------------------------------------------------------

out <- estim_param(
  grid = grid_sf,
  modele = modele,
  periode = periode,
  inits = inits
  # inits = inits_list
)

save(out, file = paste0("out_modele", modele, "_5km2_2010-2024.RData"))


# Étude des résultats obtenus --------------------------------------------------

library(MCMCvis)
library(patchwork)

out$WAIC

(resum_coeffs <- MCMCsummary(out$samples))
MCMCtrace(out$samples, pdf = FALSE, ind = TRUE, params = "beta0")
MCMCplot(out$samples)
MCMCsummary(out$samples, params = c("beta0"))

save(resum_coeffs, file = paste0("resume_M", modele, ".RData"))

# cartes

int_eff <- grid_sf %>%
  select(grid)

# Intensité

res <- rbind(out$samples2$chain1, out$samples2$chain2)
mask <- str_detect(colnames(res), "lambda")
res_lambda <- res[,mask]
res_lambda <- apply(res_lambda, 2, median)
max_lambda <- max(res_lambda)
res_lambda <- matrix(res_lambda, ncol = length(periode), byrow = FALSE)
for (a in 1:length(periode)) {
  int_eff[[paste0("lambda_med", periode[a]) ]] <- res_lambda[,a]
}

plot_l <- do.call(
  wrap_plots,
  list(lapply(periode,
              function(x) {
                ggplot() +
                  geom_sf(data = int_eff, color = NA, aes(fill = get(paste0("lambda_med", x)))) +
                  labs(fill = "Intensité") +
                  scale_fill_viridis_c(begin = 0, end = 1,
                                       limits = c(0, max_lambda)
                  ) + 
                  labs(title = x) +
                  theme_light()
              })  
       , guides = 'collect'
  )
) +
  plot_annotation(title = "Intensités - Modèle multi-années, intercepts indépendants") 


plot_l

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


# Effort

mask <- str_detect(colnames(res), "b[^d]")
res_b <- res[,mask]
res_b <- apply(res_b, 2, median)
res_b <- matrix(res_b, ncol = length(periode), byrow = FALSE)
for (a in 1:length(periode)) {
  int_eff[[paste0("b_med", periode[a]) ]] <- res_b[,a]
}

plot_b <- do.call(
  wrap_plots,
  list(lapply(periode,
              function(x) {
                ggplot() +
                  geom_sf(data = int_eff, color = NA, aes(fill = get(paste0("b_med", x)))) +
                  labs(fill = "Effort") +
                  scale_fill_viridis_c(begin = 0, end = 1, limits = c(0, 1)) + # pour fixer couleurs
                  labs(title = x) +
                  theme_light()
              }),
       guides = 'collect')
) +
  plot_annotation(title = "Effort - Modèle multi-années, intercepts indépendants") 

plot_b

# probabilité de présence

mask <- str_detect(colnames(res), "p")
res_p <- res[,mask]
res_p <- apply(res_p, 2, median)
res_p <- matrix(res_p, ncol = length(periode), byrow = FALSE)
for (a in 1:length(periode)) {
  int_eff[[paste0("p_med", periode[a]) ]] <- res_p[,a]
}

plot_p <- do.call(
  wrap_plots,
  list(lapply(periode,
              function(x) {
                ggplot() +
                  geom_sf(data = int_eff, color = NA, aes(fill = get(paste0("p_med", x)))) +
                  labs(fill = "Probabilité") +
                  scale_fill_viridis_c(begin = 0, end = 1, limits = c(0, 1)) + # pour fixer couleurs
                  labs(title = x) +
                  theme_light()
              }),
       guides = 'collect')
) +
  plot_annotation(title = "Probabilité de présence - Modèle multi-années, intercepts indépendants") 

plot_p


