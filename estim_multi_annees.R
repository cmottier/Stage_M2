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

modele <- 3
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

# library(MCMCvis)
# library(patchwork)
# 
# out$WAIC
# 
# (resum_coeffs <- MCMCsummary(out$samples))
# MCMCtrace(out$samples, pdf = FALSE, ind = TRUE, params = "beta0")
# MCMCplot(out$samples)
# MCMCsummary(out$samples, params = c("beta0"))
# 
# # cartes
# 
# int_eff <- grid_sf %>%
#   select(grid)
# 
# # Intensité
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
# 
# plot_l <- do.call(
#   wrap_plots,
#   lapply(periode,
#          function(x) {
#            ggplot() +
#              geom_sf(data = int_eff, color = NA, aes(fill = get(paste0("lambda_med", x)))) +
#              labs(fill = "Intensité") +
#              scale_fill_viridis_c(begin = 0, end = 1,
#                                   # limits = c(0, max_lambda)
#                                   ) + # pour fixer couleurs
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
# 
# mask <- str_detect(colnames(res), "b[^d]")
# res_b <- res[,mask]
# res_b <- apply(res_b, 2, median)
# res_b <- matrix(res_b, ncol = length(periode), byrow = FALSE)
# for (a in 1:length(periode)) {
#   int_eff[[paste0("b_med", periode[a]) ]] <- res_b[,a]
# }
# 
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
# 
# ## Comparaison avec modèles annuels ############################################
# 
# # coefficients gros_modele
# resum_coeffs <- resum_coeffs%>%
#   rename(
#     "lower" = "2.5%",
#     "median" = "50%",
#     "upper" = "97.5%"
#   )
# 
# # coefficients année par année
# resume <- NULL
# for (annee in periode) {
#   # out <- get(paste0("outMCMC_", annee))
#   load(paste0("Resultats_MCMC/5km2/Avec_lasso/out_multi_gbif_l_iep", annee, ".RData"))
#   resume_out <- MCMCsummary(out$samples) %>%
#     rename(
#       "lower" = "2.5%",
#       "median" = "50%",
#       "upper" = "97.5%"
#     ) %>%
#     mutate(annee = annee) %>%
#     rownames_to_column("param")
#   resume <- rbind(resume, resume_out)
# }
# 
# # comparaison des écart-types (gros modèle - année par année)
# # sur l'intercept d'intensité
# resum_coeffs[22:36,] %>% select(sd) < resume %>% filter(param == "beta[1]") %>% select(sd)
# # sur l'intercept d'effort
# resum_coeffs[1:15,] %>% select(sd) < resume %>% filter(param == "alpha[1]") %>% select(sd)
# 
# # intercept (intensité)
# ggplot() +
#   # gros modele
#   geom_ribbon(data = resum_coeffs[22:36,],
#               aes(ymin = lower, ymax = upper, x = 2010:2024),
#               fill = "royalblue",
#               alpha = 0.2) +
#   geom_line(data = resum_coeffs[22:36,],
#             color = "royalblue",
#             aes(y = median, x = 2010:2024)) +
#   # année par année
#   geom_ribbon(data = resume %>% filter(param == "beta[1]"),
#               aes(ymin = lower, ymax = upper, x = 2010:2024),
#               fill = "chartreuse4",
#               alpha = 0.2) +
#   geom_line(data = resume %>% filter(param == "beta[1]"),
#             aes(y = median, x = 2010:2024),
#             color = "chartreuse4")
# 
# # intercept effort
# ggplot() +
#   # gros modèle
#   geom_ribbon(data = resum_coeffs[1:15,],
#               aes(ymin = lower, ymax = upper, x = 2010:2024),
#               fill = "royalblue",
#               alpha = 0.2) +
#   geom_line(data = resum_coeffs[1:15,],
#             color = "royalblue",
#             aes(y = median, x = 2010:2024)) +
#   # modèle année par année
#   geom_ribbon(data = resume %>% filter(param == "alpha[1]"),
#               aes(ymin = lower, ymax = upper, x = 2010:2024),
#               fill = "chartreuse4",
#               alpha = 0.2) +
#   geom_line(data = resume %>% filter(param == "alpha[1]"),
#             aes(y = median, x = 2010:2024),
#             color = "chartreuse4")
# 
# resum_coeffs[19,]
# resume %>% filter(param == "beta[4]")
