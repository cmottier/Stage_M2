################################################################################
#                          Comparaison des modèles                            #
################################################################################

# Librairies utiles ------------------------------------------------------------

# library(nimble)
library(sf)
library(tidyverse)

# Chargement et travail préalable sur grid -------------------------------------

# Chargement de la grille complète
load("RData/5km2/grid_sf_5km2_bis.RData")

# on garde les cellules qui intersectent à au moins 5% (densité...)
grid_sf <- grid_sf %>%
  filter(as.numeric(area) > 0.05*max(as.numeric(area)))

# on renumérote les cellules de la sélection
grid_sf$grid_id <- 1:nrow(grid_sf)

# Période ----------------------------------------------------------------------

periode = 2010:2024
nb_annees <- length(periode)


# Comparaison ------------------------------------------------------------------

load("Resultats_MCMC/5km2/Annee_par_annee/resume_multi_gbif_l.RData")
# resume contient les coefficients des modèles année par année


## Comparaison des coefficients fixes ##########################################

# sélection des coefficients dans les modèles année par année
resume_global_coeffs_fixes <- resume %>%
  filter(param %in% c(paste0("beta[", c(2:6),"]"), "alpha[2]"))

# sélection des coefficients dans les gros modèles
for (modele in c(0,1,3)) {
  load(paste0("Resultats_MCMC/5km2/Gros_modele/resume_M", modele, ".RData"))
  
  resum_coeffs <- resum_coeffs%>%
    rename(
      "lower" = "2.5%",
      "median" = "50%",
      "upper" = "97.5%"
    )
  
  res_M <- resum_coeffs[16:21,]
  res_M$param <- c("alpha[2]", paste0("beta[", c(2:6),"]"))
  res_M$annee <- rep(paste0("M", modele), 6)
  
  # on regroupe tout
  resume_global_coeffs_fixes <- rbind(resume_global_coeffs_fixes, res_M)
}

# choix des couleurs 
palette_auto <- scales::hue_pal()(nb_annees)
names(palette_auto) <- periode

palette_auto["M0"] <- "darkmagenta"
palette_auto["M1"] <- "darkblue"
palette_auto["M3"] <- "darkgreen"

# plot 
p <- ggplot(data = resume_global_coeffs_fixes,
            aes(
              y = param,
              x = median,
              xmin = lower,
              xmax = upper,
              color = as.factor(annee)
            )) +
  geom_vline(aes(xintercept = 0)) +
  geom_pointrange(position = position_dodge(width = .8)) +
  labs(
    title = "Comparaison des coefficients des différents modèles",
    x = "",
    y = "",
    color = "Année"
  ) +
  scale_y_discrete(
    labels = c(
      "densité GBIF",
      "dist_eau",
      "logdensite",
      "agri",
      "preci_cum",
      "temp_min"
    )
  ) + 
  scale_color_manual(values = palette_auto)
p


## Intercepts d'intensité ######################################################

# sélection des intercepts dans les modèles année par année
resume_global_intercept_int <- resume %>%
  select(c("upper", "median", "lower", "annee", "param")) %>%
  filter(param == "beta[1]") %>%
  select(-param)

resume_global_intercept_int$modele <- "annee par annee"

# sélection du coefficient dans le modèle à intercep constant
load("Resultats_MCMC/5km2/Gros_modele/resume_M0.RData")

resum_beta0 <- resum_coeffs %>%
  rownames_to_column(var = "id") %>%       # transforme les rownames en une colonne "id"
  filter(id == "beta0") %>%           # filtre sur la colonne "id"
  column_to_rownames(var = "id") %>%
  rename(
    "lower" = "2.5%",
    "median" = "50%",
    "upper" = "97.5%"
  ) %>%
  select(c("lower", "median", "upper"))

for (annee in periode) {
  resum_beta0$annee <- annee
  resum_beta0$modele <- "M0"
  resume_global_intercept_int <- rbind(resume_global_intercept_int, resum_beta0)
}

# sélection des intercepts du modèle à intercepts indépendants
load("Resultats_MCMC/5km2/Gros_modele/resume_M1.RData")

resum_beta1 <- resum_coeffs %>%
  rownames_to_column(var = "id") %>%       # transforme les rownames en une colonne "id"
  filter(str_starts(id, "beta0")) %>%   
  column_to_rownames(var = "id") %>%
  rename(
    "lower" = "2.5%",
    "median" = "50%",
    "upper" = "97.5%"
  )  %>%
  select(c("lower", "median", "upper"))

resum_beta1$annee <- periode
resum_beta1$modele <- "M1"
resume_global_intercept_int <- rbind(resume_global_intercept_int, resum_beta1)

# sélection des intercepts dans le modèle à effet linéaire
load("Resultats_MCMC/5km2/Gros_modele/resume_M3.RData")

resum_M3 <- resum_coeffs %>%
  rownames_to_column(var = "id") %>%       # transforme les rownames en une colonne "id"
  filter(id %in% c("beta0", "c")) %>%           # filtre sur la colonne "id"
  column_to_rownames(var = "id") %>%
  rename(
    "lower" = "2.5%",
    "median" = "50%",
    "upper" = "97.5%"
  ) %>%
  select(c("lower", "median", "upper"))

for (a in 1:nb_annees) {
  res <- resum_M3[1,] + a * resum_M3[2,]
  res$annee <- periode[a]
  res$modele <- "M3"
  res$upper <- res$median
  res$lower <- res$median
  resume_global_intercept_int <- rbind(resume_global_intercept_int, res)
}

# couleur de fond
bandes_fond <- data.frame(
  annee = periode,
  ymin = seq_along(periode) - 0.5,
  ymax = seq_along(periode) + 0.5,
  fill = rep(c("gray70", "gray90"), length.out = nb_annees)
)

# plot 
p <- ggplot(data = resume_global_intercept_int,
            aes(
              y = as.factor(annee),
              x = median,
              xmin = lower,
              xmax = upper,
              color = modele
            )) +
  geom_rect(data = bandes_fond,
            aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = fill),
            inherit.aes = FALSE,
            alpha = 0.4,
            show.legend = FALSE) +
  geom_pointrange(position = position_dodge(width = .8)) +
  scale_fill_identity() +   
  labs(
    title = "Evolution des intercepts d'intensité",
    x = "",
    y = "",
    color = "Modèle"
  ) +
  scale_color_discrete(
    labels = c(
      "M0" = "intercept constant",
      "M1" = "intercepts indépendants",
      "M3" = "effet linéaire"
    )
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),  # enlève les grandes lignes horizontales
    panel.grid.minor.y = element_blank()   # enlève les petites lignes horizontales
  )
p




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
