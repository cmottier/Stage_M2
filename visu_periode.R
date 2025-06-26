
library(MCMCvis)
library(ggplot2)
library(patchwork)
library(sf)
library(tidyverse)

periode <- 2010:2024

## Observation d'une variable #################################

v <- "tmin_" 
v <- "pcum_"
v <- "dgbif_"

do.call(
  wrap_plots,
  lapply(periode,
         function(x) {
           ggplot() +
             geom_sf(data = grid_sf, color = NA, aes(fill = as.numeric(get(paste0(v, x))))) +
             labs(fill = v) + 
             scale_fill_viridis_c(begin = 0, end = 1) +# limits = c(0, 10)) pour fixer couleurs
             labs(title = x) +
             theme_light()
         })
)

v <- "agri_cover"  
v <- "density"
v <- "logdensity"   
v <- "dist_acces"     
v <- "dist_eau" 

ggplot() +
  geom_sf(data = grid_sf, color = NA, aes(fill = as.numeric(.data[[v]]))) +
  labs(v) + 
  scale_fill_viridis_c(begin = 0, end = 1) +# limits = c(0, 10)) pour fixer couleurs
  theme_light()


## Evolution des coefficients ###############################

# extraction des coefficients
resume <- NULL
for (annee in periode) {
  # out <- get(paste0("outMCMC_", annee))
  load(paste0("Resultats_MCMC/5km2/Annee_par_annee/out_multi_gbif_l_", annee, ".RData"))
  resume_out <- MCMCsummary(out$samples) %>%
    rename(
      "lower" = "2.5%",
      "median" = "50%",
      "upper" = "97.5%"
    ) %>%
    mutate(annee = annee) %>%
    rownames_to_column("param")
  resume <- rbind(resume, resume_out)
}

# save(resume, file = "Resultats_MCMC/5km2/Annee_par_annee/resume_multi_gbif_l.RData")

# plot 
p <- ggplot(data = resume,
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
    title = "Evolution des coefficients",
    subtitle = "Modèle à multiples détections, effort : données GBIF",
    x = "",
    y = "",
    color = "Année"
  ) +
  scale_y_discrete(
    labels = c(
      "intercept_eff",
      "densité GBIF",
      "intercept_int",
      "dist_eau",
      "logdensite",
      "agri",
      "preci_cum",
      "temp_min"
    )
  )
p
# ggsave(plot = p, "Image/periode_uni_prox.png", dpi = 600)

# pour rajouter la somme des intercepts : 

somme_df <- resume %>%
  filter(param %in% c("alpha[1]", "beta[1]")) %>%
  group_by(annee) %>%
  summarise(
    param = "somme intercepts",
    median = sum(median),
    lower = sum(median), # pas de valeur
    upper = sum(median), # pas de valeur
    .groups = "drop"
  )

# Ajouter les lignes au data.frame principal
resume_extended <- bind_rows(resume, somme_df)

p <- ggplot(data = resume_extended,
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
    title = "Evolution des coefficients",
    subtitle = "Modèle à multiples détections, effort : données GBIF",
    x = "",
    y = "",
    color = "Année"
  ) +
  scale_y_discrete(
    labels = c(
      "intercept eff",
      "densité GBIF",
      "intercept int",
      "dist eau",
      "logdensite",
      "agri",
      "preci cum",
      "temp min", 
      "somme intercepts"
    )
  )
p

## Distributions des intensités, efforts, probabilités ##############

annee = 2019

out <- get(paste0("outMCMC_", annee))

MCMCsummary(out, param=c("alpha", "beta"))
MCMCtrace(out, pdf = FALSE, ind = TRUE, params = c("alpha", "beta"))




## Plot des intensités et probabilités #######################

# for (a in periode) {
#   coeffs <- resume %>%
#     filter(annee == a) %>%
#     filter(str_detect(param, "beta")) %>%
#     select(param, median)
#   assign(paste0("lambda_", a),
#          exp(coeffs$median[1] +
#              coeffs$median[2] * scale(grid_sf$dist_eau)[,1] +
#              coeffs$median[3] * scale(grid_sf$logdensity)[,1] +
#              coeffs$median[4] * scale(grid_sf$agri_cover)[,1] +
#              coeffs$median[5] * scale(grid_sf[[paste0("pcum_", a)]])[,1] +
#              coeffs$median[6] * scale(grid_sf[[paste0("tmin_", a)]])[,1] +
#              log(as.numeric(units::set_units(grid_sf$area,"km^2")))))
#   assign(paste0("p_", a),
#          1-exp(-get(paste0("lambda_",a))))
# }

iep <- grid_sf  %>%
  filter(as.numeric(area) > 0.05*max(as.numeric(area))) %>%
  select(grid)

for (annee in periode) {
  # out <- get(paste0("outMCMC_", annee))
  load(paste0("Resultats_MCMC/5km2/Annee_par_annee/out_multi_gbif_l_", annee, ".RData"))
  res <- rbind(out$samples2$chain1, out$samples2$chain2)
  mask <- str_detect(colnames(res), "lambda")
  res_lambda <- res[,mask]
  iep[[paste0("lambda_med", annee) ]] <- apply(res_lambda, 2, median)
  iep[[paste0("lambda_sd", annee) ]] <- apply(res_lambda, 2, sd)
  mask <- str_detect(colnames(res), "b[^ed]")
  res_b <- res[,mask]
  iep[[paste0("b_med", annee) ]] <- apply(res_b, 2, median)
  iep[[paste0("b_sd", annee) ]] <- apply(res_b, 2, sd)
  mask <- str_detect(colnames(res), "p")
  res_p <- res[,mask]
  iep[[paste0("p_med", annee) ]] <- apply(res_p, 2, median)
  iep[[paste0("p_sd", annee) ]] <- apply(res_p, 2, sd)
}

# save(iep, file = "Resultats_MCMC/5km2/Annee_par_annee/iep_multi_gbif_l.RData")

# Intensité 

lambda_max <- iep %>% select(starts_with("lambda_med")) %>% st_drop_geometry %>% max()

plot_l <- do.call(
  wrap_plots,
  list(lapply(2010:2024,
         function(x) {
           ggplot() +
             geom_sf(data = iep, color = NA, aes(fill = get(paste0("lambda_med", x)))) +
             labs(fill = "Intensité") + 
             scale_fill_viridis_c(begin = 0, end = 1) + #, limits = c(0, lambda_max)) + 
             labs(title = x) +
             theme_light()
         })
  # , guides = 'collect'
  )
) +
  plot_annotation(title = "Intensités - année par année") 

plot_l

sd_l_max <- iep %>% select(starts_with("lambda_sd")) %>% st_drop_geometry %>% max()

plot_l_sd <- do.call(
  wrap_plots,
  lapply(periode,
         function(x) {
           ggplot() +
             geom_sf(data = iep, color = NA, aes(fill = get(paste0("lambda_sd", x)))) +
             labs(fill = "sd Int") + 
             scale_fill_viridis_c(begin = 0, end = 1) +# limits = c(0, 10)) pour fixer couleurs
           labs(title = x) +
             theme_light()
         })
)

plot_l_sd

# effort 

b_max <- iep %>% select(starts_with("b_med")) %>% st_drop_geometry %>% max()

plot_b <- do.call(
  wrap_plots,
  list(lapply(periode,
              function(x) {
                ggplot() +
                  geom_sf(data = iep, color = NA, aes(fill = get(paste0("b_med", x)))) +
                  labs(fill = "Effort") + 
                  scale_fill_viridis_c(begin = 0, end = 1, limits = c(0,1)) + #, limits = c(0, lambda_max)) + 
                  labs(title = x) +
                  theme_light()
              })
       , guides = 'collect'
  )
) +
  plot_annotation(title = "Effort - année par année") 

plot_b


# probabilité de présence

plot_p <- do.call(
  wrap_plots,
  list(lapply(periode,
              function(x) {
                ggplot() +
                  geom_sf(data = iep, color = NA, aes(fill = get(paste0("p_med", x)))) +
                  labs(fill = "Probabilité") + 
                  scale_fill_viridis_c(begin = 0, end = 1, limits = c(0,1)) + 
                  labs(title = x) +
                  theme_light()
              })
       , guides = 'collect'
  )
) +
  plot_annotation(title = "Probabilité de présence - année par année") 

plot_p


