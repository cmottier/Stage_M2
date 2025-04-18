
library(MCMCvis)
library(ggplot2)
library(patchwork)
library(sf)
library(tidyverse)

## Vérification d'une année #################################

annee = 2024

out <- get(paste0("outMCMC_", annee))

MCMCsummary(out, param=c("alpha", "beta"))
MCMCtrace(out, pdf = FALSE, ind = TRUE, params = c("alpha", "beta"))


## Evolution des coefficients ###############################

periode <- 2010:2024

# extraction des coefficients
resume <- NULL
for (annee in periode) {
  out <- get(paste0("outMCMC_", annee))
  resume_out <- MCMCsummary(out) %>%
    rename(
      "lower" = "2.5%",
      "median" = "50%",
      "upper" = "97.5%"
    ) %>%
    mutate(annee = annee) %>%
    rownames_to_column("param")
  resume <- rbind(resume, resume_out)
}

save(resume, file = "Resultats_MCMC/5km2/coeff_periode_uni_prox.RData")

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
    subtitle = "Modèle à une détection, effort : distance aux accès",
    x = "",
    y = "",
    color = "Année"
  ) +
  scale_y_discrete(
    labels = c(
      "intercept_eff",
      "dist_acces",
      "intercept_int",
      "dist_eau",
      "logdensite",
      "agri",
      "preci_cum",
      "temp_min"
    )
  )

# ggsave(plot = p, "Image/periode_uni_prox.png", dpi = 600)




## Plot des intensités et probabilités #######################

# Intensité 
plot_l <- do.call(
  wrap_plots,
  lapply(2010:2018,
         function(x) {
           ggplot() +
             geom_sf(data = grid_sf, color = NA, aes(fill = get(paste0("lambda_", x)))) +
             labs(fill = "Intensité") + 
             scale_fill_viridis_c(begin = 0, end = 1) +
             labs(title = x) +
             theme_light()
         })
)

plot_l


plot_p <- do.call(
  wrap_plots,
  lapply(2010:2016,
         function(x) {
           ggplot() +
             geom_sf(data = grid_sf, color = NA, aes(fill = get(paste0("p_", x)))) +
             labs(fill = "Probabilité") + 
             scale_fill_viridis_c(begin = 0, end = 1) +
             labs(title = x) +
             theme_light()
         })
)

plot_p

