################################################################################
#                Intensité, effort, probabilité de présence                    #
################################################################################

# Chargement des librairies ----------------------------------------------------
library(sf)
library(tidyverse)

# Fonctions utiles -------------------------------------------------------------

# fonction d'intensité
lambda <- function(coeff, var, logarea) {
  return(exp(sum(coeff*var) + logarea))
}

# fonction d'effort
effort <- function(coeff, var) {
  return(plogis(sum(coeff*var)))
}

# Transformation des MCMC

#' Indicateurs statistiques de l'intensité, l'effort et 
#' la probabilité de présence pour chaque cellule de la grille
#'
#' @param grid # table contenant la grille et toutes les variables associées
#' @param out # contient les 2 chaînes de Markov obtenue pour une année
#' @param annee # l'année associée aux chaînes de Markov
#' 
intensite_effort_proba <- function(grid, out, annee) {
  
  grid <- grid %>%
    mutate(dist_eau_scaled = scale(dist_eau)[,1]) %>%
    mutate(logdensity_scaled = scale(logdensity)[,1]) %>%
    mutate(agri_cover_scaled = scale(agri_cover)[,1]) %>%
    mutate(across(starts_with("pcum"), ~scale(.x)[,1], .names = "scaled_{.col}")) %>%
    mutate(across(starts_with("tmin"), ~scale(.x)[,1], .names = "scaled_{.col}")) %>%
    mutate(across(starts_with("dgbif"), ~scale(.x)[,1], .names = "scaled_{.col}")) %>%
    mutate(logarea = log(as.numeric(units::set_units(area,"km^2")))) 
  
  grid_out <- grid %>%
    select(c(grid_id, grid))

  # Variables utiles
  var_int <- cbind(rep(1, nrow(grid)),
               grid$dist_eau_scaled,
               grid$logdensity_scaled,
               grid$agri_cover_scaled,
               grid[[paste0("scaled_pcum_", annee)]],
               grid[[paste0("scaled_tmin_", annee)]])
  
  var_eff <- cbind(rep(1, nrow(grid)),
              grid[[paste0("scaled_dgbif_", annee)]])
  
  # MCMC 
  res <- rbind(out$chain1, out$chain2)
  
  # coefficients beta
  mask_beta <- str_detect(colnames(res), "beta")
  coeff_beta <- res[,mask_beta]
  
  # coefficients alpha
  mask_alpha <- str_detect(colnames(res), "alpha")
  coeff_alpha <- res[,mask_alpha]
  
  for (cell in 1:nrow(grid)) {
    print(cell)

    # intensité
    lambda_simu_cell <- apply(coeff_beta, 1, lambda, var = var_int[cell,], logarea = grid$logarea[cell])
    grid_out$lambda_med[cell] <- median(lambda_simu_cell)
    grid_out$lambda_mean[cell] <- mean(lambda_simu_cell)
    grid_out$lambda_sd[cell] <- sd(lambda_simu_cell)
    # grid_out$lambda_2.5[cell] <- quantile(lambda_simu_cell, 0.25)
    # grid_out$lambda_97.5[cell] <- quantile(lambda_simu_cell, 0.975)

    # effort
    b_simu_cell <- apply(coeff_alpha, 1, effort, var = var_eff[cell,])
    grid_out$b_med[cell] <- median(b_simu_cell)
    grid_out$b_mean[cell] <- mean(b_simu_cell)
    grid_out$b_sd[cell] <- sd(b_simu_cell)
    # grid_out$b_2.5[cell] <- quantile(b_simu_cell, 0.25)
    # grid_out$b_97.5[cell] <- quantile(b_simu_cell, 0.975)

    # probabilité d'occupation
    prob_simu_cell <- 1-exp(-lambda_simu_cell)
    grid_out$prob_med[cell] <- median(prob_simu_cell)
    grid_out$prob_mean[cell] <- mean(prob_simu_cell)
    grid_out$prob_sd[cell] <- sd(prob_simu_cell)
    # grid_out$prob_2.5[cell] <- quantile(prob_simu_cell, 0.25)
    # grid_out$prob_97.5[cell] <- quantile(prob_simu_cell, 0.975)
    }

  return(grid_out)
}

## Application -----------------------------------------------------------------

load("out_multi_gbif_2016.RData")
load("grid_sf_5km2.RData")

start <- Sys.time()
IEP_2016 <- intensite_effort_proba(grid_sf, out, 2016)
end <- Sys.time()
print(end - start)

save(IEP_2016, file = "IEP_multi_gbif_l_2016.RData")

# ggplot() +
#   geom_sf(data = IEP_2021, aes(fill = lambda_med))


