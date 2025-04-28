################################################################################
#                Intensité, effort, probabilité de présence                    #
################################################################################
library(sf)
library(tidyverse)

lambda <- function(coeff, var, logarea) {
  return(exp(sum(coeff*var) + logarea))
}

intensite <- function(grid, out, annee) {
  
  grid <- grid %>%
    mutate(dist_eau_scaled = scale(dist_eau)[,1]) %>%
    mutate(logdensity_scaled = scale(logdensity)[,1]) %>%
    mutate(agri_cover_scaled = scale(agri_cover)[,1]) %>%
    mutate(across(starts_with("pcum"), ~scale(.x)[,1], .names = "scaled_{.col}")) %>%
    mutate(across(starts_with("tmin"), ~scale(.x)[,1], .names = "scaled_{.col}")) %>%
    mutate(logarea = log(as.numeric(units::set_units(area,"km^2"))))
  
  grid_out <- grid %>%
    select(c(grid_id, grid))

  var <- cbind(rep(1, nrow(grid)),
               grid$dist_eau_scaled,
               grid$logdensity_scaled,
               grid$agri_cover_scaled,
               grid[[paste0("scaled_pcum_", annee)]],
               grid[[paste0("scaled_tmin_", annee)]])
  
  res <- rbind(out$chain1, out$chain2)
  
  # coefficients beta
  mask_beta <- str_detect(colnames(res), "beta")
  coeff_beta <- res[,mask_beta]
  
  for (cell in 1:nrow(grid)) {
    print(cell)
    lambda_simu_cell <- apply(coeff_beta, 1, lambda, var = var[cell,], logarea = grid$logarea[cell])
    grid_out$lambda_med[cell] <- median(lambda_simu_cell)
    grid_out$lambda_2.5[cell] <- quantile(lambda_simu_cell, 0.25)
    grid_out$lambda_97.5[cell] <- quantile(lambda_simu_cell, 0.975)
    }
  return(grid_out)
}

test <- intensite(grid_sf[1:500,], outMCMC_2021, 2021)

ggplot() +
  geom_sf(data = test, aes(fill = lambda_med))

# 
# 
# for (i in 1:nrow(coeff_beta)) {
#   lambda_cell <- NULL
#   lambda_cell[i] <- exp(
#     coeff_beta[i, 1] +
#       coeff_beta[i, 2] * grid$dist_eau_scaled[cell] +
#       coeff_beta[i, 3] * grid$logdensity_scaled[cell] +
#       coeff_beta[i, 4] * grid$agri_cover_scaled[cell] +
#       coeff_beta[i, 5] * grid[[paste0("scaled_pcum_", annee)]][cell] +
#       coeff_beta[i, 6] * grid[[paste0("scaled_tmin_", annee)]][cell] +
#       grid$logarea[cell]
#   )
