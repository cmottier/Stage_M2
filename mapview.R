library(mapview)

mapview(nutria)

mapview(nutria, cex = 0.5) + mapview(grid_sf, zcol = "temp_max", alpha.region = 0.5, alpha = 0)

mapview(nutria, cex = 0.5) + mapview(grid_selec, burst = TRUE, hide = TRUE, alpha.region = 0.5, alpha = 0, legend)
