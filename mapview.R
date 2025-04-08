library(mapview)

mapview(nutria)

mapview(nutria, cex = 0.5) + mapview(grid_sf, zcol = "temp_max", alpha.region = 0.5, alpha = 0)

mapview(nutria, cex = 0.5) + mapview(grid_selec, burst = TRUE, hide = TRUE, alpha.region = 0.5, alpha = 0, legend)

# cellules particulières nb obs ragondins ou GBIF
mapview(grid_sf[grid_sf$grid_id %in% c(9533,9401,13188),])
