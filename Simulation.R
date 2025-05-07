################################################################################
#                                Simulations                                   #
################################################################################

# Librairies utiles ------------------------------------------------------------
library(ggplot2)
library(tidyverse)
# library(sf)
# library(spatstat)
# library(AHMbook)
library(nimble)
# library(MCMCvis)
library(plot.matrix)




# Avec AHMbook -----------------------------------------------------------------

# Simulation d'un PPP aminci
# On utilise ici la variable d'environnement x et la variable d'effort w 
# déterminées dans l'annexe de Koshkina

set.seed(123)

dat <- simDataDK(
  sqrt.npix = 100, 
  alpha = c(-2, -1), # effort : b = -2-1*w
  beta = c(6, 1), # intensité : l = 6+1*x
  drop.out.prop.pb = 0, 
  quadrat.size = 4,
  # gamma = c(0, -1.5), 
  # nquadrats = 250, 
  # nsurveys = 0,
  show.plot = FALSE
  )

str(data, 1)

# Our landscape is inhabited by a total of ... individuals (dat$N.ipp)
dat$N.ipp

# of which ... are detected
dat$N.det

#Their locations represent our presence-only data: 
# loc.det contains the coordinates 
head(dat$loc.det) # Actual coordinates of points detected

# and pixel.id.det the pixel ID for each detected point.
head(dat$pixel.id.det) # Pixel ID of 302 points detected

# le thinning est donc de :
1-dat$N.det/dat$N.ipp

# surface des cellules
(logarea <- log(dat$s.area / dat$npix))



ggplot() +
  geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = dat$xcov))+
  geom_point(data = dat$loc.ipp, aes(x = x, y = y), col = "white") +
  scale_fill_viridis_c() +
  labs(fill = "x", x = "", y = "")

ggplot()+
  geom_raster(data = dat$s.loc, aes(x = x, y = y, fill = dat$wcov))+
  geom_point(data = dat$loc.det, aes(x = x, y = y), col = "white") +
  scale_fill_viridis_c() +
  labs(fill = "h", x = "", y = "")

