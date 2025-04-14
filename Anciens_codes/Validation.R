library("spatstat")

X <- as.ppp(st_coordinates(nutria), st_bbox(nutria))
plot(X)

K <- envelope(X, fun  = Kinhom, nsim = 1001, nrank = 26)
plot(K)
