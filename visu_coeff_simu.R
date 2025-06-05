# Pour visualiser les coefficients obtenus lors des 50 simulations d'un même IPP

library(ggplot2)

num <- 2

load(paste0("RData/Simulations/simulation_scenario",num,".RData"))

# coefficients utilisés et nombre de variables
if (num ==1) {
  # scenario1
  alpha <- c(0, -5) 
  beta <- c(6, 1)
} else if (num == 2) {
  # scenario2
  alpha <- c(-2, 2) 
  beta <- c(-4, 4, 1)
} else if (num == 3) {
  # scenario 3
  alpha <- c(-4, 1) 
  beta <- c(-2, 4, 1)
} else if (num == 4) {
  # scenario 4
  alpha <- c(-2, 1) 
  beta <- c(-2, 1, 1)
}

nvar <- length(beta)-1
    

# on transforme en dataframe
df <- NULL
res <- unlist(resultat)
df <- data.frame(
  low = res[names(res) %in% c("avec1", "avec4", "avec7", "avec10", "avec13", "sans1", "sans4", "sans7")],
  med = res[names(res) %in% c("avec2", "avec5", "avec8", "avec11", "avec14", "sans2", "sans5", "sans8")], 
  up = res[names(res) %in% c("avec3", "avec6", "avec9", "avec12", "avec15", "sans3", "sans6", "sans9")],
  vrai = rep(c(alpha, beta, beta), 50)
)

if (nvar == 1) {
  df$param <- rep(c("a1", "a2", "b1", "b2", "b1", "b2"), 50)
  df$effort <- rep(c(rep("avec", 4), rep("sans", 2)), 50)
} else {
  df$param <- rep(c("a1", "a2", "b1", "b2", "b3", "b1", "b2", "b3"), 50)
  df$effort <- rep(c(rep("avec", 5), rep("sans", 3)), 50)
}

ggplot(data = df, aes(y = 1:(50*(2+(nvar+1)*2)), x = med, xmin = low, xmax = up)) +
  facet_grid(effort ~ param, scales='free') + #
  geom_pointrange(position = position_dodge(width = .8)) +
  geom_vline(aes(xintercept = vrai), col = "red") +
  labs(x = "", y = "", title = paste0("scenario ", num))
             
        