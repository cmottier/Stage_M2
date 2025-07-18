# Pour visualiser les coefficients obtenus lors des 50 simulations d'un même IPP
library(tidyverse)
library(ggplot2)

num <- 3

load(paste0("RData/Simulations/simulation_scenario",num,"_bis.RData"))

n <- 3 # nb de simu

# coefficients utilisés et nombre de variables
if (num ==1) {
  # scenario1
  alpha <- c(0, -5) 
  beta <- c(6, 1)
} else if (num == "1bis") {
  # scenario2
  alpha <- c(-1, -5) 
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
} else if (num == "4ter") {
  # scenario 4
  alpha <- c(-2,2) 
  beta <- c(-6,2,2)
}

nvar <- length(beta)-1
    

# on transforme en dataframe
df <- NULL
res <- unlist(resultat)
df <- data.frame(
  low = res[names(res) %in% c("avec1", "avec4", "avec7", "avec10", "avec13", "sans1", "sans4", "sans7")],
  med = res[names(res) %in% c("avec2", "avec5", "avec8", "avec11", "avec14", "sans2", "sans5", "sans8")], 
  up = res[names(res) %in% c("avec3", "avec6", "avec9", "avec12", "avec15", "sans3", "sans6", "sans9")],
  vrai = rep(c(alpha, beta, beta), n)
)

if (nvar == 1) {
  df$param <- rep(c("a0", "a1", "b0", "b1", "b0", "b1"), n)
  df$effort <- rep(c(rep("avec effort", 4), rep("sans effort", 2)), n)
} else {
  df$param <- rep(c("a0", "a1", "b0", "b1", "b2", "b0", "b1", "b2"), n)
  df$effort <- rep(c(rep("avec effort", 5), rep("sans effort", 3)), n)
}

ggplot(data = df, aes(y = 1:(n*(2+(nvar+1)*2)), x = med, xmin = low, xmax = up)) +
  facet_grid(effort ~ param, scales='free') + #
  geom_pointrange(position = position_dodge(width = .8), cex = 0.2) +
  geom_vline(aes(xintercept = vrai), col = "red") +
  labs(x = "", y = "", 
       title = "Comparaison des modèles avec et sans effort"
       # title = paste0("scenario ", num)
       ) +
  theme(
    axis.text.y = element_blank(),     # Enlève les textes de l'axe Y
    axis.ticks.y = element_blank(),    # Enlève les ticks de l'axe Y
    axis.title.y = element_blank(),    # Enlève le titre de l'axe Y
    panel.grid.major.y = element_blank(), # Enlève les lignes de grille horizontales majeures
    panel.grid.minor.y = element_blank(),  # Enlève les lignes de grille horizontales mineures
    panel.grid.minor.x = element_blank()
  )
      
# Complément de graphique 
# proportion d'IC contenant le vrai paramètre
prop_ok <- df %>%
  mutate(in_interval = (vrai >= low & vrai <= up)) %>% # crée une colonne logique
  group_by(param, effort) %>%                          # groupe par param et effort
  summarise(proportion = mean(in_interval), .groups = "drop") # moyenne = proportion

# largeur moyenne d'IC
largeur_IC <- df %>%
  mutate(IC = up-low) %>% 
  group_by(param, effort) %>%                     
  summarise(largeur = mean(IC), .groups = "drop") 


# df_plot <- left_join(df, prop_ok, by = c("param", "effort"))

ggplot(data = df, aes(y = 1:(n*(2+(nvar+1)*2)), x = med, xmin = low, xmax = up)) +
  facet_grid(effort ~ param,
             scales='free',
             space = 'free'
             ) + #
  geom_pointrange(position = position_dodge(width = .8), cex = 0.2) +
  geom_vline(aes(xintercept = vrai), col = "red") +
  labs(x = "", y = "", 
       title = "Estimation des coefficients, cas d'une très forte corrélation"
       # title = paste0("scenario ", num)
  ) +
  theme(
    axis.text.y = element_blank(),     # Enlève les textes de l'axe Y
    axis.ticks.y = element_blank(),    # Enlève les ticks de l'axe Y
    axis.title.y = element_blank(),    # Enlève le titre de l'axe Y
    panel.grid.major.y = element_blank(), # Enlève les lignes de grille horizontales majeures
    panel.grid.minor.y = element_blank(),  # Enlève les lignes de grille horizontales mineures
    panel.grid.minor.x = element_blank()
  ) + 
  geom_text(
    data = prop_ok,
    aes(x = -Inf, y = -Inf, label = paste0(round(proportion * 100, 1), "%")),
    hjust = -0.1, vjust = -0.2,
    inherit.aes = FALSE,
    col = "chartreuse4"
  ) + 
  geom_text(
    data = largeur_IC,
    aes(x = -Inf, y = -Inf, label = round(largeur,1)),
    hjust = -.1, vjust = -1.5,
    inherit.aes = FALSE,
    col =  "royalblue"
  )
        
