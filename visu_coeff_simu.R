# Pour visualiser les coefficients obtenus lors des 50 simulations d'un même IPP

load("simulation_scenario1.RData")

# on transforme en dataframe
res <- unlist(resultat)
df <- data.frame(
  a1_l = res[names(res) == "avec1"],
  a1_m = res[names(res) == "avec2"],
  a1_u = res[names(res) == "avec3"],
  a2_l = res[names(res) == "avec4"],
  a2_m = res[names(res) == "avec5"],
  a2_u = res[names(res) == "avec6"],
  b1_l = res[names(res) == "avec7"],
  b1_m = res[names(res) == "avec8"],
  b1_u = res[names(res) == "avec9"],
  b2_l = res[names(res) == "avec10"],
  b2_m = res[names(res) == "avec11"],
  b2_u = res[names(res) == "avec12"],
  b1_s_l = res[names(res) == "sans1"],
  b1_s_m = res[names(res) == "sans2"],
  b1_s_u = res[names(res) == "sans3"],
  b2_s_l = res[names(res) == "sans4"],
  b2_s_m = res[names(res) == "sans5"],
  b2_s_u = res[names(res) == "sans6"]
)



# plot 
p <- ggplot(data = df,
            aes(x = 1:50,
              # y = param,
              y = b1_s_m,
              ymin = b1_s_l,
              ymax = b1_s_u
            )) +
  # geom_hline(aes(yintercept = 0)) +
  geom_pointrange(position = position_dodge(width = .8))
  # labs(
  #   title = "Evolution des coefficients",
  #   subtitle = "Modèle à une détection, effort : données GBIF",
  #   x = "",
  #   y = "",
  #   color = "Année"
  # ) +
  # scale_y_discrete(
  #   labels = c(
  #     "intercept_eff",
  #     "densité GBIF",
  #     "intercept_int",
  #     "dist_eau",
  #     "logdensite",
  #     "agri",
  #     "preci_cum",
  #     "temp_min"
  #   )
  # )
p
