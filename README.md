# Projet Ragondin

### Description du contenu

* `Ragondins_2019.R` contient le code permettant de créer la table grid_sf contenant toutes les données nécessaires sur l'année 2019 : observations de ragondins et covariables (par unité d'aire).
Le nettoyage des cellules génantes de la frontière n'est pas fait dans ce code. `grid_sf` est disponible dans le fichier `grid_sf.RData`.

* `model_ragondins_2019.R` contient le code pour fitter le modèle associé à l'année 2019. Les deux approches de l'effort y sont présentes. 

* `buffer_eau.R` permet de construire le buffer autour des rivières et plans d'eau d'Occitanie. Les résultats obtenus sont disponibles dans `buffer.RData`, `buffer_plan.RData`, `buffer_riv.RData`.

* `gbif_RData` contient les données GBIF de l'année 2019 sur une zone englobant l'Occitanie
