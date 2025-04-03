# Projet Ragondin

## Description du contenu

### Les scripts

-   **`Ragondins_2019.R`** contient le code permettant de créer la table `grid_sf` contenant toutes les données nécessaires sur l'année 2019 : observations de ragondins et covariables (par unité d'aire). Le nettoyage des cellules génantes de la frontière n'est pas fait dans ce code.

-   **`model_ragondins_2019.R`** contient le code pour fitter les différents modèles associés à l'année 2019 : en fonction de l'approche de l'effort retenue et du nombre d'observations considéré par cellule.

-   **`buffer_eau.R`** permet de construire le buffer autour des rivières et plans d'eau d'Occitanie.

-   **`mapview.R`** contient les essais d'utilisation du package mapview.

-   Le dossier *Anciens_codes* contient les scripts de travail portant sur les différents éléments.

### RData

Ce dossier contient les fichiers .RData permettant de charger des variables de l'environnement R sans avoir à tout recompiler.

-   `buffer...` : fournit les buffers autour des plans d'eau et des rivières

-   `gbif_data` : les données GBIF importées qui sont utilisées pour créer la variable d'effort GBIF

-   `grid_sf...` : les différentes grilles construites en fonction de la taille de maille

-   `out...` : les résultats obtenus après l'estimation des différents modèles.

    -   env / gbif : si l'effort est décrit à partir des routes et chemins ou à partir des données GBIF

    -   mult : si toutes les détections par cellule sont utilisées

    -   dist : si la présence d'eau est décrite à partir des distances aux plans d'eau et aux rivières

### Images

Ce dossier contient les cartes d'intensité, d'effort et de probabilité de présence obtenues à partir des différents modèles.
