# Projet Ragondin

## Description du contenu

### Les scripts

-   **`Donnees_brutes.R`** contient le code permettant d'extraire les données nécessaires et de les ramener à l'échelle de l'Occitanie (avec ou sans buffer selon le cas). Les variables construites sont enregistrées dans le dossier `Data/A_charger`.

-   **`Construction_grille.R`** contient le code permettant de construire la table `grid_sf` donnant la grille à l'échelle voulue et contenant toutes les données associées sur la période 2010-2024 : observations de ragondins et covariables.

-   **`model_periode.R`** contient le code permettant de fitter les différents modèles associés à la période voulue : en fonction de l'approche de l'effort retenue et du nombre d'observations considéré par cellule.

-   **`mapview.R`** contient les essais d'utilisation du package mapview.

-   **`visu.R`** contient le code nécessaire pour les vérifications et études visuelles des modèles obtenus.

-   **`CLC_evolution12-18.R`** et **`Script_CLC.R`** contiennent le code permettant de mettre en évidence les changements minimes de couverture agricole entre les données de 2012 et de 2018 de la base Corine Land Cover.

-   **`Script_pop.R`** contient le code permettant de mettre en évidence les changements de densité de population entre les données eurostat de 2011, 2018, 2021 et leurs influences sur les modèles.

-   **`Distrib_IEP.R`** contient le code pour transformer les simulations MCMC des coefficients du modèle en simulations de l'intensité, de l'effort et de la probabilité de présence et d'en déduire les estimateurs associés ainsi que leurs écart-types. Ceci permet d'obtenir les cartes associées.

-   **`model_periode_avec_IEP.R`** est une modification de **`model_periode.R`** permettant de garder les simulations de l'intensité, de l'effort et de la probabilité de présence, avec un thinning élevé (pour alléger la mémoire nécessaire). Donne des résultats proches de **`Distrib_IEP.R`** en conservant moins de simulations.

-   **`Ajout_dat_poly.R`** permet de comparer les modèles obtenus en utilisant les données ponctuelles seules ou les données ponctuelles et les centroïdes des polygones. La décision a été prise de continuer à travailler avec les données ponctuelles seules.

-   **`Validation.R`** contient les éléments de validation du modèle de processus ponctuel de Poisson.

-   Le dossier *Anciens_codes* contient les scripts de travail portant sur les différents éléments.

### Les données

Le dossier `Data` contient les bases de données utilisées, ainsi que le dossier `A_charger` qui rassemble les sorties de **`Donnees_brutes.R`** résumant les données utiles à la construction de la grille.

### RData

(à mettre à jour)

Ce dossier contient les fichiers .RData permettant de charger des variables de l'environnement R sans avoir à tout recompiler.

-   `buffer...` : fournit les buffers autour des plans d'eau et des rivières

-   `gbif_data` : les données GBIF importées qui sont utilisées pour créer la variable d'effort GBIF

-   `grid_sf...` : les différentes grilles construites en fonction de la taille de maille

-   `out...` : les résultats obtenus après l'estimation des différents modèles.

    -   env / gbif : si l'effort est décrit à partir des routes et chemins ou à partir des données GBIF

    -   mult : si toutes les détections par cellule sont utilisées

    -   dist : si la présence d'eau est décrite à partir des distances aux plans d'eau et aux rivières

### Images

(à mettre à jour)

Ce dossier contient les cartes d'intensité, d'effort et de probabilité de présence obtenues à partir des différents modèles.
