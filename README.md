# Projet Ragondin

## Description du contenu

### Les scripts

-   **`Donnees_brutes.R`** contient le code permettant d'extraire les données nécessaires et de les ramener à l'échelle de l'Occitanie (avec ou sans buffer selon le cas). Les variables construites sont enregistrées dans le dossier `Data/A_charger`.

-   **`Construction_grille.R`** contient le code permettant de construire la table `grid_sf` donnant la grille à l'échelle voulue et contenant toutes les données associées sur la période 2010-2024 : observations de ragondins et covariables.

-   **`model_periode.R`** contient le code permettant de fitter les différents modèles associés à la période voulue (année par année), en fonction de l'approche de l'effort retenue et du nombre d'observations considéré par cellule.

-   **`mapview.R`** contient les essais d'utilisation du package mapview.

-   **`visu.R`** contient le code nécessaire pour les vérifications et études visuelles des modèles obtenus.

-   **`CLC_evolution12-18.R`** et **`Script_CLC.R`** contiennent le code permettant de mettre en évidence les changements minimes de couverture agricole entre les données de 2012 et de 2018 de la base Corine Land Cover.

-   **`Script_pop.R`** contient le code permettant de mettre en évidence les changements de densité de population entre les données eurostat de 2011, 2018, 2021 et leurs influences sur les modèles.

-   **`Distrib_IEP.R`** contient le code pour transformer les simulations MCMC des coefficients du modèle en simulations de l'intensité, de l'effort et de la probabilité de présence et d'en déduire les estimateurs associés ainsi que leurs écart-types. Ceci permet d'obtenir les cartes associées.

-   **`model_periode_avec_IEP.R`** est une modification de **`model_periode.R`** permettant de garder les simulations de l'intensité, de l'effort et de la probabilité de présence, avec un thinning élevé (pour alléger la mémoire nécessaire). Donne des résultats proches de **`Distrib_IEP.R`** en conservant moins de simulations.

-   **`Ajout_dat_poly.R`** permet de comparer les modèles obtenus en utilisant les données ponctuelles seules ou les données ponctuelles et les centroïdes des polygones. La décision a été prise de continuer à travailler avec les données ponctuelles seules.

-   **`Validation.R`** contient les éléments de validation du modèle de processus ponctuel de Poisson.

-   Le dossier *Anciens_codes* contient les scripts de travail portant sur les différents éléments.

-   **`code_multi_annees.R`**, **`inits_multi_annees.R`** et **`estim_multi_annees.R`** regroupent les éléments de codes utiles pour les modèles prenant en compte toutes les années : les codes, la recherche d'initialisation adéquate et enfin l'ajustement des modèles.

-   **`simu_scenario....R`** permettent les simulations d'un IPPP selon différents scénarios (et à partir de variables construites différemment). **`simdatadk_script.R`** reprend les fonctions utiles du package AHMbook. **`Simulation.R`** et **`Simulation_2cov.R`** regroupent les différents essais effectués (avec une ou deux covariables d'intensité).

-   **`visu_coeff_simu.R`** permet de visualiser les coefficients fittés et les intervalles de crédibilité associés pour les différents scénarios de simulation.

### Les données

Le dossier `Data` contient les bases de données utilisées, ainsi que le dossier `A_charger` qui rassemble les sorties de **`Donnees_brutes.R`** résumant les données utiles à la construction de la grille.

### RData

(à mettre à jour)

Ce dossier contient les fichiers .RData permettant de charger des variables de l'environnement R sans avoir à tout recompiler.

-   Les dossiers **`5km2`** et **`50km2`** contiennent les grilles construites et les variables associées.

-   Le dossier **`Inits`** contient les initiations des coefficients pour les différents gros modèles (afin d'éviter les logprob infinis)

-   Le dossier **`Simulations`** contient les résultats des différentes simulations effectuées

-   `buffer...` : fournit les buffers autour des plans d'eau et des rivières

-   `gbif_data` : les données GBIF importées qui sont utilisées pour créer la variable d'effort GBIF

### Images

(à mettre à jour)

Ce dossier contient les cartes d'intensité, d'effort et de probabilité de présence obtenues à partir des différents modèles.
