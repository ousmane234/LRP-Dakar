# LRP — Collecte de Déchets à Dakar

Application Shiny de résolution d'un **Location-Routing Problem (LRP)** pour optimiser la collecte des déchets ménagers à Dakar.

## Fonctionnalités

- **Chargement des données** : import CSV des ménages (id, longitude, latitude, poids_dechets)
- **Clustering K-means** : détection automatique du k optimal par la méthode du coude
- **Optimisation MILP** : modèle LRP résolu avec `ompr` + solveur GLPK/CBC
- **Résultats** : carte des tournées, KPIs, analyse de sensibilité, export CSV

## Structure

```
optimisation/
├── ui.R        # Interface utilisateur (shinydashboard)
├── server.R    # Logique serveur (clustering + optimisation)
└── README.md
```

## Dépendances R

```r
install.packages(c(
  "shiny", "shinydashboard", "leaflet", "DT", "ggplot2", "dplyr",
  "ompr", "ompr.roi", "ROI.plugin.glpk", "geosphere"
))
```

## Lancer l'application localement

```r
shiny::runApp("optimisation")
```

## Auteur

Ousmane LO — Semestre 2
