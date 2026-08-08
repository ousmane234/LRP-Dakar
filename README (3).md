# Optimisation de la collecte des déchets à Dakar

Application R Shiny pour l'optimisation des tournées de collecte des déchets à Dakar, combinant modélisation en optimisation (MIP) et données géospatiales (réseau routier OpenStreetMap).

## Prérequis

* **R** (version 4.x recommandée) — [télécharger ici](https://cran.r-project.org/)
* **RStudio** (recommandé mais pas obligatoire) — [télécharger ici](https://posit.co/download/rstudio-desktop/)
* Une connexion internet (pour l'installation initiale des packages et le téléchargement des données OSM)

Ce projet utilise [`renv`](https://rstudio.github.io/renv/) pour figer les versions exactes des packages R utilisés. Vous n'avez **pas besoin d'installer les packages manuellement** : `renv` s'en charge automatiquement.

## Installation

1. **Cloner ou télécharger le projet**

```bash
   git clone <https://github.com/ousmane234/LRP-Dakar.git

>
   cd optimisation--collecte-des-dechets-a-dakar
   ```

(ou téléchargez et dézippez l'archive du projet)

2. **Ouvrir le projet dans RStudio**

   * Double-cliquez sur le fichier `.Rproj`, ou
   * Ouvrez RStudio puis `File > Open Project...` et sélectionnez le dossier du projet
3. **Laisser `renv` s'activer automatiquement**
À l'ouverture du projet, `renv` devrait s'activer tout seul (grâce au fichier `.Rprofile`). Vous verrez un message dans la console du type :

```
   Project 'nom-du-projet' loaded. \\\[renv x.x.x]
   ```

4. **Installer les packages nécessaires**
Dans la console R, lancez :

```r
   renv::restore()
   ```

Cette commande installe automatiquement **toutes les bonnes versions** des packages listés dans `renv.lock`. Cela peut prendre quelques minutes la première fois. Répondez `y` (yes) si une confirmation est demandée.

5. **Vérifier que tout est en ordre**

```r
   renv::status()
   ```

Le message `No issues found -- the project is in a consistent state.` confirme que l'installation s'est bien passée.

## Lancer l'application

Une fois les packages installés, lancez l'application Shiny :

```r
shiny::runApp()
```

Ou, si vous êtes dans RStudio, ouvrez `app.R` et cliquez sur le bouton **"Run App"** en haut de l'éditeur.

L'application devrait s'ouvrir dans une fenêtre RStudio ou dans votre navigateur par défaut.

## Structure du projet

```
.
├── app.R                          # Point d'entrée principal de l'application Shiny
├── data/
│   ├── buildings.geojson          # Bâtiments (données géospatiales)
│   ├── routes\\\_dakar.geojson       # Réseau routier de Dakar
│   ├── graphe\\\_dakar.Rds           # Graphe routier pré-calculé (utilisé pour le routage)
│   ├── menages\\\_gps.csv            # Jeu de données de test (voir note ci-dessous)
│   └── ...                        # Autres fichiers de données
├── modules/
│   ├── cache\\\_graphe.R             # Gestion du cache du graphe routier
│   ├── chargement\\\_matrice\\\_osm.R   # Chargement de la matrice de distances OSM
│   ├── clustering.R               # Clustering des points de collecte
│   ├── distances.R                # Calcul des distances/temps de trajet
│   ├── modele\\\_mip.R               # Modèle d'optimisation (programmation en nombres entiers mixtes)
│   └── resultats.R                # Module d'affichage et de traitement des résultats
├── scripts/
│   └── preparer\\\_graphe\\\_dakar.R    # Script de préparation/génération du graphe routier de Dakar
├── utils/
│   └── osm\\\_routing.R              # Fonctions utilitaires pour le routage via OpenStreetMap
├── rsconnect/                      # Config de déploiement shinyapps.io/Posit Connect (non versionné, contient des tokens)
├── renv/                           # Dossier de gestion des dépendances (géré par renv)
├── renv.lock                        # Liste figée des packages R et leurs versions exactes
├── manifest.json                    # Manifeste de déploiement (généré par rsconnect)
├── .Rprofile                        # Active renv automatiquement à l'ouverture du projet
├── AGENTS.md                         # \\\[À COMPLÉTER : précisez le contenu/rôle de ce fichier]
└── README.md                         # Ce fichier
```

> \\\*\\\*Note :\\\*\\\* les données brutes (` `\\\*.geojson`, `\\\*.xlsx`) ne sont pas versionnées dans le dépôt Git (exclues via `.gitignore`) pour des raisons de taille/confidentialité. Le dossier 

## Tester l'application

Pour tester rapidement l'application sans avoir à générer ou récupérer un jeu de données complet, utilisez le fichier **`data/menages\\\_gps.csv`**. Il contient un échantillon de données prêt à l'emploi, suffisant pour vérifier que l'application fonctionne correctement de bout en bout (chargement, calcul du modèle d'optimisation, affichage des résultats).

## Dépendances principales

Le projet s'appuie notamment sur :

* **Shiny** — framework de l'application web interactive
* **osmdata**, **dodgr**, **geodist** — récupération et analyse du réseau routier OpenStreetMap, calcul de distances/itinéraires
* **rvest**, **xml2**, **httr** — récupération de données web
* Un solveur d'optimisation pour le modèle MIP \[ `Rglpk`, `ompr`]

L'ensemble des packages et leurs versions exactes sont consultables dans `renv.lock`.

## Problèmes courants

* **`renv::restore()` échoue ou reste bloqué** : vérifiez votre connexion internet, puis relancez la commande.
* **Une erreur mentionne un package manquant malgré `renv::restore()`** : lancez `renv::status()` pour diagnostiquer, puis `renv::snapshot()` si un package a été ajouté après coup.
* **L'application ne se lance pas / erreur de chemin** : assurez-vous d'avoir bien ouvert le projet via le fichier `.Rproj` (et non un simple dossier), afin que le working directory soit correctement positionné à la racine du projet.

## Contact

tel : 783155011

email : lousmane2021@gmail.com

