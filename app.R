# ============================================================
#  app.R — LRP Collecte de Déchets à Dakar
#  Auteur : Ousmane LO
#
#  Point d'entrée unique. Lance l'application avec :
#    shiny::runApp(".")            # depuis le dossier du projet
#    shiny::runApp("dakar_lrp/")  # depuis le dossier parent
#
#  Structure du projet :
#    app.R
#    ui.R
#    server.R
#    modules/
#      clustering.R     K-means contraint + détection du coude
#      distances.R      Matrice étendue dépôt + candidats + décharge
#      modele_mip.R     Formulation et résolution LRP (ompr)
#    utils/
#      osm_routing.R    Routage routier réel via dodgr + OSM
#
#  Packages requis (installer une seule fois) :
#    install.packages(c(
#      "shiny", "shinydashboard", "leaflet", "DT",
#      "ggplot2", "dplyr", "geosphere",
#      "ompr", "ompr.roi", "ROI.plugin.glpk",
#      "dodgr", "osmdata", "sf"
#    ))
# ============================================================

# ── Packages ────────────────────────────────────────────────
library(shiny)
library(shinydashboard)
library(leaflet)
library(DT)
library(ggplot2)
library(dplyr)
library(geosphere)
library(ompr)
library(ompr.roi)
library(ROI.plugin.glpk)
library(dodgr)
library(osmdata)
library(sf)

# ── Chargement des modules métier ───────────────────────────
source("modules/distances.R")
source("modules/modele_mip.R")
source("modules/clustering.R")
source("utils/osm_routing.R")
source("modules/resultats.R")

# ── UI et Server ────────────────────────────────────────────
source("ui.R")
source("server.R")

# ── Lancement ───────────────────────────────────────────────
shinyApp(ui = ui, server = server)