# ============================================================
#  app.R — LRP Collecte de Déchets à Dakar
#  Auteur : Ousmane LO
#
#  Point d'entrée unique. Lance l'application avec :
#    shiny::runApp(".")
# ============================================================

# ── Packages ────────────────────────────────────────────────
library(shiny)
library(shinydashboard)
library(leaflet)
library(DT)
library(ggplot2)
library(dplyr)
library(geosphere)
library(digest)
library(ompr)
library(ompr.roi)
library(ROI.plugin.glpk)
library(dodgr)
library(osmdata)
library(sf)

# ── Chargement des modules ─────────────────────────────────
source("modules/clustering.R")
source("modules/distances.R")
source("utils/osm_routing.R")
source("modules/cache_graphe.R")      # NOUVEAU
source("modules/modele_mip.R")
source("modules/resultats.R")

# ── UI et Server ────────────────────────────────────────────
source("ui.R")
source("server.R")

# ── Lancement ───────────────────────────────────────────────
shinyApp(ui = ui, server = server)