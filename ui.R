# ============================================================
#  Application Shiny — LRP Collecte de Déchets à Dakar
#  Auteur : Ousmane LO
#  Structure : app.R unique (ui + server)
# ============================================================

library(shiny)
library(shinydashboard)
library(leaflet)
library(DT)
library(ggplot2)
library(dplyr)

# -- Modules séparés (sourcer dans un vrai projet) -----------
# source("modules/clustering.R")
# source("modules/optimisation.R")
# source("modules/carte.R")

# ============================================================
#  UI
# ============================================================
ui <- dashboardPage(
    skin = "blue",
    
    dashboardHeader(
        title = "LRP — Collecte Dakar"
    ),
    
    dashboardSidebar(
        sidebarMenu(
            menuItem("Données",        tabName = "donnees",   icon = icon("database")),
            menuItem("Clustering",     tabName = "cluster",   icon = icon("object-group")),
            menuItem("Optimisation",   tabName = "optim",     icon = icon("route")),
            menuItem("Résultats",      tabName = "resultats", icon = icon("chart-bar"))
        )
    ),
    
    dashboardBody(
        tabItems(
            
            # ── ONGLET 1 : DONNÉES ──────────────────────────────
            tabItem(tabName = "donnees",
                    fluidRow(
                        box(title = "Chargement des données", width = 4, status = "primary",
                            fileInput("fichier_menages", "Fichier ménages (CSV/GeoJSON)",
                                      accept = c(".csv", ".geojson")),
                            helpText("Colonnes attendues : id, longitude, latitude, poids_dechets"),
                            hr(),
                            numericInput("cap_point",  "Capacité point de collecte Cⱼ (kg)", value = 500, min = 100),
                            numericInput("cap_camion", "Capacité camion Qᵥ (kg)",            value = 1000, min = 100),
                            numericInput("d_max",      "Distance max ménage→point Dmax (m)", value = 500,  min = 50),
                            numericInput("n_max",      "Budget : nb max points Nmax",        value = 10,   min = 1),
                            hr(),
                            # Coordonnées des nœuds spéciaux
                            h5("Dépôt de départ (o)"),
                            numericInput("depot_lon", "Longitude", value = -17.4441),
                            numericInput("depot_lat", "Latitude",  value = 14.6937),
                            h5("Site de décharge (s)"),
                            numericInput("decharge_lon", "Longitude", value = -17.3667),
                            numericInput("decharge_lat", "Latitude",  value = 14.7500),
                            hr(),
                            actionButton("charger", "Charger les données", class = "btn-primary btn-block")
                        ),
                        box(title = "Aperçu — Carte des ménages", width = 8, status = "primary",
                            leafletOutput("carte_menages", height = 450),
                            hr(),
                            verbatimTextOutput("resume_donnees")
                        )
                    )
            ),
            
            # ── ONGLET 2 : CLUSTERING ───────────────────────────
            tabItem(tabName = "cluster",
                    fluidRow(
                        box(title = "Paramètres K-means contraint", width = 4, status = "warning",
                            sliderInput("k_min", "k minimum", min = 2, max = 20, value = 5),
                            sliderInput("k_max", "k maximum", min = 5, max = 50, value = 20),
                            helpText("L'algorithme teste toutes les valeurs de k dans cet intervalle."), ##########" modifier le texte
                            hr(),
                            checkboxInput("contrainte_charge",
                                          "Appliquer contrainte de charge (K-means contraint)", value = TRUE),
                            checkboxInput("projection_osm",
                                          "Projeter centroïdes sur réseau routier (OSM)", value = TRUE),
                            hr(),
                            numericInput("n_init", "Nb d'initialisations (K-means++)", value = 10, min = 1),
                            actionButton("lancer_cluster", "Lancer le clustering", class = "btn-warning btn-block")
                        ),
                        box(title = "Résultats du clustering", width = 8, status = "warning",
                            tabsetPanel(
                                tabPanel("Carte des clusters",
                                         leafletOutput("carte_clusters", height = 380)
                                ),
                                tabPanel("Courbe d'inertie",
                                         plotOutput("courbe_inertie", height = 380),
                                         helpText("Méthode du coude : choisir k là où la courbe s'infléchit.")
                                ),
                                tabPanel("Points candidats J",
                                         DTOutput("table_candidats")
                                )
                            )
                        )
                    )
            ),
            
            # ── ONGLET 3 : OPTIMISATION ─────────────────────────
            tabItem(tabName = "optim",
                    fluidRow(
                        box(title = "Paramètres du modèle LRP", width = 4, status = "danger",
                            h5("Flotte de véhicules"),
                            numericInput("nb_vehicules", "Nombre de véhicules |V|", value = 3, min = 1),
                            hr(),
                            h5("Solveur"),
                            selectInput("solveur", "Solveur",
                                        choices = c("GLPK (gratuit)" = "glpk",
                                                    "CBC (gratuit)"  = "cbc",
                                                    "Gurobi"         = "gurobi")),
                            numericInput("time_limit", "Limite de temps (secondes)", value = 120, min = 10),
                            hr(),
                            h5("Fonction objectif"),
                            selectInput("obj", "Minimiser",
                                        choices = c("Distance totale" = "distance",
                                                    "Nombre de points ouverts" = "nb_points",
                                                    "Distance + pénalité ouverture" = "mixte")),
                            hr(),
                            actionButton("lancer_optim", "Lancer l'optimisation", class = "btn-danger btn-block"),
                            br(),
                            verbatimTextOutput("statut_optim")
                        ),
                        box(title = "Formulation du modèle", width = 8, status = "danger",
                            tabsetPanel(
                                tabPanel("Variables",
                                         helpText("Variables de décision du modèle LRP aller simple :"),
                                         tags$ul(
                                             tags$li("yⱼ ∈ {0,1} : point j ouvert"),
                                             tags$li("xᵢⱼ ∈ {0,1} : ménage i affecté au point j"),
                                             tags$li("zᵛᵢⱼ ∈ {0,1} : véhicule v emprunte l'arc i→j"),
                                             tags$li("uᵛⱼ ≥ 0 : charge cumulée du véhicule v en j")
                                         )
                                ),
                                tabPanel("Contraintes actives",
                                         verbatimTextOutput("contraintes_actives")
                                ),
                                tabPanel("Log solveur",
                                         verbatimTextOutput("log_solveur")
                                )
                            )
                        )
                    )
            ),
            
            # ── ONGLET 4 : RÉSULTATS ────────────────────────────
            tabItem(tabName = "resultats",
                    fluidRow(
                        # KPIs
                        valueBoxOutput("kpi_distance",  width = 3),
                        valueBoxOutput("kpi_points",    width = 3),
                        valueBoxOutput("kpi_vehicules", width = 3),
                        valueBoxOutput("kpi_couverture",width = 3)
                    ),
                    fluidRow(
                        box(title = "Carte des tournées", width = 8, status = "success",
                            leafletOutput("carte_tournees", height = 450),
                            helpText("Chaque couleur correspond à un véhicule. ▲ = points ouverts, ★ = dépôt, ■ = décharge.")
                        ),
                        box(title = "Détail des tournées", width = 4, status = "success",
                            DTOutput("table_tournees"),
                            hr(),
                            downloadButton("export_resultats", "Exporter les résultats (CSV)")
                        )
                    ),
                    fluidRow(
                        box(title = "Analyse de sensibilité", width = 12, status = "info",
                            fluidRow(
                                column(3,
                                       sliderInput("sens_nmax", "Faire varier Nmax", min = 1, max = 30, value = c(5,20))
                                ),
                                column(9,
                                       plotOutput("graphe_sensibilite", height = 280)
                                )
                            )
                        )
                    )
            )
        )
    )
)



