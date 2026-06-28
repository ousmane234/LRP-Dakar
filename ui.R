# ============================================================
#  ui.R — LRP Collecte de Déchets à Dakar
#  Auteur : Ousmane LO
# ============================================================

library(shiny)
library(shinydashboard)
library(leaflet)
library(DT)

ui <- dashboardPage(
    skin = "blue",
    
    # ── En-tête ────────────────────────────────────────────────
    dashboardHeader(
        title = "LRP — Collecte Dakar"
    ),
    
    # ── Sidebar ────────────────────────────────────────────────
    dashboardSidebar(
        sidebarMenu(
            id = "menu_actif",
            menuItem("Données",      tabName = "donnees",   icon = icon("database")),
            menuItem("Clustering",   tabName = "cluster",   icon = icon("object-group")),
            menuItem("Optimisation", tabName = "optim",     icon = icon("route")),
            menuItem("Résultats",    tabName = "resultats", icon = icon("chart-bar"))
        ),
        hr(),
        # Indicateur de progression globale visible depuis toutes les pages
        div(style = "padding: 0 15px;",
            h5("Progression", style = "color:#aaa; font-size:12px; margin-bottom:4px;"),
            uiOutput("badge_etapes")
        )
    ),
    
    # ── Corps ──────────────────────────────────────────────────
    dashboardBody(
        
        # CSS léger pour les badges de progression et les points surchargés
        tags$head(tags$style(HTML("
      .badge-ok   { background:#1D9E75; color:#fff; border-radius:3px; padding:2px 7px; font-size:11px; }
      .badge-wait { background:#aaa;    color:#fff; border-radius:3px; padding:2px 7px; font-size:11px; }
      .surcharge  { color:#D85A30; font-weight:bold; }
      .leaflet-container { cursor: crosshair; }
    "))),
        
        tabItems(
            
            # ════════════════════════════════════════════════════════
            #  ONGLET 1 — DONNÉES
            # ════════════════════════════════════════════════════════
            tabItem(tabName = "donnees",
                    fluidRow(
                        
                        # Panneau de gauche : paramètres
                        box(title = "Chargement & paramètres", width = 4, status = "primary",
                            solidHeader = TRUE,
                            
                            fileInput("fichier_menages", "Fichier ménages (CSV)",
                                      accept = c(".csv"),
                                      placeholder = "id, longitude, latitude, poids_dechets"),
                            helpText("Colonnes obligatoires : id, longitude, latitude, poids_dechets"),
                            hr(),
                            
                            h5(icon("sliders"), " Paramètres capacitaires"),
                            numericInput("cap_point",  "Capacité point de collecte Cⱼ (kg)",
                                         value = 500, min = 100, step = 50),
                            numericInput("cap_camion", "Capacité camion Qᵥ (kg)",
                                         value = 1000, min = 100, step = 50),
                            numericInput("d_max",      "Distance max ménage → point Dmax (m)",
                                         value = 500, min = 50, step = 50),
                            numericInput("n_max",      "Budget : nb max de points Nmax",
                                         value = 10, min = 1),
                            hr(),
                            
                            h5(icon("map-marker"), " Nœuds spéciaux"),
                            fluidRow(
                                column(6, numericInput("depot_lon",    "Dépôt — Lon",  value = -17.4441, step = 0.0001)),
                                column(6, numericInput("depot_lat",    "Dépôt — Lat",  value =  14.6937, step = 0.0001))
                            ),
                            fluidRow(
                                column(6, numericInput("decharge_lon", "Décharge — Lon", value = -17.3667, step = 0.0001)),
                                column(6, numericInput("decharge_lat", "Décharge — Lat", value =  14.7500, step = 0.0001))
                            ),
                            helpText("Cliquez sur la carte pour repérer les coordonnées."),
                            hr(),
                            
                            actionButton("charger", "Charger les données",
                                         class = "btn-primary btn-block", icon = icon("upload"))
                        ),
                        
                        # Panneau de droite : carte + résumé
                        box(title = "Carte des ménages", width = 8, status = "primary",
                            solidHeader = TRUE,
                            leafletOutput("carte_menages", height = 430),
                            hr(),
                            verbatimTextOutput("resume_donnees")
                        )
                    )
            ),
            
            # ════════════════════════════════════════════════════════
            #  ONGLET 2 — CLUSTERING
            # ════════════════════════════════════════════════════════
            tabItem(tabName = "cluster",
                    fluidRow(
                        
                        box(title = "Paramètres K-means contraint", width = 4, status = "warning",
                            solidHeader = TRUE,
                            
                            sliderInput("k_min", "k minimum", min = 2,  max = 20, value = 3),
                            sliderInput("k_max", "k maximum", min = 5,  max = 50, value = 20),
                            helpText("L'algorithme teste toutes les valeurs de k dans cet intervalle,",
                                     "puis sélectionne le coude optimal par la méthode de la droite."),
                            hr(),
                            
                            numericInput("n_init", "Nb d'initialisations K-means++",
                                         value = 10, min = 1, max = 50),
                            hr(),
                            
                            checkboxInput("contrainte_charge",
                                          "Appliquer contrainte de charge (split clusters surchargés)",
                                          value = TRUE),
                            helpText("Si un cluster dépasse Cⱼ, il est automatiquement divisé en 2."),
                            hr(),
                            
                            # Note : projection_osm est lue aussi lors de l'optimisation
                            checkboxInput("projection_osm",
                                          "Utiliser le réseau routier OSM (distances & itinéraires réels)",
                                          value = TRUE),
                            helpText("Désactiver pour un calcul rapide (distances à vol d'oiseau)."),
                            hr(),
                            
                            actionButton("lancer_cluster", "Lancer le clustering",
                                         class = "btn-warning btn-block", icon = icon("play"))
                        ),
                        
                        box(title = "Résultats du clustering", width = 8, status = "warning",
                            solidHeader = TRUE,
                            tabsetPanel(
                                tabPanel("Carte des clusters",
                                         leafletOutput("carte_clusters", height = 380)
                                ),
                                tabPanel("Courbe d'inertie",
                                         plotOutput("courbe_inertie", height = 340),
                                         helpText("Le point le plus éloigné de la droite (coude) donne k optimal.")
                                ),
                                tabPanel("Points candidats J",
                                         DTOutput("table_candidats"),
                                         helpText("⚠ Surchargé = charge > Cⱼ. Augmenter Cⱼ ou Nmax si besoin.")
                                )
                            )
                        )
                    )
            ),
            
            # ════════════════════════════════════════════════════════
            #  ONGLET 3 — OPTIMISATION
            # ════════════════════════════════════════════════════════
            tabItem(tabName = "optim",
                    fluidRow(
                        
                        box(title = "Paramètres du modèle LRP", width = 4, status = "danger",
                            solidHeader = TRUE,
                            
                            h5(icon("truck"), " Flotte"),
                            numericInput("nb_vehicules", "Nombre de véhicules |V|",
                                         value = 3, min = 1, max = 20),
                            hr(),
                            
                            h5(icon("cogs"), " Solveur"),
                            selectInput("solveur", "Solveur MIP",
                                        choices = c("GLPK (gratuit)" = "glpk",
                                                    "CBC (gratuit)"  = "cbc",
                                                    "Gurobi"         = "gurobi"),
                                        selected = "glpk"),
                            numericInput("time_limit", "Limite de temps (secondes)",
                                         value = 120, min = 10, max = 3600),
                            hr(),
                            
                            h5(icon("bullseye"), " Fonction objectif"),
                            selectInput("obj", "Minimiser",
                                        choices = c("Distance totale (km)"        = "distance",
                                                    "Nombre de points ouverts"    = "nb_points",
                                                    "Distance + pénalité ouverture" = "mixte"),
                                        selected = "distance"),
                            hr(),
                            
                            actionButton("lancer_optim", "Lancer l'optimisation",
                                         class = "btn-danger btn-block", icon = icon("play")),
                            br(),
                            verbatimTextOutput("statut_optim")
                        ),
                        
                        box(title = "Détail du modèle", width = 8, status = "danger",
                            solidHeader = TRUE,
                            tabsetPanel(
                                tabPanel("Variables",
                                         tags$ul(style = "margin-top:10px; line-height:2;",
                                                 tags$li(tags$b("yⱼ ∈ {0,1}"),      " : point candidat j ouvert"),
                                                 tags$li(tags$b("xᵢⱼ ∈ {0,1}"),     " : ménage i affecté au point j"),
                                                 tags$li(tags$b("zᵛᵢⱼ ∈ {0,1}"),    " : véhicule v emprunte l'arc i → j"),
                                                 tags$li(tags$b("uᵛⱼ ≥ 0"),          " : charge cumulée du véhicule v en j (MTZ)")
                                         ),
                                         hr(),
                                         tags$p(tags$b("Graphe étendu :"),
                                                " nœud 1 = Dépôt (O), nœuds 2..(n+1) = candidats, nœud (n+2) = Décharge (S)")
                                ),
                                tabPanel("Contraintes",
                                         verbatimTextOutput("contraintes_actives")
                                ),
                                tabPanel("Log solveur",
                                         verbatimTextOutput("log_solveur")
                                )
                            )
                        )
                    )
            ),
            
            # ════════════════════════════════════════════════════════
            #  ONGLET 4 — RÉSULTATS
            # ════════════════════════════════════════════════════════
            tabItem(tabName = "resultats",
                    
                    # KPIs
                    fluidRow(
                        valueBoxOutput("kpi_distance",   width = 3),
                        valueBoxOutput("kpi_points",     width = 3),
                        valueBoxOutput("kpi_vehicules",  width = 3),
                        valueBoxOutput("kpi_couverture", width = 3)
                    ),
                    
                    # Carte + table des tournées
                    fluidRow(
                        box(title = "Carte des tournées", width = 8, status = "success",
                            solidHeader = TRUE,
                            leafletOutput("carte_tournees", height = 430),
                            helpText(
                                icon("circle", style = "color:#378ADD"), " Véhicule 1  ",
                                icon("circle", style = "color:#D85A30"), " Véhicule 2  ",
                                icon("circle", style = "color:#1D9E75"), " Véhicule 3  ",
                                " | ▲ = point de collecte  ★ = dépôt  ✖ = décharge"
                            )
                        ),
                        box(title = "Détail des tournées", width = 4, status = "success",
                            solidHeader = TRUE,
                            DTOutput("table_tournees"),
                            hr(),
                            downloadButton("export_resultats", "Exporter (CSV)",
                                           class = "btn-success btn-block")
                        )
                    ),
                    
                    # Distance par véhicule
                    fluidRow(
                        box(title = "Distance par véhicule", width = 6, status = "info",
                            solidHeader = TRUE,
                            plotOutput("graphe_vehicules", height = 250)
                        ),
                        
                        # Analyse de sensibilité Nmax
                        box(title = "Sensibilité au budget Nmax", width = 6, status = "info",
                            solidHeader = TRUE,
                            fluidRow(
                                column(4,
                                       sliderInput("sens_nmax", "Plage Nmax",
                                                   min = 1, max = 30, value = c(5, 20))
                                ),
                                column(8,
                                       plotOutput("graphe_sensibilite", height = 220)
                                )
                            )
                        )
                    )
            )
        )
    )
)