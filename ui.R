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
  
  dashboardHeader(title = "LRP — Collecte Dakar"),
  
  dashboardSidebar(
    sidebarMenu(
      id = "menu_actif",
      menuItem("Données",      tabName = "donnees",   icon = icon("database")),
      menuItem("Clustering",   tabName = "cluster",   icon = icon("object-group")),
      menuItem("Optimisation", tabName = "optim",     icon = icon("route")),
      menuItem("Résultats",    tabName = "resultats", icon = icon("chart-bar"))
    ),
    hr(),
    div(style = "padding: 0 15px;",
        h5("Progression", style = "color:#aaa; font-size:12px; margin-bottom:4px;"),
        uiOutput("badge_etapes")
    )
  ),
  
  dashboardBody(
    
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
                    actionButton("charger", "Charger les données",
                                 class = "btn-primary btn-block", icon = icon("upload"))
                ),
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
                    helpText("L'algorithme teste toutes les valeurs de k dans cet intervalle."),
                    hr(),
                    numericInput("n_init", "Nb d'initialisations K-means++",
                                 value = 10, min = 1, max = 50),
                    hr(),
                    checkboxInput("contrainte_charge",
                                  "Appliquer contrainte de charge (split clusters surchargés)",
                                  value = TRUE),
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
                               plotOutput("courbe_inertie", height = 340)
                      ),
                      tabPanel("Points candidats J",
                               DTOutput("table_candidats")
                      )
                    )
                )
              )
      ),
      
      # ════════════════════════════════════════════════════════
      #  ONGLET 3 — OPTIMISATION (MODIFIÉ)
      # ════════════════════════════════════════════════════════
      tabItem(tabName = "optim",
              fluidRow(
                
                box(title = "Paramètres du modèle LRP", width = 4, status = "danger",
                    solidHeader = TRUE,
                    
                    # ── Coordonnées utilisateur ──
                    h5(icon("map-marker"), " Coordonnées utilisateur"),
                    fluidRow(
                      column(6, numericInput("depot_lon_opt", "Dépôt — Lon", 
                                             value = -17.4441, step = 0.0001)),
                      column(6, numericInput("depot_lat_opt", "Dépôt — Lat", 
                                             value = 14.6937, step = 0.0001))
                    ),
                    fluidRow(
                      column(6, numericInput("decharge_lon_opt", "Décharge — Lon", 
                                             value = -17.3667, step = 0.0001)),
                      column(6, numericInput("decharge_lat_opt", "Décharge — Lat", 
                                             value = 14.7500, step = 0.0001))
                    ),
                    helpText("Ces coordonnées sont utilisées pour le calcul des distances sur le réseau OSM."),
                    hr(),
                    
                    # ── Calcul de la matrice ──
                    actionButton("calculer_matrice", 
                                 "Calculer la matrice OSM",
                                 class = "btn-primary btn-block", 
                                 icon = icon("route")),
                    br(),
                    verbatimTextOutput("statut_matrice"),
                    hr(),
                    
                    # ── Paramètres du modèle ──
                    h5(icon("truck"), " Flotte"),
                    numericInput("nb_vehicules", "Nombre de véhicules |V|",
                                 value = 3, min = 1, max = 20),
                    hr(),
                    
                    h5(icon("cogs"), " Solveur"),
                    selectInput("solveur", "Solveur MIP",
                                choices = c("GLPK (gratuit)" = "glpk",
                                            "CBC (gratuit)"  = "cbc"),
                                selected = "glpk"),
                    numericInput("time_limit", "Limite de temps (secondes)",
                                 value = 120, min = 10, max = 3600),
                    hr(),
                    
                    actionButton("lancer_optim", 
                                 "Lancer l'optimisation",
                                 class = "btn-danger btn-block", 
                                 icon = icon("play")),
                    br(),
                    verbatimTextOutput("statut_optim")
                ),
                
                box(title = "Détail du modèle", width = 8, status = "danger",
                    solidHeader = TRUE,
                    tabsetPanel(
                      tabPanel("Contraintes",
                               verbatimTextOutput("contraintes_actives")
                      ),
                      tabPanel("Log solveur",
                               verbatimTextOutput("log_solveur")
                      ),
                      tabPanel("Aperçu matrice",
                               DTOutput("apercu_matrice")
                      )
                    )
                )
              )
      ),
      
      # ════════════════════════════════════════════════════════
      #  ONGLET 4 — RÉSULTATS
      # ════════════════════════════════════════════════════════
      tabItem(tabName = "resultats",
              fluidRow(
                valueBoxOutput("kpi_distance",   width = 3),
                valueBoxOutput("kpi_points",     width = 3),
                valueBoxOutput("kpi_vehicules",  width = 3),
                valueBoxOutput("kpi_couverture", width = 3)
              ),
              fluidRow(
                box(title = "Carte des tournées", width = 8, status = "success",
                    solidHeader = TRUE,
                    leafletOutput("carte_tournees", height = 430)
                ),
                box(title = "Détail des tournées", width = 4, status = "success",
                    solidHeader = TRUE,
                    DTOutput("table_tournees"),
                    hr(),
                    downloadButton("export_resultats", "Exporter (CSV)",
                                   class = "btn-success btn-block")
                )
              ),
              fluidRow(
                box(title = "Distance par véhicule", width = 6, status = "info",
                    solidHeader = TRUE,
                    plotOutput("graphe_vehicules", height = 250)
                ),
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