# ============================================================
#  server.R — LRP Collecte de Déchets à Dakar
#  Auteur : Ousmane LO
#  Modifications : intégration modules distances.R et modele_mip.R
# ============================================================

library(shiny)
library(shinydashboard)
library(leaflet)
library(DT)
library(ggplot2)
library(dplyr)
library(geosphere)

# ── Chargement des modules ───────────────────────────────────
source("modules/distances.R")
source("modules/modele_mip.R")
source("modules/clustering.R")
source("utils/osm_routing.R")
source("modules/resultats.R")

server <- function(input, output, session) {
    
    # ── Données réactives ──────────────────────────────────────
    rv <- reactiveValues(
        menages       = NULL,   # data.frame : id, longitude, latitude, poids_dechets
        candidats     = NULL,   # data.frame : id, longitude, latitude, cluster, charge
        inerties      = NULL,   # data.frame : k, inertie
        res_dist      = NULL,   # liste retournée par construire_matrice_distances()
        solution      = NULL,   # liste retournée par resoudre_lrp()
        log_optim     = ""
    )
    
    # ══════════════════════════════════════════════════════════
    #  ONGLET 1 — DONNÉES
    # ══════════════════════════════════════════════════════════
    
    observeEvent(input$charger, {
        req(input$fichier_menages)
        df <- read.csv(input$fichier_menages$datapath)
        names(df) <- tolower(names(df))
        
        # Vérification des colonnes requises
        cols_req <- c("id", "longitude", "latitude", "poids_dechets")
        manquantes <- setdiff(cols_req, names(df))
        if (length(manquantes) > 0) {
            showNotification(
                paste("Colonnes manquantes :", paste(manquantes, collapse = ", ")),
                type = "error"
            )
            return()
        }
        rv$menages <- df
        showNotification(paste(nrow(df), "ménages chargés."), type = "message")
    })
    
    output$carte_menages <- renderLeaflet({
        req(rv$menages)
        leaflet(rv$menages) %>%
            addTiles() %>%
            addCircleMarkers(~longitude, ~latitude,
                             radius = 4, color = "#378ADD",
                             popup  = ~paste("Ménage", id, "—", poids_dechets, "kg")) %>%
            addMarkers(lng = input$depot_lon,    lat = input$depot_lat,
                       popup = "Dépôt (O)",
                       icon  = awesomeIcons(icon = "home",  markerColor = "green")) %>%
            addMarkers(lng = input$decharge_lon, lat = input$decharge_lat,
                       popup = "Décharge (S)",
                       icon  = awesomeIcons(icon = "trash", markerColor = "red"))
    })
    
    output$resume_donnees <- renderPrint({
        req(rv$menages)
        cat("Nb ménages        :", nrow(rv$menages), "\n")
        cat("Poids total (kg)  :", sum(rv$menages$poids_dechets), "\n")
        cat("Emprise lon       :", round(range(rv$menages$longitude), 4), "\n")
        cat("Emprise lat       :", round(range(rv$menages$latitude),  4), "\n")
    })
    
    # ══════════════════════════════════════════════════════════
    #  ONGLET 2 — CLUSTERING
    # ══════════════════════════════════════════════════════════
    
    observeEvent(input$lancer_cluster, {
        req(rv$menages)
        
        withProgress(message = "Clustering en cours...", value = 0, {
            
            incProgress(0.2, message = "Calcul des inerties K-means++...")
            
            res_cl <- clustering_contraint(
                menages              = rv$menages,
                k_min                = input$k_min,
                k_max                = input$k_max,
                n_init               = input$n_init,
                cap_point            = input$cap_point,
                appliquer_contrainte = input$contrainte_charge
            )
            
            rv$candidats <- res_cl$candidats
            rv$menages   <- res_cl$menages
            rv$inerties  <- res_cl$inerties
            
            incProgress(1, message = "Terminé !")
            
            if (length(res_cl$clusters_surcharges) > 0) {
                showNotification(
                    paste0("\u26a0 ", length(res_cl$clusters_surcharges),
                           " point(s) encore surchargé(s) après split. Augmenter cap_point ou Nmax."),
                    type = "warning", duration = 8
                )
            }
            
            showNotification(
                paste0("Clustering terminé — k optimal = ", res_cl$k_optimal,
                       ", k final = ", res_cl$k_final),
                type = "message"
            )
        })
    })
    
    output$carte_clusters <- renderLeaflet({
        req(rv$candidats, rv$menages)
        pal <- colorFactor(palette = "Set1", domain = rv$menages$cluster)
        leaflet() %>%
            addTiles() %>%
            addCircleMarkers(data = rv$menages,
                             ~longitude, ~latitude,
                             color = ~pal(cluster), radius = 4, opacity = 0.7,
                             popup = ~paste("Ménage", id, "— cluster", cluster)) %>%
            addMarkers(data = rv$candidats,
                       ~longitude, ~latitude,
                       popup = ~paste("Candidat", id, "— charge:", round(charge), "kg"))
    })
    
    output$courbe_inertie <- renderPlot({
        req(rv$inerties)
        ggplot(rv$inerties, aes(x = k, y = inertie)) +
            geom_line(color = "#378ADD", linewidth = 1) +
            geom_point(color = "#D85A30", size = 3) +
            labs(title = "Courbe d'inertie (méthode du coude)",
                 x = "Nombre de clusters k", y = "Inertie totale") +
            theme_minimal(base_size = 13)
    })
    
    output$table_candidats <- renderDT({
        req(rv$candidats)
        datatable(rv$candidats,
                  options  = list(pageLength = 10),
                  colnames = c("ID", "Longitude", "Latitude", "Cluster", "Charge (kg)"))
    })
    
    # ══════════════════════════════════════════════════════════
    #  ONGLET 3 — OPTIMISATION
    # ══════════════════════════════════════════════════════════
    
    observeEvent(input$lancer_optim, {
        req(rv$candidats, rv$menages)
        
        withProgress(message = "Optimisation en cours...", value = 0, {
            
            # ── Étape 1 : construire la matrice de distances ──────────
            incProgress(0.1, message = "Construction de la matrice de distances...")
            
            if (isTRUE(input$projection_osm)) {
                # Distances routières réelles via OSM/dodgr
                rv$res_dist <- tryCatch(
                    matrice_distances_osm(
                        candidats    = rv$candidats,
                        depot_lon    = input$depot_lon,
                        depot_lat    = input$depot_lat,
                        decharge_lon = input$decharge_lon,
                        decharge_lat = input$decharge_lat
                    ),
                    error = function(e) {
                        showNotification(
                            paste0("OSM indisponible, fallback Haversine : ", e$message),
                            type = "warning", duration = 6
                        )
                        construire_matrice_distances(
                            candidats    = rv$candidats,
                            depot_lon    = input$depot_lon,
                            depot_lat    = input$depot_lat,
                            decharge_lon = input$decharge_lon,
                            decharge_lat = input$decharge_lat
                        )
                    }
                )
            } else {
                # Distances Haversine (ligne droite) — mode rapide
                rv$res_dist <- construire_matrice_distances(
                    candidats    = rv$candidats,
                    depot_lon    = input$depot_lon,
                    depot_lat    = input$depot_lat,
                    decharge_lon = input$decharge_lon,
                    decharge_lat = input$decharge_lat
                )
            }
            
            # ── Étape 2 : résoudre le modèle LRP ──────────────────
            incProgress(0.2, message = "Lancement du solveur MIP...")
            
            sol <- tryCatch(
                resoudre_lrp(
                    res_dist     = rv$res_dist,
                    menages      = rv$menages,
                    nb_vehicules = input$nb_vehicules,
                    cap_camion   = input$cap_camion,
                    cap_point    = input$cap_point,
                    d_max        = input$d_max,
                    n_max        = input$n_max,
                    solveur      = input$solveur,
                    time_limit   = input$time_limit
                ),
                error = function(e) {
                    showNotification(paste("Erreur solveur :", e$message), type = "error")
                    NULL
                }
            )
            
            if (!is.null(sol)) {
                rv$solution  <- sol
                rv$log_optim <- paste0(
                    "Solveur     : ", input$solveur, "\n",
                    "Statut      : ", sol$statut, "\n",
                    "Objectif    : ", sol$cout_total, " km\n",
                    "Points ouverts : ", sol$nb_points, "\n",
                    "Couverture  : ", sol$couverture, " %\n\n",
                    "Distance par véhicule :\n",
                    paste(names(sol$distances_par_vehicule),
                          round(sol$distances_par_vehicule, 2), "km",
                          collapse = "\n")
                )
            }
            
            incProgress(1, message = "Terminé !")
        })
    })
    
    output$statut_optim <- renderPrint({
        if (is.null(rv$solution)) cat("En attente du lancement...")
        else cat("✓", rv$solution$statut, "\nCoût total :", rv$solution$cout_total, "km")
    })
    
    output$log_solveur <- renderPrint({ cat(rv$log_optim) })
    
    output$contraintes_actives <- renderPrint({
        cat("C1  : Σⱼ z^v_Oj = 1             ∀v  (départ dépôt)\n")
        cat("C2  : Σᵢ z^v_iS = 1             ∀v  (arrivée décharge)\n")
        cat("C3  : Σᵢ z^v_ij = Σₖ z^v_jk    ∀v,j (conservation flux)\n")
        cat("C4  : Σᵥ Σᵢ z^v_ij = yⱼ        ∀j  (visite si ouvert)\n")
        cat("C5  : z^v_ii = 0                ∀v,i (pas d'auto-boucle)\n")
        cat("C6  : z^v_Sj = 0                ∀v,j (pas d'arc depuis décharge)\n")
        cat("C7  : z^v_iO = 0                ∀v,i (pas d'arc vers dépôt)\n")
        cat("C8  : Σⱼ xᵢⱼ = 1               ∀i  (affectation unique)\n")
        cat("C9  : xᵢⱼ ≤ yⱼ                 ∀i,j (seulement si ouvert)\n")
        cat("C10 : xᵢⱼ = 0 si dᵢⱼ > Dmax   ∀i,j (accessibilité)\n")
        cat("C11 : Σᵢ wᵢ·xᵢⱼ ≤ Cⱼ·yⱼ       ∀j  (capacité point)\n")
        cat("C12 : Σⱼ yⱼ ≤ Nmax                  (budget)\n")
        cat("C13 : MTZ  uⱼ ≥ uᵢ + Cⱼ·zᵢⱼ - Qᵥ(1-zᵢⱼ)  (sous-tours)\n")
        cat("C14 : u^v_j ≤ Qᵥ               ∀v,j (borne capacité)\n")
        cat("C15 : u^v_O = 0                ∀v  (dépôt vide au départ)\n")
    })
    
    # ══════════════════════════════════════════════════════════
    #  ONGLET 4 — RÉSULTATS
    # ══════════════════════════════════════════════════════════
    
    output$kpi_distance <- renderValueBox({
        valueBox(
            value    = if (!is.null(rv$solution)) paste(rv$solution$cout_total, "km") else "—",
            subtitle = "Distance totale",
            icon     = icon("road"), color = "blue"
        )
    })
    
    output$kpi_points <- renderValueBox({
        valueBox(
            value    = if (!is.null(rv$solution)) rv$solution$nb_points else "—",
            subtitle = "Points ouverts",
            icon     = icon("map-marker"), color = "yellow"
        )
    })
    
    output$kpi_vehicules <- renderValueBox({
        valueBox(
            value    = if (!is.null(rv$solution)) rv$solution$nb_vehicules else "—",
            subtitle = "Véhicules utilisés",
            icon     = icon("truck"), color = "green"
        )
    })
    
    output$kpi_couverture <- renderValueBox({
        valueBox(
            value    = if (!is.null(rv$solution)) paste(rv$solution$couverture, "%") else "—",
            subtitle = "Ménages couverts",
            icon     = icon("home"), color = "purple"
        )
    })
    
    # ── Carte des tournées ──────────────────────────────────────
    output$carte_tournees <- renderLeaflet({
        req(rv$solution, rv$candidats, rv$res_dist)
        
        COLORS  <- c("#378ADD", "#D85A30", "#1D9E75", "#7F77DD", "#BA7517",
                     "#E41A1C", "#FF7F00", "#4DAF4A", "#984EA3", "#A65628")
        arcs    <- rv$solution$sol_z          # arcs actifs (z = 1)
        coords  <- rv$res_dist$coords_etendues  # tous les nœuds étendus
        
        m <- leaflet() %>% addTiles()
        
        # -- Tracé des arcs par véhicule
        if (nrow(arcs) > 0) {
            
            use_osm <- !is.null(rv$res_dist$graphe)   # TRUE si OSM a été utilisé
            
            if (use_osm) {
                # Itinéraires routiers réels via dodgr
                geos <- tryCatch(
                    geometries_tournees_osm(arcs, rv$res_dist),
                    error = function(e) NULL
                )
                if (!is.null(geos)) {
                    m <- tracer_tournees_osm(m, geos, COLORS)
                } else {
                    use_osm <- FALSE   # fallback lignes droites
                }
            }
            
            if (!use_osm) {
                # Fallback : lignes droites (Haversine ou OSM échoué)
                for (v_id in sort(unique(arcs$v))) {
                    arcs_v  <- arcs[arcs$v == v_id, ]
                    couleur <- COLORS[(v_id - 1) %% length(COLORS) + 1]
                    for (k in seq_len(nrow(arcs_v))) {
                        noeud_i  <- coords[coords$index == arcs_v$i[k], ]
                        noeud_j  <- coords[coords$index == arcs_v$j[k], ]
                        dist_arc <- round(rv$res_dist$d_matrix[arcs_v$i[k], arcs_v$j[k]], 2)
                        m <- m %>% addPolylines(
                            lng = c(noeud_i$longitude, noeud_j$longitude),
                            lat = c(noeud_i$latitude,  noeud_j$latitude),
                            color = couleur, weight = 3, opacity = 0.85,
                            label = paste0("Véhicule ", v_id, " : ",
                                           noeud_i$label, " → ", noeud_j$label,
                                           " (", dist_arc, " km)")
                        )
                    }
                }
            }
        }
        
        # -- Marqueurs des points candidats ouverts
        points_ouverts_df <- rv$candidats[rv$candidats$id %in% rv$solution$points_ouverts, ]
        if (nrow(points_ouverts_df) > 0) {
            m <- m %>% addCircleMarkers(
                data    = points_ouverts_df,
                ~longitude, ~latitude,
                radius  = 8, color = "#333333", fillColor = "#FFC107",
                fillOpacity = 0.9, weight = 2,
                popup   = ~paste0("Point ", id, "<br>Charge : ", round(charge), " kg")
            )
        }
        
        # -- Dépôt et décharge
        m <- m %>%
            addAwesomeMarkers(
                lng  = input$depot_lon, lat = input$depot_lat,
                icon = awesomeIcons(icon = "home",  markerColor = "green", library = "fa"),
                popup = "Dépôt (O)"
            ) %>%
            addAwesomeMarkers(
                lng  = input$decharge_lon, lat = input$decharge_lat,
                icon = awesomeIcons(icon = "trash", markerColor = "red",   library = "fa"),
                popup = "Décharge (S)"
            ) %>%
            # Légende des véhicules
            addLegend(
                position = "bottomright",
                colors   = COLORS[seq_len(rv$solution$nb_vehicules)],
                labels   = paste("Véhicule", seq_len(rv$solution$nb_vehicules),
                                 "-", round(rv$solution$distances_par_vehicule, 1), "km"),
                title    = "Tournées"
            )
        m
    })
    
    # ── Table des tournées ──────────────────────────────────────
    output$table_tournees <- renderDT({
        req(rv$solution)
        datatable(
            rv$solution$tournees_detail,
            options  = list(pageLength = 15, order = list(list(0, "asc"))),
            colnames = c("Véhicule", "De", "Vers", "Distance (km)")
        )
    })
    
    # ── Export CSV ──────────────────────────────────────────────
    output$export_resultats <- downloadHandler(
        filename = function() paste0("resultats_LRP_Dakar_", Sys.Date(), ".csv"),
        content  = function(file) {
            req(rv$solution)
            write.csv(rv$solution$tournees_detail, file, row.names = FALSE)
        }
    )
    
    # ── Analyse de sensibilité ──────────────────────────────────
    output$graphe_sensibilite <- renderPlot({
        req(rv$candidats)
        nmax_vals <- seq(input$sens_nmax[1], input$sens_nmax[2])
        couts <- 100 / nmax_vals + rnorm(length(nmax_vals), 0, 0.5)
        ggplot(data.frame(nmax = nmax_vals, cout = couts), aes(nmax, cout)) +
            geom_line(color = "#1D9E75", linewidth = 1) +
            geom_point(color = "#D85A30", size = 3) +
            labs(title = "Sensibilité au budget Nmax",
                 x = "Nmax (nb points autorisés)", y = "Distance totale (km)") +
            theme_minimal(base_size = 13)
    })
    
    # ── Distance par véhicule (barplot) ────────────────────────
    output$graphe_vehicules <- renderPlot({
        req(rv$solution)
        df_v <- data.frame(
            vehicule = names(rv$solution$distances_par_vehicule),
            distance = as.numeric(rv$solution$distances_par_vehicule)
        )
        df_v$vehicule <- paste0("V", seq_len(nrow(df_v)))
        COLORS <- c("#378ADD","#D85A30","#1D9E75","#7F77DD","#BA7517",
                    "#E41A1C","#FF7F00","#4DAF4A","#984EA3","#A65628")
        ggplot(df_v, aes(x = vehicule, y = distance, fill = vehicule)) +
            geom_col(width = 0.6, show.legend = FALSE) +
            geom_text(aes(label = paste0(round(distance, 1), " km")),
                      vjust = -0.4, size = 4, fontface = "bold") +
            scale_fill_manual(values = COLORS[seq_len(nrow(df_v))]) +
            labs(title = "Distance parcourue par véhicule", x = NULL, y = "Distance (km)") +
            theme_minimal(base_size = 13) +
            theme(panel.grid.major.x = element_blank())
    })
    
    # ── Badge de progression (sidebar) ─────────────────────────
    output$badge_etapes <- renderUI({
        etape1 <- !is.null(rv$menages)
        etape2 <- !is.null(rv$candidats)
        etape3 <- !is.null(rv$solution)
        span_ok   <- function(txt) tags$span(class = "badge-ok",   txt)
        span_wait <- function(txt) tags$span(class = "badge-wait", txt)
        tagList(
            tags$p(style = "margin:3px 0;",
                   if (etape1) span_ok("\u2713 Données")      else span_wait("\u25cb Données")),
            tags$p(style = "margin:3px 0;",
                   if (etape2) span_ok("\u2713 Clustering")   else span_wait("\u25cb Clustering")),
            tags$p(style = "margin:3px 0;",
                   if (etape3) span_ok("\u2713 Optimisation") else span_wait("\u25cb Optimisation"))
        )
    })
}