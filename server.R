# ============================================================
#  server.R — LRP Collecte de Déchets à Dakar
#  Auteur : Ousmane LO
#  Version FINALE : UNIQUEMENT OSM, pas de Haversine
#  + corrections cache et diagnostic
# ============================================================

library(shiny)
library(shinydashboard)
library(leaflet)
library(DT)
library(ggplot2)
library(dplyr)
library(geosphere)
library(sf)

# ── Chargement des modules ───────────────────────────────────
source("modules/clustering.R")
source("modules/distances.R")
source("utils/osm_routing.R")
source("modules/cache_graphe.R")
source("modules/modele_mip.R")
source("modules/resultats.R")

server <- function(input, output, session) {
    
    # ── Données réactives ──────────────────────────────────────
    rv <- reactiveValues(
        menages = NULL,
        candidats = NULL,
        inerties = NULL,
        res_dist = NULL,
        solution = NULL,
        log_optim = "",
        matrice_chargee = FALSE,
        graphe = NULL
    )
    
    # ── Chargement du graphe pré-téléchargé ────────────────────
    observe({
        if (graphe_en_cache()) {
            rv$graphe <- charger_graphe_cache()
            message("[Serveur] ✅ Graphe OSM chargé")
            message("[Serveur] Arêtes : ", nrow(rv$graphe))
            message("[Serveur] Nœuds : ", nrow(dodgr_vertices(rv$graphe)))
        } else {
            showNotification(
                "⚠️ Graphe OSM non trouvé. Assurez-vous que data/graphe_dakar.Rds existe.",
                type = "error",
                duration = 10
            )
        }
    })
    
    # ══════════════════════════════════════════════════════════
    #  ONGLET 1 — DONNÉES
    # ══════════════════════════════════════════════════════════
    
    observeEvent(input$charger, {
        req(input$fichier_menages)
        df <- read.csv(input$fichier_menages$datapath)
        names(df) <- tolower(names(df))
        
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
                             popup = ~paste("Ménage", id, "—", poids_dechets, "kg"))
    })
    
    output$resume_donnees <- renderPrint({
        req(rv$menages)
        cat("Nb ménages :", nrow(rv$menages), "\n")
        cat("Poids total :", sum(rv$menages$poids_dechets), "kg\n")
        cat("Emprise lon :", round(range(rv$menages$longitude), 4), "\n")
        cat("Emprise lat :", round(range(rv$menages$latitude), 4), "\n")
    })
    
    # ══════════════════════════════════════════════════════════
    #  ONGLET 2 — CLUSTERING
    # ══════════════════════════════════════════════════════════
    
    # ── Suggestion automatique de k_min / k_max ─────────────────
    # k_min : plancher physique — nb minimal de clusters pour respecter
    #         la capacité par point (cap_point), sinon la contrainte de
    #         charge est mathématiquement infaisable quel que soit k.
    # k_max : le plus petit des trois plafonds suivants :
    #   (a) 3x le budget max envisagé (input$n_max), pour garder de la
    #       marge de choix au MILP (cf. discussion Nmax/couverture) ;
    #   (b) taille d'échantillon / taille minimale raisonnable d'un
    #       cluster (20 ménages, ajustable) — évite des clusters
    #       microscopiques sans réalité opérationnelle ;
    #   (c) plafond de tractabilité du solveur exact (par défaut 50).
    K_MAX_SOLVEUR       <- 50   # plafond de tractabilité MILP
    LAMBDA_MARGE_BUDGET <- 3    # k_max >= lambda * n_max souhaité
    
    observe({
        req(rv$menages, input$cap_point)
        
        poids_total <- sum(rv$menages$poids_dechets)
        m_menages   <- nrow(rv$menages)
        
        # Taille minimale de cluster : pilotable via input$taille_min_cluster
        # si le contrôle existe dans l'UI, sinon repli sur une valeur adaptée
        # à la taille de l'échantillon (evite qu'un seuil fixe de 20 écrase
        # k_max sur de petits jeux de données, ex. 135 menages -> k_max=6).
        taille_min_cluster <- if (!is.null(input$taille_min_cluster)) {
            input$taille_min_cluster
        } else {
            max(5, round(m_menages / 15))  # ~15 clusters de taille egale par defaut
        }
        
        # (k_min) plancher physique lié à la capacité
        k_min_sugg <- ceiling(poids_total / input$cap_point)
        k_min_sugg <- max(k_min_sugg, 2)  # au moins 2 clusters
        
        # (k_max) trois plafonds combinés
        plafond_taille <- floor(m_menages / taille_min_cluster)
        plafond_budget <- if (!is.null(input$n_max) && input$n_max > 0) {
            LAMBDA_MARGE_BUDGET * input$n_max
        } else {
            NA_integer_
        }
        plafond_solveur <- K_MAX_SOLVEUR
        
        k_max_sugg <- min(c(plafond_taille, plafond_budget, plafond_solveur),
                          na.rm = TRUE)
        
        # Garde-fou : si les contraintes se contredisent (k_max < k_min),
        # on élargit k_max plutôt que de bloquer l'utilisateur, avec alerte.
        if (k_max_sugg < k_min_sugg) {
            showNotification(
                paste0("⚠️ cap_point trop faible pour la taille min. de cluster (",
                       taille_min_cluster, " ménages) : k_max ajusté à k_min + 5."),
                type = "warning"
            )
            k_max_sugg <- k_min_sugg + 5
        }
        
        message("[Clustering] taille_min_cluster utilisee = ", taille_min_cluster)
        message("[Clustering] k_min suggéré = ", k_min_sugg,
                " | k_max suggéré = ", k_max_sugg,
                " (plafonds : taille=", plafond_taille,
                ", budget=", plafond_budget,
                ", solveur=", plafond_solveur, ")")
        
        updateSliderInput(session, "k_min",
                          min = 2, max = k_max_sugg, value = k_min_sugg)
        updateSliderInput(session, "k_max",
                          min = k_min_sugg, max = max(k_max_sugg, K_MAX_SOLVEUR),
                          value = k_max_sugg)
    })
    
    observeEvent(input$lancer_cluster, {
        req(rv$menages)
        
        withProgress(message = "Clustering en cours...", value = 0, {
            
            incProgress(0.2, message = "Calcul des inerties K-means++...")
            
            res_cl <- clustering_contraint(
                menages = rv$menages,
                k_min = input$k_min,
                k_max = input$k_max,
                n_init = input$n_init,
                cap_point = input$cap_point,
                appliquer_contrainte = input$contrainte_charge
            )
            
            rv$candidats <- res_cl$candidats
            rv$menages <- res_cl$menages
            rv$inerties <- res_cl$inerties
            
            incProgress(1, message = "Terminé !")
            
            if (length(res_cl$clusters_surcharges) > 0) {
                showNotification(
                    paste0("\u26a0 ", length(res_cl$clusters_surcharges),
                           " point(s) surchargé(s)"),
                    type = "warning"
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
                  options = list(pageLength = 10),
                  colnames = c("ID", "Longitude", "Latitude", "Cluster", "Charge (kg)"))
    })
    
    # ══════════════════════════════════════════════════════════
    #  ONGLET 3 — OPTIMISATION (OSM UNIQUEMENT)
    # ══════════════════════════════════════════════════════════
    
    # ── Calcul de la matrice OSM ──────────────────────────────
    observeEvent(input$calculer_matrice, {
        req(rv$candidats)
        
        if (is.null(rv$graphe)) {
            showNotification(
                "❌ Graphe OSM non disponible. Vérifiez que data/graphe_dakar.Rds existe.",
                type = "error"
            )
            return()
        }
        
        withProgress(message = "Calcul de la matrice OSM...", value = 0, {
            
            incProgress(0.1, message = "Vérification du cache...")
            
            # Essayer de charger depuis le cache utilisateur
            res_dist <- charger_matrice_utilisateur(
                candidats = rv$candidats,
                depot_lon = input$depot_lon_opt,
                depot_lat = input$depot_lat_opt,
                decharge_lon = input$decharge_lon_opt,
                decharge_lat = input$decharge_lat_opt
            )
            
            # 🔧 CORRECTION : Si chargé depuis cache, vérifier/réparer coords_snap
            if (!is.null(res_dist)) {
                # Vérifier que coords_snap existe
                if (is.null(res_dist$coords_snap)) {
                    message("[Serveur] ⚠️ coords_snap manquant dans le cache, tentative de réparation...")
                    if (!is.null(res_dist$coords_etendues)) {
                        res_dist$coords_snap <- res_dist$coords_etendues
                        res_dist$coords_snap$snap_lon <- res_dist$coords_snap$longitude
                        res_dist$coords_snap$snap_lat <- res_dist$coords_snap$latitude
                        res_dist$coords_snap$node_id <- paste0("node_", res_dist$coords_snap$index)
                        message("[Serveur] ✅ coords_snap réparé à partir de coords_etendues")
                    } else {
                        message("[Serveur] ❌ Impossible de réparer coords_snap")
                    }
                }
                # Vérifier que le graphe est présent
                if (is.null(res_dist$graphe)) {
                    message("[Serveur] ⚠️ graphe manquant dans le cache, ajout du graphe courant")
                    res_dist$graphe <- rv$graphe
                }
                # S'assurer que la source est osm
                if (is.null(res_dist$source)) {
                    res_dist$source <- "osm"
                }
            }
            
            # Si pas en cache ou réparation échouée, calculer avec OSM
            if (is.null(res_dist) || is.null(res_dist$coords_snap)) {
                incProgress(0.3, message = "Calcul OSM en cours...")
                
                res_dist <- matrice_distances_osm(
                    candidats = rv$candidats,
                    depot_lon = input$depot_lon_opt,
                    depot_lat = input$depot_lat_opt,
                    decharge_lon = input$decharge_lon_opt,
                    decharge_lat = input$decharge_lat_opt,
                    graphe = rv$graphe
                )
                
                # Stocker le graphe dans res_dist pour les géométries
                res_dist$graphe <- rv$graphe
                
                # Sauvegarder dans le cache
                incProgress(0.1, message = "Sauvegarde en cache...")
                sauvegarder_matrice_utilisateur(
                    res_dist = res_dist,
                    candidats = rv$candidats,
                    depot_lon = input$depot_lon_opt,
                    depot_lat = input$depot_lat_opt,
                    decharge_lon = input$decharge_lon_opt,
                    decharge_lat = input$decharge_lat_opt
                )
                
                # Vérifier que c'est bien OSM
                if (is.null(res_dist$source) || res_dist$source != "osm") {
                    showNotification(
                        "❌ ERREUR : La matrice n'est pas OSM !",
                        type = "error",
                        duration = 10
                    )
                    return()
                }
                
                showNotification("✅ Matrice OSM calculée et mise en cache !", type = "message")
                
            } else {
                # Vérifier que le cache contient bien OSM
                if (!is.null(res_dist$source) && res_dist$source == "osm") {
                    showNotification("📦 Matrice OSM chargée depuis le cache !", type = "message")
                } else {
                    showNotification(
                        "❌ ERREUR : Le cache contient une matrice Haversine !",
                        type = "error",
                        duration = 10
                    )
                    return()
                }
            }
            
            rv$res_dist <- res_dist
            rv$matrice_chargee <- TRUE
            
            incProgress(1, message = "Terminé !")
        })
    })
    
    # ── État de la matrice ──────────────────────────────────────
    output$statut_matrice <- renderPrint({
        if (rv$matrice_chargee) {
            if (!is.null(rv$res_dist$source) && rv$res_dist$source == "osm") {
                cat("✅ Matrice OSM chargée\n")
                cat("  - Points candidats :", rv$res_dist$n, "\n")
                cat("  - Taille :", nrow(rv$res_dist$d_matrix), "x", ncol(rv$res_dist$d_matrix), "\n")
                cat("  - Source :", rv$res_dist$source, " (DISTANCES RÉELLES)\n")
                cat("  - coords_snap présent : ", ifelse("coords_snap" %in% names(rv$res_dist), "✅", "❌"), "\n")
                if ("coords_snap" %in% names(rv$res_dist)) {
                    cat("  - coords_snap : ", nrow(rv$res_dist$coords_snap), " lignes\n")
                }
                cat("  - Dépôt :", input$depot_lon_opt, ",", input$depot_lat_opt, "\n")
                cat("  - Décharge :", input$decharge_lon_opt, ",", input$decharge_lat_opt, "\n")
            } else {
                cat("❌ ERREUR : Matrice non OSM !\n")
                cat("  - Source :", rv$res_dist$source %||% "inconnue", "\n")
            }
        } else {
            cat("⏳ En attente du calcul de la matrice OSM...\n")
            cat("📌 Entrez les coordonnées et cliquez sur 'Calculer la matrice OSM'")
        }
    })
    
    # ── Aperçu de la matrice ────────────────────────────────────
    output$apercu_matrice <- renderDT({
        req(rv$matrice_chargee)
        
        df <- as.data.frame(rv$res_dist$d_matrix)
        rownames(df) <- rv$res_dist$coords_etendues$label
        colnames(df) <- rv$res_dist$coords_etendues$label
        
        datatable(round(df, 2),
                  options = list(pageLength = 10, scrollX = TRUE),
                  caption = "Matrice de distances OSM (km)")
    })
    
    # ── Lancer l'optimisation ──────────────────────────────────
    observeEvent(input$lancer_optim, {
        req(rv$menages, rv$candidats)
        
        if (!rv$matrice_chargee) {
            showNotification(
                "⚠️ Veuillez d'abord calculer la matrice OSM !",
                type = "warning"
            )
            return()
        }
        
        # Vérifier que c'est bien OSM
        if (is.null(rv$res_dist$source) || rv$res_dist$source != "osm") {
            showNotification(
                "❌ ERREUR : La matrice n'est pas OSM ! Recalculez la matrice.",
                type = "error"
            )
            return()
        }
        
        withProgress(message = "Optimisation en cours...", value = 0, {
            
            incProgress(0.3, message = "Lancement du solveur MIP...")
            
            sol <- tryCatch(
                resoudre_lrp(
                    res_dist = rv$res_dist,
                    menages = rv$menages,
                    nb_vehicules = input$nb_vehicules,
                    cap_camion = input$cap_camion,
                    cap_point = input$cap_point,
                    d_max = input$d_max,
                    n_max = input$n_max,
                    solveur = input$solveur,
                    time_limit = input$time_limit,
                    tol_couverture = if (!is.null(input$tol_couverture)) {
                        input$tol_couverture
                    } else {
                        0  # valeur par défaut si le slider n'est pas encore dans l'UI
                    }
                ),
                error = function(e) {
                    showNotification(paste("Erreur solveur :", e$message), type = "error")
                    NULL
                }
            )
            
            if (!is.null(sol)) {
                rv$solution <- sol
                rv$log_optim <- paste0(
                    "Solveur : ", input$solveur, "\n",
                    "Statut : ", sol$statut, "\n",
                    "Distance des tournées (étape B) : ", sol$cout_total, " km\n",
                    "Couverture max. atteignable (étape A, plafond) : ",
                    sol$couverture_max_atteignable, " %\n",
                    "Tolérance couverture utilisée : ", sol$tol_couverture_utilisee, " ménage(s)\n",
                    "Points ouverts : ", sol$nb_points, "\n",
                    "Véhicules utilisés : ", sol$nb_vehicules_utilises, " / ", sol$nb_vehicules, "\n",
                    "Couverture réelle obtenue : ", sol$couverture, " % (",
                    sol$nb_non_couverts, " ménage(s) au-delà de D_max)\n\n",
                    "Distance par véhicule :\n",
                    paste(names(sol$distances_par_vehicule),
                          round(sol$distances_par_vehicule, 2), "km",
                          collapse = "\n")
                )
                showNotification("✅ Optimisation terminée avec succès !", type = "message")
            }
            
            incProgress(1, message = "Terminé !")
        })
    })
    
    output$statut_optim <- renderPrint({
        if (is.null(rv$solution)) {
            if (!rv$matrice_chargee) {
                cat("⏳ Calculez d'abord la matrice OSM")
            } else {
                cat("En attente du lancement...")
            }
        } else {
            cat("✓", rv$solution$statut, "\nCoût total :", rv$solution$cout_total, "km")
        }
    })
    
    output$log_solveur <- renderPrint({ cat(rv$log_optim) })
    
    output$contraintes_actives <- renderPrint({
        cat("── ÉTAPE A (couverture maximale, modèle léger) ──────────\n")
        cat("CA1 : Σⱼ xᵢⱼ = 1                    ∀i  (affectation unique)\n")
        cat("CA2 : xᵢⱼ ≤ yⱼ                      ∀i,j (seulement si ouvert)\n")
        cat("CA3 : Σⱼ yⱼ ≤ Nmax                       (budget)\n")
        cat("CA4 : Σᵢ wᵢ·xᵢⱼ ≤ Cⱼ·yⱼ              ∀j  (capacité point)\n")
        cat("CA5 : cᵢ ≤ Σ_{j: dᵢⱼ≤Dmax} xᵢⱼ        ∀i  (couverture)\n")
        cat("      objectif : max Σᵢ cᵢ  →  couverture_max\n\n")
        cat("── ÉTAPE B (routage optimal, modèle complet) ────────────\n")
        cat("C1  : Σⱼ z^v_Oj ≤ 1              ∀v  (départ optionnel)\n")
        cat("C2  : Σᵢ z^v_iS ≤ 1              ∀v  (arrivée optionnelle)\n")
        cat("C2b : Σⱼ z^v_Oj = Σᵢ z^v_iS      ∀v  (cohérence départ/arrivée)\n")
        cat("C3  : Σᵢ z^v_ij = Σₖ z^v_jk    ∀v,j (conservation flux)\n")
        cat("C4  : Σᵥ Σᵢ z^v_ij = yⱼ        ∀j  (visite si ouvert)\n")
        cat("C5  : z^v_ii = 0                ∀v,i (pas d'auto-boucle)\n")
        cat("C6  : z^v_Sj = 0                ∀v,j (pas d'arc depuis décharge)\n")
        cat("C7  : z^v_iO = 0                ∀v,i (pas d'arc vers dépôt)\n")
        cat("C8  : Σⱼ xᵢⱼ = 1               ∀i  (affectation unique — toujours obligatoire,\n")
        cat("      un ménage est toujours collecté, même au-delà de Dmax)\n")
        cat("C9  : xᵢⱼ ≤ yⱼ                 ∀i,j (seulement si ouvert)\n")
        cat("C10 : cᵢ ≤ Σ_{j: dᵢⱼ≤Dmax} xᵢⱼ   ∀i  (indicateur couverture, idem étape A)\n")
        cat("C10b: Σᵢ cᵢ ≥ couverture_max - tol_couverture  (couverture garantie,\n")
        cat("      ancrée sur l'optimum de l'étape A — remplace l'ancienne pénalité\n")
        cat("      pondérée poids_desserte, abandonnée : cf. document du modèle)\n")
        cat("C11 : Σᵢ wᵢ·xᵢⱼ ≤ Cⱼ·yⱼ       ∀j  (capacité point)\n")
        cat("C12 : Σⱼ yⱼ ≤ Nmax                  (budget)\n")
        cat("C13 : MTZ  u^v_j ≥ u^v_i + Cⱼ·z^v_ij - Qᵥ(1-z^v_ij)  (sous-tours)\n")
        cat("C14 : u^v_j ≤ Qᵥ               ∀v,j (borne capacité)\n")
        cat("C15 : u^v_O = 0                ∀v  (dépôt vide au départ)\n\n")
        cat("Objectif étape B : min Σᵥ Σᵢ Σⱼ dᵢⱼ·z^v_ij  (distance des tournées,\n")
        cat("uniquement — la couverture est déjà garantie par C10b, pas besoin\n")
        cat("de la mélanger dans l'objectif)\n")
    })
    
    # ══════════════════════════════════════════════════════════
    #  ONGLET 4 — RÉSULTATS (OSM UNIQUEMENT)
    # ══════════════════════════════════════════════════════════
    
    output$kpi_distance <- renderValueBox({
        valueBox(
            value = if (!is.null(rv$solution)) paste(rv$solution$cout_total, "km") else "—",
            subtitle = "Distance totale (OSM)",
            icon = icon("road"), color = "blue"
        )
    })
    
    output$kpi_points <- renderValueBox({
        valueBox(
            value = if (!is.null(rv$solution)) rv$solution$nb_points else "—",
            subtitle = "Points ouverts",
            icon = icon("map-marker"), color = "yellow"
        )
    })
    
    output$kpi_vehicules <- renderValueBox({
        valueBox(
            value = if (!is.null(rv$solution)) {
                paste0(rv$solution$nb_vehicules_utilises, " / ", rv$solution$nb_vehicules)
            } else "—",
            subtitle = "Véhicules utilisés",
            icon = icon("truck"), color = "green"
        )
    })
    
    output$kpi_couverture <- renderValueBox({
        valueBox(
            value = if (!is.null(rv$solution)) paste(rv$solution$couverture, "%") else "—",
            subtitle = if (!is.null(rv$solution)) {
                paste0("Ménages couverts (max atteignable sous ce budget : ",
                       rv$solution$couverture_max_atteignable, "%)")
            } else {
                "Ménages couverts"
            },
            icon = icon("home"), color = "purple"
        )
    })
    
    # ── Carte des tournées OSM ─────────────────────────────────
    output$carte_tournees <- renderLeaflet({
        req(rv$solution, rv$candidats, rv$res_dist)
        
        # 🔍 DIAGNOSTIC
        message("\n🔍 DIAGNOSTIC CARTE TOURNEES")
        message("─────────────────────────────")
        message("source : ", rv$res_dist$source %||% "inconnue")
        message("coords_snap présent : ", "coords_snap" %in% names(rv$res_dist))
        message("graphe présent : ", "graphe" %in% names(rv$res_dist))
        
        if ("coords_snap" %in% names(rv$res_dist)) {
            message("coords_snap : ", nrow(rv$res_dist$coords_snap), " lignes")
            message("Premiers points snappés :")
            print(head(rv$res_dist$coords_snap[, c("index", "label", "snap_lon", "snap_lat")]))
        }
        
        if (!is.null(rv$solution$sol_z)) {
            message("sol_z : ", nrow(rv$solution$sol_z), " arcs")
            message("Premiers arcs :")
            print(head(rv$solution$sol_z))
        }
        message("─────────────────────────────\n")
        
        # Vérifier que c'est OSM
        if (is.null(rv$res_dist$source) || rv$res_dist$source != "osm") {
            showNotification(
                "❌ ERREUR : Les distances ne sont pas OSM !",
                type = "error"
            )
            return(leaflet() %>% addTiles() %>% 
                       addControl("❌ Erreur : Distances non OSM", position = "topcenter"))
        }
        
        COLORS <- c("#378ADD", "#D85A30", "#1D9E75", "#7F77DD", "#BA7517",
                    "#E41A1C", "#FF7F00", "#4DAF4A", "#984EA3", "#A65628")
        
        m <- leaflet() %>% addTiles()
        
        # Récupérer les arcs actifs
        if (!is.null(rv$solution$sol_z) && nrow(rv$solution$sol_z) > 0) {
            
            arcs <- rv$solution$sol_z
            
            message("[Carte] 🗺️ Génération des itinéraires OSM pour ", nrow(arcs), " arcs")
            
            # Générer les géométries OSM
            geos <- tryCatch(
                geometries_tournees_osm(arcs, rv$res_dist),
                error = function(e) {
                    message("[Carte] ❌ Erreur geometries_tournees_osm : ", e$message)
                    NULL
                }
            )
            
            if (!is.null(geos) && length(geos) > 0) {
                message("[Carte] ✅ ", length(geos), " géométries OSM générées")
                
                # Tracer les itinéraires OSM
                for (geo in geos) {
                    if (is.null(geo)) next
                    
                    v_id <- geo$v
                    couleur <- COLORS[(v_id - 1) %% length(COLORS) + 1]
                    
                    tryCatch({
                        coords_geo <- st_coordinates(geo$geometry)
                        
                        if (nrow(coords_geo) > 1) {
                            if (!any(is.na(coords_geo))) {
                                m <- m %>% addPolylines(
                                    lng = coords_geo[, 1],
                                    lat = coords_geo[, 2],
                                    color = couleur,
                                    weight = 4,
                                    opacity = 0.9,
                                    label = paste0("Véhicule ", v_id, " : ",
                                                   geo$label_i, " → ", geo$label_j,
                                                   " (", geo$distance_km, " km)")
                                )
                            }
                        }
                    }, error = function(e) {
                        message("[Carte] ⚠️ Erreur tracé OSM : ", e$message)
                    })
                }
            } else {
                message("[Carte] ❌ Aucune géométrie OSM générée")
                showNotification(
                    "⚠️ Impossible de générer les itinéraires OSM. Vérifiez le graphe.",
                    type = "error",
                    duration = 10
                )
            }
        }
        
        # Points candidats ouverts
        if (!is.null(rv$solution$points_ouverts) && length(rv$solution$points_ouverts) > 0) {
            points_ouverts_df <- rv$candidats[rv$candidats$id %in% rv$solution$points_ouverts, ]
            if (nrow(points_ouverts_df) > 0) {
                m <- m %>% addCircleMarkers(
                    data = points_ouverts_df,
                    ~longitude, ~latitude,
                    radius = 8,
                    color = "#333333",
                    fillColor = "#FFC107",
                    fillOpacity = 0.9,
                    weight = 2,
                    popup = ~paste0("Point ", id, "<br>Charge : ", round(charge), " kg")
                )
            }
        }
        
        # Dépôt et décharge
        m <- m %>%
            addAwesomeMarkers(
                lng = input$depot_lon_opt, lat = input$depot_lat_opt,
                icon = awesomeIcons(icon = "home", markerColor = "green", library = "fa"),
                popup = "Dépôt (O)"
            ) %>%
            addAwesomeMarkers(
                lng = input$decharge_lon_opt, lat = input$decharge_lat_opt,
                icon = awesomeIcons(icon = "trash", markerColor = "red", library = "fa"),
                popup = "Décharge (S)"
            )
        
        # Légende
        if (!is.null(rv$solution$nb_vehicules) && rv$solution$nb_vehicules > 0) {
            nb_v <- min(rv$solution$nb_vehicules, length(COLORS))
            m <- m %>% addLegend(
                position = "bottomright",
                colors = COLORS[1:nb_v],
                labels = paste("Véhicule", 1:nb_v),
                title = "Tournées OSM"
            )
        }
        
        m
    })
    
    # ── Table des tournées ──────────────────────────────────────
    output$table_tournees <- renderDT({
        req(rv$solution)
        datatable(
            rv$solution$tournees_detail,
            options = list(pageLength = 15, order = list(list(0, "asc"))),
            colnames = c("Véhicule", "De", "Vers", "Distance (km)")
        )
    })
    
    # ── Export CSV ──────────────────────────────────────────────
    output$export_resultats <- downloadHandler(
        filename = function() paste0("resultats_LRP_Dakar_", Sys.Date(), ".csv"),
        content = function(file) {
            req(rv$solution)
            write.csv(rv$solution$tournees_detail, file, row.names = FALSE)
        }
    )
    
    # ── Analyse de sensibilité ──────────────────────────────────
    output$graphe_sensibilite <- renderPlot({
        req(rv$solution, rv$res_dist, rv$menages, rv$candidats)
        
        nb_candidats <- nrow(rv$candidats)
        nmax_demande <- seq(input$sens_nmax[1], input$sens_nmax[2])
        
        # Au-delà de nb_candidats, Σy_j <= n_max cesse d'être contraignante
        # (il n'y a tout simplement plus de points à ouvrir) -- la courbe
        # plafonne nécessairement à partir de là. Avertir plutôt que de
        # laisser croire à un bug si la plage demandée dépasse ce seuil.
        if (max(nmax_demande) > nb_candidats) {
            showNotification(
                paste0("ℹ️ Seulement ", nb_candidats, " points candidats disponibles : ",
                       "la courbe sera plate au-delà de Nmax=", nb_candidats,
                       " (relancer le clustering avec un k_max plus élevé pour tester plus loin)."),
                type = "message", duration = 8
            )
        }
        nmax_vals <- nmax_demande[nmax_demande <= max(nb_candidats, nmax_demande[1])]
        if (length(nmax_vals) == 0) nmax_vals <- nmax_demande
        
        n_total <- length(nmax_vals)
        withProgress(message = "Analyse de sensibilité (Nmax)...", value = 0, {
            df_sens <- sensibilite_nmax(
                sol = rv$solution,
                nmax_vals = nmax_vals,
                res_dist = rv$res_dist,
                menages = rv$menages,
                resoudre_fn = resoudre_lrp,
                progress_cb = function(k, n, nm) {
                    incProgress(1 / n_total,
                                detail = paste0(k, "/", n_total, " (Nmax=", nm, ")"))
                },
                nb_vehicules = input$nb_vehicules,
                cap_camion = input$cap_camion,
                cap_point = input$cap_point,
                d_max = input$d_max,
                solveur = input$solveur,
                time_limit = 30,             # Temps réduit pour la sensibilité (étape B)
                time_limit_couverture = 15,  # Temps réduit pour l'étape A
                tol_couverture = if (!is.null(input$tol_couverture)) input$tol_couverture else 0
            )
        })
        
        graphe_sensibilite_nmax(df_sens, nmax_ref = rv$solution$nb_points)
    })
    
    # ── Distance par véhicule ──────────────────────────────────
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
            labs(title = "Distance parcourue par véhicule (OSM)", x = NULL, y = "Distance (km)") +
            theme_minimal(base_size = 13) +
            theme(panel.grid.major.x = element_blank())
    })
    
    # ── Badge de progression ────────────────────────────────────
    output$badge_etapes <- renderUI({
        etape1 <- !is.null(rv$menages)
        etape2 <- !is.null(rv$candidats)
        etape3 <- !is.null(rv$solution)
        etape4 <- rv$matrice_chargee && !is.null(rv$res_dist$source) && rv$res_dist$source == "osm"
        
        span_ok <- function(txt) tags$span(class = "badge-ok", txt)
        span_wait <- function(txt) tags$span(class = "badge-wait", txt)
        span_error <- function(txt) tags$span(class = "badge-error", txt)
        
        tagList(
            tags$p(style = "margin:3px 0;",
                   if (etape1) span_ok("\u2713 Données") else span_wait("\u25cb Données")),
            tags$p(style = "margin:3px 0;",
                   if (etape2) span_ok("\u2713 Clustering") else span_wait("\u25cb Clustering")),
            tags$p(style = "margin:3px 0;",
                   if (etape4) span_ok("\u2713 Matrice OSM") else if (rv$matrice_chargee) span_error("\u2717 Matrice non OSM") else span_wait("\u25cb Matrice OSM")),
            tags$p(style = "margin:3px 0;",
                   if (etape3) span_ok("\u2713 Optimisation") else span_wait("\u25cb Optimisation"))
        )
    })
}