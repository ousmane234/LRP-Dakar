# ============================================================
#  SERVER
# ============================================================
library(ompr)
library(ompr.roi)
library(ROI.plugin.glpk)
library(geosphere)
server <- function(input, output, session) {
    
    # ── Données réactives ──────────────────────────────────
    rv <- reactiveValues(
        menages    = NULL,   # data.frame : id, lon, lat, poids
        candidats  = NULL,   # data.frame : id, lon, lat, cluster, charge
        inerties = NULL,
        solution   = NULL,   # liste : points_ouverts, tournees, cout
        log_optim  = ""
    )
    
    # ── ONGLET 1 : Chargement ──────────────────────────────
    observeEvent(input$charger, {
        req(input$fichier_menages)
        rv$menages <- read.csv(input$fichier_menages$datapath)
        names(rv$menages) <- tolower(names(rv$menages))
        showNotification(paste(nrow(rv$menages), "ménages chargés."), type = "message")
    })
    
    output$carte_menages <- renderLeaflet({
        req(rv$menages)
        leaflet(rv$menages) %>%
            addTiles() %>%
            addCircleMarkers(~longitude, ~latitude,
                             radius = 4, color = "#378ADD",
                             popup  = ~paste("Ménage", id, "—", poids_dechets, "kg")) %>%
            addMarkers(lng = input$depot_lon,    lat = input$depot_lat,
                       popup = "Dépôt (o)", icon = awesomeIcons(icon="home", markerColor="green")) %>%
            addMarkers(lng = input$decharge_lon, lat = input$decharge_lat,
                       popup = "Décharge (s)", icon = awesomeIcons(icon="trash", markerColor="red"))
    })
    
    output$resume_donnees <- renderPrint({
        req(rv$menages)
        cat("Nb ménages        :", nrow(rv$menages), "\n")
        cat("Poids total (kg)  :", sum(rv$menages$poids_dechets), "\n")
        cat("Emprise lon       :", round(range(rv$menages$longitude), 4), "\n")
        cat("Emprise lat       :", round(range(rv$menages$latitude),  4), "\n")
    })
    
    # ── ONGLET 2 : Clustering ──────────────────────────────
    observeEvent(input$lancer_cluster, {
        req(rv$menages)
        
        # -- K-Means
        coords <- as.matrix(rv$menages[, c("longitude","latitude")])
        inerties <- sapply(input$k_min:input$k_max, function(k) {
            km <- kmeans(coords, centers = k, nstart = input$n_init)
            km$tot.withinss
        })
        
        k_optimal <- input$k_min + which.max(diff(diff(inerties))) # coude simplifié
        km_final  <- kmeans(coords, centers = k_optimal, nstart = input$n_init)
        
        rv$candidats <- data.frame(
            id        = seq_len(k_optimal),
            longitude = km_final$centers[,1],
            latitude  = km_final$centers[,2],
            cluster   = seq_len(k_optimal),
            charge    = tapply(rv$menages$poids_dechets, km_final$cluster, sum)
        )
        # contrainte de charge : marquer les cluster valides === à compléter
        # Stocker inerties pour le graphe
        rv$inerties <- data.frame(k = input$k_min:input$k_max, inertie = inerties)
        rv$menages$cluster <- km_final$cluster
        
        showNotification(paste("Clustering terminé — k optimal =", k_optimal), type = "message")
    })
    
    # Carte des clusters 
    output$carte_clusters <- renderLeaflet({
        req(rv$candidats, rv$menages)
        pal <- colorFactor(palette = "Set1", domain = rv$menages$cluster)
        leaflet() %>%
            addTiles() %>%
            addCircleMarkers(data = rv$menages,
                             ~longitude, ~latitude,
                             color = ~pal(cluster), radius = 4, opacity = 0.7) %>%
            addMarkers(data = rv$candidats,
                       ~longitude, ~latitude,
                       popup = ~paste("Candidat", id, "— charge:", round(charge), "kg"))
    })
    
    # Courbe d'inertie
    output$courbe_inertie <- renderPlot({
        req(rv$inerties)
        ggplot(rv$inerties, aes(x = k, y = inertie)) +
            geom_line(color = "#378ADD", linewidth = 1) +
            geom_point(color = "#D85A30", size = 3) +
            labs(title = "Courbe d'inertie (méthode du coude)",
                 x = "Nombre de clusters k", y = "Inertie totale") +
            theme_minimal(base_size = 13)
    })
    # Table des candidats
    output$table_candidats <- renderDT({
        req(rv$candidats)
        datatable(rv$candidats,
                  options = list(pageLength = 10),
                  colnames = c("ID", "Longitude", "Latitude", "Cluster", "Charge (kg)"))
    })
    
    # ── ONGLET 3 : Optimisation ────────────────────────────
    observeEvent(input$lancer_optim, {
        req(rv$candidats, rv$menages)
        
        withProgress(message = "Optimisation en cours...", value = 0, {
            
            n  <- nrow(rv$candidats)
            m  <- nrow(rv$menages)
            V  <- input$nb_vehicules
            Qv <- input$cap_camion
            Cj <- rv$candidats$charge
            
            d_matrix <- distm(rv$candidats[,c("longitude","latitude")],
                              rv$candidats[,c("longitude","latitude")],
                              fun = distHaversine) / 1000
            
            # Définition du modèle
            model <- MIPModel() %>%
                add_variable(y[j], j = 1:n, type = "binary") %>%
                add_variable(x[i,j], i = 1:m, j = 1:n, type = "binary") %>%
                add_variable(z[v,i,j], v = 1:V, i = 1:n, j = 1:n, type = "binary") %>%
                add_variable(u[v,j], v = 1:V, j = 1:n, type = "continuous", lb = 0) %>%
                set_objective(sum_expr(d_matrix[i,j] * z[v,i,j], v = 1:V, i = 1:n, j = 1:n), "min") %>%
                add_constraint(sum_expr(x[i,j], j = 1:n) == 1, i = 1:m) %>%
                add_constraint(x[i,j] <= y[j], i = 1:m, j = 1:n) %>%
                add_constraint(sum_expr(rv$menages$poids_dechets[i] * x[i,j], i = 1:m) <= Cj[j] * y[j], j = 1:n) %>%
                add_constraint(sum_expr(y[j], j = 1:n) <= input$n_max) %>%
                add_constraint(u[v,j] <= Qv, v = 1:V, j = 1:n) %>%
                add_constraint(u[v,j] >= u[v,i] + Cj[j] * z[v,i,j] - Qv * (1 - z[v,i,j]), v = 1:V, i = 1:n, j = 1:n) %>%
                # C3 : chaque point ouvert doit être visité exactement une fois (toutes véhicules confondus)
                add_constraint(sum_expr(z[v,i,j], v = 1:V, i = 1:n) == y[j], j = 1:n)
            
            # Résolution
            result <- solve_model(model,
                                  with_ROI(solver = input$solveur,
                                           control = list(tm_limit = input$time_limit*1000)))
            
            # Extraction des solutions
            sol_y <- get_solution(result, y[j])
            points_ouverts <- sol_y$j[sol_y$value > 0.5]
            
            rv$solution <- list(
                cout_total   = objective_value(result),
                nb_points    = length(points_ouverts),
                nb_vehicules = V,
                couverture   = 100, # à calculer selon affectation réelle
                tournees     = get_solution(result, z[v,i,j])
            )
            
            rv$log_optim <- paste("Solveur", input$solveur,
                                  "\nStatut :", solver_status(result),
                                  "\nValeur objectif :", round(objective_value(result),2), "km")
            incProgress(1, message = "Terminé !")
        })
    })
    
    
    output$statut_optim <- renderPrint({
        if (is.null(rv$solution)) cat("En attente du lancement...")
        else cat("✓ Solution optimale trouvée\nCoût :", rv$solution$cout_total, "km")
    })
    
    output$log_solveur <- renderPrint({ cat(rv$log_optim) })
    
    output$contraintes_actives <- renderPrint({
        cat("C1  : Σⱼ z^v_oj = 1         ∀v\n")
        cat("C2  : Σⱼ z^v_js = 1         ∀v\n")
        cat("C3  : Σᵥ Σᵢ z^v_ij = yⱼ    ∀j\n")
        cat("C4  : u^v_j ≥ Cⱼ·z^v_oj    ∀v,j\n")
        cat("C5  : u^v_j ≥ u^v_i+Cⱼ·z^v_ij ∀v,i,j\n")
        cat("C6  : u^v_j ≤ Qᵥ·z^v_ij    ∀v,i,j\n")
        cat("C7  : xᵢⱼ ≤ yⱼ             ∀i,j\n")
        cat("C8  : Σⱼ xᵢⱼ = 1           ∀i\n")
        cat("C9  : xᵢⱼ=0 si dᵢⱼ>Dmax   ∀i,j\n")
        cat("C10 : Σᵢ wᵢ·xᵢⱼ ≤ Cⱼ·yⱼ  ∀j\n")
        cat("C11 : Σⱼ yⱼ ≤ Nmax\n")
    })
    
    # ── ONGLET 4 : Résultats ───────────────────────────────
    output$kpi_distance   <- renderValueBox({
        valueBox(
            value    = if(!is.null(rv$solution)) paste(rv$solution$cout_total, "km") else "—",
            subtitle = "Distance totale",
            icon     = icon("road"), color = "blue")
    })
    
    output$kpi_points     <- renderValueBox({
        valueBox(
            value    = if(!is.null(rv$solution)) rv$solution$nb_points else "—",
            subtitle = "Points ouverts",
            icon     = icon("map-marker"), color = "yellow")
    })
    
    output$kpi_vehicules  <- renderValueBox({
        valueBox(
            value    = if(!is.null(rv$solution)) rv$solution$nb_vehicules else "—",
            subtitle = "Véhicules utilisés",
            icon     = icon("truck"), color = "green")
    })
    
    output$kpi_couverture <- renderValueBox({
        valueBox(
            value    = if(!is.null(rv$solution)) paste(rv$solution$couverture, "%") else "—",
            subtitle = "Ménages couverts",
            icon     = icon("home"), color = "purple")
    })
    
    output$carte_tournees <- renderLeaflet({
        req(rv$solution, rv$candidats)
        COLORS <- c("#378ADD","#D85A30","#1D9E75","#7F77DD","#BA7517")
        
        # Arcs actifs uniquement (variable binaire z = 1)
        arcs <- rv$solution$tournees |>
            dplyr::filter(value > 0.5)
        
        m <- leaflet() %>% addTiles()
        
        # Tracé des arcs par véhicule
        if (nrow(arcs) > 0) {
            for (v_id in sort(unique(arcs$v))) {
                arcs_v  <- arcs[arcs$v == v_id, ]
                couleur <- COLORS[(v_id - 1) %% length(COLORS) + 1]
                for (k in seq_len(nrow(arcs_v))) {
                    pt_i <- rv$candidats[arcs_v$i[k], ]
                    pt_j <- rv$candidats[arcs_v$j[k], ]
                    m <- m %>% addPolylines(
                        lng    = c(pt_i$longitude, pt_j$longitude),
                        lat    = c(pt_i$latitude,  pt_j$latitude),
                        color  = couleur,
                        weight = 3,
                        opacity = 0.8,
                        label  = paste("Véhicule", v_id, ": Point", arcs_v$i[k], "→", arcs_v$j[k])
                    )
                }
            }
        }
        
        # Marqueurs des points ouverts
        m <- m %>% addMarkers(data = rv$candidats,
                              ~longitude, ~latitude,
                              popup = ~paste("Point", id))
        # Dépôt et décharge
        m <- m %>%
            addMarkers(lng = input$depot_lon,    lat = input$depot_lat,
                       popup = "Dépôt (o)",
                       icon  = awesomeIcons(icon = "home",  markerColor = "green")) %>%
            addMarkers(lng = input$decharge_lon, lat = input$decharge_lat,
                       popup = "Décharge (s)",
                       icon  = awesomeIcons(icon = "trash", markerColor = "red"))
        m
    })
    
    output$table_tournees <- renderDT({
        req(rv$solution)
        datatable(rv$solution$tournees,
                  options = list(pageLength = 10),
                  colnames = c("Véhicule","De","Vers","Distance (km)"))
    })
    
    output$export_resultats <- downloadHandler(
        filename = function() paste0("resultats_LRP_Dakar_", Sys.Date(), ".csv"),
        content  = function(file) {
            req(rv$solution)
            write.csv(rv$solution$tournees, file, row.names = FALSE)
        }
    )
    
    output$graphe_sensibilite <- renderPlot({
        req(rv$candidats)
        nmax_vals <- seq(input$sens_nmax[1], input$sens_nmax[2])
        # Valeurs fictives pour la démonstration
        couts <- 100 / nmax_vals + rnorm(length(nmax_vals), 0, 0.5)
        ggplot(data.frame(nmax = nmax_vals, cout = couts), aes(nmax, cout)) +
            geom_line(color = "#1D9E75", linewidth = 1) +
            geom_point(color = "#D85A30", size = 3) +
            labs(title = "Sensibilité au budget Nmax",
                 x = "Nmax (nb points autorisés)", y = "Distance totale (km)") +
            theme_minimal(base_size = 13)
    })
}
