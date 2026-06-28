# ============================================================
#  modules/resultats.R
#  Passerelle de synchronisation des indicateurs pour server.R
# ============================================================

library(dplyr)

# ------------------------------------------------------------
#  calculer_couverture()
#  Calcule le taux de couverture réel des ménages par rapport
#  aux points de collecte effectivement ouverts par le MIP.
# ------------------------------------------------------------
calculer_couverture <- function(candidats, menages) {
    if (is.null(candidats) || is.null(menages) || nrow(candidats) == 0 || nrow(menages) == 0) {
        return(0)
    }
    
    # Un point est actif s'il a reçu de la charge
    points_actifs <- candidats %>% 
        dplyr::filter(charge > 0) %>% 
        dplyr::pull(id)
    
    menages_couverts <- menages %>% 
        dplyr::filter(cluster %in% points_actifs)
    
    total_menages <- nrow(menages)
    if (total_menages == 0) return(0)
    
    round(100 * nrow(menages_couverts) / total_menages, 1)
}


# ------------------------------------------------------------
#  extraire_metriques_vehicules()
#  Calcule la distance parcourue par chaque véhicule (pour le barplot)
# ------------------------------------------------------------
extraire_metriques_vehicules <- function(sol, res_dist) {
    if (is.null(sol) || is.null(sol$sol_z) || nrow(sol$sol_z) == 0) {
        return(numeric())
    }
    
    d_matrix <- res_dist$d_matrix
    
    # Calcul rigoureux des distances d'arcs via mapply
    distances <- sol$sol_z %>%
        dplyr::group_by(v) %>%
        dplyr::summarise(
            dist = sum(mapply(function(ii, jj) d_matrix[ii, jj], i, j)),
            .groups = "drop"
        )
    
    # Formatage en vecteur nommé attendu par server.R
    vec_distances <- distances$dist
    names(vec_distances) <- paste0("Vehicule_", distances$v)
    
    vec_distances
}


# ------------------------------------------------------------
#  enrichir_solution_lrp()
#  Fonction adaptatrice : prend la sortie brute de resoudre_lrp()
#  et crée tous les champs dont server.R a impérativement besoin.
# ------------------------------------------------------------
enrichir_solution_lrp <- function(sol, res_dist, candidats, menages) {
    if (is.null(sol)) return(NULL)
    
    # 1. Gestion des alias de noms (cout_total vs distance_totale)
    if (is.null(sol$cout_total) && !is.null(sol$distance_totale)) {
        sol$cout_total <- round(sol$distance_totale, 2)
    } else if (is.null(sol$cout_total)) {
        sol$cout_total <- 0
    }
    
    # 2. Identification des points ouverts depuis les arcs actifs (z)
    # Les candidats sont les indices intermédiaires (excluant 1=Dépôt et n+2=Décharge)
    if (!is.null(sol$sol_z) && nrow(sol$sol_z) > 0) {
        noeuds_visites <- unique(c(sol$sol_z$i, sol$sol_z$j))
        idx_decharge <- max(res_dist$coords_etendues$index)
        
        # Les points candidats réels correspondent à : index - 1
        sol$points_ouverts <- unique(noeuds_visites[noeuds_visites > 1 & noeuds_visites < idx_decharge]) - 1
        sol$nb_points <- length(sol$points_ouverts)
    } else {
        sol$points_ouverts <- integer()
        sol$nb_points <- 0
    }
    
    # 3. Calcul des distances spécifiques par véhicule
    sol$distances_par_vehicule <- extraire_metriques_vehicules(sol, res_dist)
    sol$nb_vehicules <- length(sol$distances_par_vehicule)
    
    # 4. Calcul de la couverture
    sol$couverture <- calculer_couverture(candidats, menages)
    
    # 5. Construction de la table détaillée pour l'UI (DT)
    if (!is.null(sol$sol_z) && nrow(sol$sol_z) > 0) {
        coords <- res_dist$coords_etendues
        sol$tournees_detail <- sol$sol_z %>%
            dplyr::mutate(
                De = coords$label[match(i, coords$index)],
                Vers = coords$label[match(j, coords$index)],
                Distance_km = round(mapply(function(ii, jj) res_dist$d_matrix[ii, jj], i, j), 2)
            ) %>%
            dplyr::select(Vehicule = v, De, Vers, Distance_km) %>%
            dplyr::arrange(Vehicule)
    } else {
        sol$tournees_detail <- data.frame(Vehicule = integer(), De = character(), Vers = character(), Distance_km = numeric())
    }
    
    return(sol)
}