# ============================================================
#  modules/distances.R
#  Calcul de la matrice de distances étendue (Haversine)
#  
#  ⚠️ Ce n'est qu'un FALLBACK si OSM n'est pas disponible
#  La méthode PRINCIPALE est OSM dans osm_routing.R
# ============================================================

library(geosphere)

construire_matrice_distances <- function(candidats,
                                         depot_lon, depot_lat,
                                         decharge_lon, decharge_lat) {
    
    n <- nrow(candidats)
    
    coords_etendues <- data.frame(
        index     = c(1, seq(2, n + 1), n + 2),
        role      = c("depot", rep("candidat", n), "decharge"),
        id_local  = c(0, seq_len(n), n + 1),
        longitude = c(depot_lon,    candidats$longitude, decharge_lon),
        latitude  = c(depot_lat,    candidats$latitude,  decharge_lat),
        label     = c("Dépôt (O)", paste0("Point_", seq_len(n)), "Décharge (S)")
    )
    
    coords_mat <- as.matrix(coords_etendues[, c("longitude", "latitude")])
    d_metres <- distm(coords_mat, coords_mat, fun = distHaversine)
    d_matrix <- d_metres / 1000
    diag(d_matrix) <- 0
    
    list(
        d_matrix        = d_matrix,
        coords_etendues = coords_etendues,
        n               = n,
        idx_depot       = 1L,
        idx_decharge    = as.integer(n + 2)
    )
}