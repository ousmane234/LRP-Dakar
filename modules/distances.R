# ============================================================
#  modules/distances.R
#  Calcul de la matrice de distances étendue
#
#  Nœuds du graphe étendu :
#    - Nœud 0       : Dépôt (O)
#    - Nœuds 1..n   : Points candidats (résultat du clustering)
#    - Nœud n+1     : Décharge (S)
#
#  La matrice retournée est de taille (n+2) x (n+2)
#  avec des distances en kilomètres (Haversine pour l'instant,
#  remplacé par dodgr à l'étape 5).
#
#  Fonction principale : construire_matrice_distances()
#  Retourne            : liste(d_matrix, coords_etendues, n_noeuds)
# ============================================================

library(geosphere)

# ------------------------------------------------------------
#  construire_matrice_distances()
#
#  Arguments :
#    candidats   : data.frame avec colonnes longitude, latitude
#    depot_lon   : longitude du dépôt
#    depot_lat   : latitude du dépôt
#    decharge_lon: longitude de la décharge
#    decharge_lat: latitude de la décharge
#
#  Retourne une liste :
#    $d_matrix        : matrice (n+2)x(n+2) de distances en km
#    $coords_etendues : data.frame de tous les nœuds avec leur rôle
#    $n               : nombre de points candidats
#    $idx_depot       : index du dépôt dans la matrice (toujours 1)
#    $idx_decharge    : index de la décharge dans la matrice (toujours n+2)
# ------------------------------------------------------------

construire_matrice_distances <- function(candidats,
                                         depot_lon, depot_lat,
                                         decharge_lon, decharge_lat) {
    
    n <- nrow(candidats)
    
    # ── 1. Construire le tableau de tous les nœuds dans l'ordre ──
    # Ordre : [dépôt | candidat_1 ... candidat_n | décharge]
    # Index R : 1 = dépôt, 2..(n+1) = candidats, (n+2) = décharge
    
    coords_etendues <- data.frame(
        index     = c(1, seq(2, n + 1), n + 2),
        role      = c("depot", rep("candidat", n), "decharge"),
        id_local  = c(0, seq_len(n), n + 1),   # 0 = dépôt, 1..n = candidats, n+1 = décharge
        longitude = c(depot_lon,    candidats$longitude, decharge_lon),
        latitude  = c(depot_lat,    candidats$latitude,  decharge_lat),
        label     = c("Dépôt (O)", paste0("Point_", seq_len(n)), "Décharge (S)")
    )
    
    # ── 2. Matrice de distances Haversine (en mètres → convertir en km) ──
    coords_mat <- as.matrix(coords_etendues[, c("longitude", "latitude")])
    
    d_metres <- distm(coords_mat, coords_mat, fun = distHaversine)
    d_matrix <- d_metres / 1000   # → kilomètres
    
    # Diagonale à 0 (distance d'un nœud à lui-même)
    diag(d_matrix) <- 0
    
    # ── 3. Vérifications de cohérence ──
    stopifnot(
        nrow(d_matrix) == n + 2,
        ncol(d_matrix) == n + 2
    )
    
    # ── 4. Retour ──
    list(
        d_matrix        = d_matrix,
        coords_etendues = coords_etendues,
        n               = n,
        idx_depot       = 1L,
        idx_decharge    = as.integer(n + 2)
    )
}


# ------------------------------------------------------------
#  distance_arc()
#  Accès rapide à la distance entre deux nœuds du graphe étendu.
#
#  Arguments :
#    res_dist : liste retournée par construire_matrice_distances()
#    de       : index du nœud de départ  (1 = dépôt, n+2 = décharge)
#    vers     : index du nœud d'arrivée
#
#  Retourne : distance en km (scalaire)
# ------------------------------------------------------------

distance_arc <- function(res_dist, de, vers) {
    res_dist$d_matrix[de, vers]
}


# ------------------------------------------------------------
#  distance_tournee()
#  Calcule la distance totale d'une séquence de nœuds (tournée).
#
#  Arguments :
#    res_dist : liste retournée par construire_matrice_distances()
#    sequence : vecteur d'index de nœuds dans l'ordre de visite
#               ex : c(1, 3, 5, 2, 7)  (1 = dépôt, 7 = décharge)
#
#  Retourne : distance totale en km (scalaire)
# ------------------------------------------------------------

distance_tournee <- function(res_dist, sequence) {
    if (length(sequence) < 2) return(0)
    total <- 0
    for (k in seq_len(length(sequence) - 1)) {
        total <- total + distance_arc(res_dist, sequence[k], sequence[k + 1])
    }
    total
}


# ------------------------------------------------------------
#  resumer_matrice()
#  Affiche un résumé lisible de la matrice pour debug.
# ------------------------------------------------------------

resumer_matrice <- function(res_dist) {
    cat("=== Matrice de distances (km) ===\n")
    cat("Taille      :", nrow(res_dist$d_matrix), "x", ncol(res_dist$d_matrix), "\n")
    cat("Index dépôt :", res_dist$idx_depot, "\n")
    cat("Index décharge:", res_dist$idx_decharge, "\n")
    cat("Nœuds       :\n")
    print(res_dist$coords_etendues[, c("index", "role", "label", "longitude", "latitude")])
    cat("\nDistance dépôt → décharge (ligne droite) :",
        round(res_dist$d_matrix[res_dist$idx_depot, res_dist$idx_decharge], 2), "km\n")
    cat("Distance min entre candidats :",
        round(min(res_dist$d_matrix[2:(res_dist$n+1), 2:(res_dist$n+1)][
            res_dist$d_matrix[2:(res_dist$n+1), 2:(res_dist$n+1)] > 0]), 2), "km\n")
    cat("Distance max entre candidats :",
        round(max(res_dist$d_matrix[2:(res_dist$n+1), 2:(res_dist$n+1)]), 2), "km\n")
}