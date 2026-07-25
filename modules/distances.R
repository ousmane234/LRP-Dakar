# ============================================================
#  modules/distances.R
#  Calcul des distances ROUTIÈRES via OSM
#  
#  ⚠️ CE FICHIER UTILISE OSM PAR DÉFAUT
#  Haversine n'est utilisé qu'en dernier recours
# ============================================================

library(geosphere)

#' Construire la matrice de distances avec OSM (PRINCIPAL)
#' 
#' Cette fonction utilise OSM pour calculer les distances réelles.
#' Haversine n'est utilisé qu'en fallback ABSOLU.
#' 
#' @param candidats Data.frame des points candidats
#' @param depot_lon Longitude du dépôt
#' @param depot_lat Latitude du dépôt
#' @param decharge_lon Longitude de la décharge
#' @param decharge_lat Latitude de la décharge
#' @param graphe Graphe OSM (obligatoire)
#' @param utiliser_osm TRUE pour utiliser OSM (défaut), FALSE pour forcer Haversine
#' @return Liste avec d_matrix, coords_etendues, n, idx_depot, idx_decharge
#' @export
construire_matrice_distances <- function(candidats,
                                         depot_lon, depot_lat,
                                         decharge_lon, decharge_lat,
                                         graphe = NULL,
                                         utiliser_osm = TRUE) {
    
    # ── Vérifications ──────────────────────────────────────────
    if (is.null(candidats) || nrow(candidats) == 0) {
        stop("❌ Aucun candidat fourni")
    }
    
    cols_requises <- c("longitude", "latitude")
    manquantes <- setdiff(cols_requises, names(candidats))
    if (length(manquantes) > 0) {
        stop("❌ Colonnes manquantes dans candidats : ", 
             paste(manquantes, collapse = ", "))
    }
    
    # Vérifier les coordonnées
    if (is.na(depot_lon) || is.na(depot_lat)) {
        stop("❌ Coordonnées du dépôt invalides")
    }
    if (is.na(decharge_lon) || is.na(decharge_lat)) {
        stop("❌ Coordonnées de la décharge invalides")
    }
    
    n <- nrow(candidats)
    message("[Distances] 🗺️ Calcul des distances pour ", n, " points")
    
    # ── Construction des coordonnées étendues ──────────────────
    coords_etendues <- data.frame(
        index     = c(1L, seq(2L, n + 1L), n + 2L),
        role      = c("depot", rep("candidat", n), "decharge"),
        id_local  = c(0L, seq_len(n), n + 1L),
        longitude = c(depot_lon,    candidats$longitude, decharge_lon),
        latitude  = c(depot_lat,    candidats$latitude,  decharge_lat),
        label     = c("Dépôt (O)", paste0("Point_", seq_len(n)), "Décharge (S)"),
        stringsAsFactors = FALSE
    )
    
    # ── Calcul avec OSM (PRINCIPAL) ────────────────────────────
    d_matrix <- NULL
    
    if (utiliser_osm && !is.null(graphe)) {
        message("[Distances] 🚗 Utilisation du réseau OSM pour les distances réelles...")
        
        tryCatch({
            # Utiliser la fonction OSM
            res_osm <- matrice_distances_osm(
                candidats = candidats,
                depot_lon = depot_lon,
                depot_lat = depot_lat,
                decharge_lon = decharge_lon,
                decharge_lat = decharge_lat,
                graphe = graphe
            )
            d_matrix <- res_osm$d_matrix
            message("[Distances] ✅ Matrice OSM calculée avec succès")
            
        }, error = function(e) {
            message("[Distances] ❌ Erreur OSM : ", e$message)
            message("[Distances] ⚠️ Fallback vers Haversine (à éviter !)")
            utiliser_osm <- FALSE  # Forcer le fallback
        })
    }
    
    # ── Fallback Haversine (CAS EXTRÊME) ──────────────────────
    if (is.null(d_matrix)) {
        message("[Distances] ⚠️⚠️⚠️ UTILISATION DE HAVERSINE (VOL D'OISEAU) ⚠️⚠️⚠️")
        message("[Distances] 📏 Ce n'est PAS réaliste pour un problème de collecte !")
        message("[Distances] 📏 Les distances seront sous-estimées !")
        
        coords_mat <- as.matrix(coords_etendues[, c("longitude", "latitude")])
        
        # Vérifier les NA
        if (any(is.na(coords_mat))) {
            warning("[Distances] ⚠️ Des coordonnées sont NA, correction...")
            coords_mat[is.na(coords_mat)] <- 0
        }
        
        d_metres <- tryCatch({
            geosphere::distm(coords_mat, coords_mat, fun = geosphere::distHaversine)
        }, error = function(e) {
            message("[Distances] ❌ Erreur Haversine : ", e$message)
            # Fallback ULTIME : distance euclidienne
            d <- as.matrix(dist(coords_mat))
            d * 111.32  # 1 degré ≈ 111.32 km
        })
        
        # Conversion en km
        d_matrix <- d_metres / 1000
        diag(d_matrix) <- 0
        d_matrix <- round(d_matrix, 6)
        
        # Ajouter un attribut pour indiquer que c'est Haversine
        attr(d_matrix, "source") <- "haversine_fallback"
        
        message("[Distances] ⚠️ Matrice Haversine utilisée (vol d'oiseau)")
    }
    
    # ── Vérification finale ────────────────────────────────────
    if (!is.matrix(d_matrix)) {
        stop("❌ d_matrix n'est pas une matrice")
    }
    
    if (nrow(d_matrix) != n + 2 || ncol(d_matrix) != n + 2) {
        stop("❌ Taille incorrecte : attendu ", n + 2, "x", n + 2)
    }
    
    message("[Distances] ✅ Matrice finale : ", nrow(d_matrix), "x", ncol(d_matrix))
    message("[Distances]   Distance min : ", round(min(d_matrix[d_matrix > 0]), 3), " km")
    message("[Distances]   Distance max : ", round(max(d_matrix), 3), " km")
    message("[Distances]   Source : ", attr(d_matrix, "source") %||% "OSM")
    
    # ── Retour ──────────────────────────────────────────────────
    list(
        d_matrix        = d_matrix,
        coords_etendues = coords_etendues,
        n               = n,
        idx_depot       = 1L,
        idx_decharge    = as.integer(n + 2),
        source          = attr(d_matrix, "source") %||% "osm"
    )
}

#' Alias pour OSM
#' @export
matrice_distances_principale <- construire_matrice_distances

#' Vérifier si une matrice est OSM ou Haversine
#' @export
est_matrice_osm <- function(res_dist) {
    if (is.null(res_dist)) return(FALSE)
    if (!is.null(res_dist$source)) {
        return(res_dist$source == "osm")
    }
    # Vérifier si la matrice a un attribut source
    if (!is.null(res_dist$d_matrix)) {
        src <- attr(res_dist$d_matrix, "source")
        return(!is.null(src) && src == "osm")
    }
    return(FALSE)
}

#' Avertir si on utilise Haversine
#' @export
verifier_source_distances <- function(res_dist) {
    if (is.null(res_dist)) {
        warning("⚠️ res_dist est NULL")
        return(FALSE)
    }
    
    if (!is.null(res_dist$source) && res_dist$source == "haversine_fallback") {
        warning("⚠️⚠️⚠️ UTILISATION DE HAVERSINE (VOL D'OISEAU) ⚠️⚠️⚠️")
        warning("Les distances ne sont PAS réalistes pour un problème de collecte !")
        return(FALSE)
    }
    
    if (!is.null(res_dist$source) && res_dist$source == "osm") {
        message("✅ Matrice OSM (distances réelles)")
        return(TRUE)
    }
    
    # Vérifier l'attribut
    if (!is.null(res_dist$d_matrix)) {
        src <- attr(res_dist$d_matrix, "source")
        if (!is.null(src) && src == "haversine_fallback") {
            warning("⚠️⚠️⚠️ UTILISATION DE HAVERSINE ⚠️⚠️⚠️")
            return(FALSE)
        }
    }
    
    message("ℹ️ Source inconnue, vérification manuelle requise")
    return(TRUE)
}

# ── Fonction utilitaire ──────────────────────────────────────

#' Opérateur %||% pour les valeurs NULL
`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}