# ============================================================
#  modules/cache_graphe.R
#  Gestion du cache du graphe OSM et des matrices utilisateur
# ============================================================

library(digest)

#' Charger le graphe OSM depuis le cache
charger_graphe_cache <- function(dossier_cache = "data/") {
    fichier_graphe <- paste0(dossier_cache, "graphe_dakar.Rds")
    
    if (file.exists(fichier_graphe)) {
        message("[Cache] Chargement du graphe depuis : ", fichier_graphe)
        return(readRDS(fichier_graphe))
    } else {
        stop("❌ Graphe non trouvé. Exécutez d'abord scripts/preparer_graphe_dakar.R")
    }
}

#' Vérifier si le graphe est en cache
graphe_en_cache <- function(dossier_cache = "data/") {
    fichier_graphe <- paste0(dossier_cache, "graphe_dakar.Rds")
    return(file.exists(fichier_graphe))
}

#' Sauvegarder la matrice de distances calculée (cache utilisateur)
sauvegarder_matrice_utilisateur <- function(res_dist, candidats, 
                                            depot_lon, depot_lat,
                                            decharge_lon, decharge_lat,
                                            dossier_cache = "cache/") {
    
    if (!dir.exists(dossier_cache)) dir.create(dossier_cache, recursive = TRUE)
    
    # Clé unique basée sur les coordonnées utilisateur
    cle <- digest::digest(paste(
        nrow(candidats),
        paste(round(candidats$longitude, 6), collapse = "_"),
        paste(round(candidats$latitude, 6), collapse = "_"),
        round(depot_lon, 6), round(depot_lat, 6),
        round(decharge_lon, 6), round(decharge_lat, 6)
    ), algo = "md5")
    
    fichier <- paste0(dossier_cache, "matrice_utilisateur_", cle, ".Rds")
    saveRDS(res_dist, file = fichier)
    message("[Cache] Matrice sauvegardée : ", fichier)
    return(fichier)
}

#' Charger une matrice utilisateur depuis le cache
charger_matrice_utilisateur <- function(candidats, depot_lon, depot_lat,
                                        decharge_lon, decharge_lat,
                                        dossier_cache = "cache/") {
    
    if (!dir.exists(dossier_cache)) return(NULL)
    
    cle <- digest::digest(paste(
        nrow(candidats),
        paste(round(candidats$longitude, 6), collapse = "_"),
        paste(round(candidats$latitude, 6), collapse = "_"),
        round(depot_lon, 6), round(depot_lat, 6),
        round(decharge_lon, 6), round(decharge_lat, 6)
    ), algo = "md5")
    
    fichier <- paste0(dossier_cache, "matrice_utilisateur_", cle, ".Rds")
    
    if (file.exists(fichier)) {
        message("[Cache] Matrice utilisateur chargée : ", fichier)
        res_dist <- readRDS(fichier)
        attr(res_dist, "source") <- "cache_utilisateur"
        return(res_dist)
    }
    return(NULL)
}

#' Vider le cache utilisateur
vider_cache_utilisateur <- function(dossier_cache = "cache/") {
    if (dir.exists(dossier_cache)) {
        fichiers <- list.files(dossier_cache, pattern = "matrice_utilisateur_.*\\.Rds$", full.names = TRUE)
        if (length(fichiers) > 0) {
            file.remove(fichiers)
            message("[Cache] ", length(fichiers), " fichiers supprimés")
            return(length(fichiers))
        }
    }
    return(0)
}