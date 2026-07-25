# ============================================================
#  modules/cache_graphe.R
#  Gestion du cache du graphe OSM et des matrices utilisateur
# ============================================================

library(digest)

# ── Constantes ───────────────────────────────────────────────
CACHE_DEFAUT <- "cache/"
DATA_DEFAUT <- "data/"

#' Charger le graphe OSM depuis le cache
#' 
#' @param dossier_cache Dossier contenant le graphe (par défaut "data/")
#' @return Le graphe chargé
#' @export
charger_graphe_cache <- function(dossier_cache = DATA_DEFAUT) {
    fichier_graphe <- file.path(dossier_cache, "graphe_dakar.Rds")
    
    if (file.exists(fichier_graphe)) {
        message("[Cache] ✅ Chargement du graphe depuis : ", fichier_graphe)
        return(readRDS(fichier_graphe))
    } else {
        # Essayer de charger depuis le dossier courant
        fichier_alt <- "graphe_dakar.Rds"
        if (file.exists(fichier_alt)) {
            message("[Cache] ✅ Chargement du graphe depuis : ", fichier_alt)
            return(readRDS(fichier_alt))
        }
        
        # Essayer depuis utils/
        fichier_utils <- "utils/graphe_dakar.Rds"
        if (file.exists(fichier_utils)) {
            message("[Cache] ✅ Chargement du graphe depuis : ", fichier_utils)
            return(readRDS(fichier_utils))
        }
        
        stop("❌ Graphe non trouvé. Exécutez d'abord scripts/preparer_graphe_dakar.R")
    }
}

#' Vérifier si le graphe est en cache
#' 
#' @param dossier_cache Dossier contenant le graphe (par défaut "data/")
#' @return TRUE si le graphe existe, FALSE sinon
#' @export
graphe_en_cache <- function(dossier_cache = DATA_DEFAUT) {
    # Vérifier dans data/
    fichier_graphe <- file.path(dossier_cache, "graphe_dakar.Rds")
    if (file.exists(fichier_graphe)) return(TRUE)
    
    # Vérifier dans le dossier courant
    if (file.exists("graphe_dakar.Rds")) return(TRUE)
    
    # Vérifier dans utils/
    if (file.exists("utils/graphe_dakar.Rds")) return(TRUE)
    
    return(FALSE)
}

#' Générer une clé unique pour le cache utilisateur
#' 
#' La clé est stable et reproductible pour une même configuration
#' 
#' @param candidats Data.frame des points candidats
#' @param depot_lon Longitude du dépôt
#' @param depot_lat Latitude du dépôt
#' @param decharge_lon Longitude de la décharge
#' @param decharge_lat Latitude de la décharge
#' @return Clé MD5 unique et stable
#' @export
cle_cache_utilisateur <- function(candidats, depot_lon, depot_lat,
                                  decharge_lon, decharge_lat) {
    
    # Vérifier que les coordonnées sont valides
    if (is.null(depot_lon) || is.na(depot_lon)) stop("depot_lon est NULL ou NA")
    if (is.null(depot_lat) || is.na(depot_lat)) stop("depot_lat est NULL ou NA")
    if (is.null(decharge_lon) || is.na(decharge_lon)) stop("decharge_lon est NULL ou NA")
    if (is.null(decharge_lat) || is.na(decharge_lat)) stop("decharge_lat est NULL ou NA")
    
    # Trier les candidats par ID pour garantir l'ordre
    if (!is.null(candidats$id)) {
        candidats <- candidats[order(candidats$id), ]
    } else {
        # Si pas d'ID, trier par longitude puis latitude
        candidats <- candidats[order(candidats$longitude, candidats$latitude), ]
    }
    
    # Arrondir à 6 décimales (précision ~ 0.1 m)
    lon_arr <- round(candidats$longitude, 6)
    lat_arr <- round(candidats$latitude, 6)
    
    # Créer une chaîne stable avec séparateur
    chaine <- paste(
        nrow(candidats),
        paste(lon_arr, collapse = "_"),
        paste(lat_arr, collapse = "_"),
        round(depot_lon, 6),
        round(depot_lat, 6),
        round(decharge_lon, 6),
        round(decharge_lat, 6),
        sep = "|||"
    )
    
    # Générer le hash
    digest::digest(chaine, algo = "md5")
}

#' Sauvegarder la matrice de distances calculée (cache utilisateur)
#' 
#' @param res_dist Résultat de matrice_distances_osm()
#' @param candidats Data.frame des points candidats
#' @param depot_lon Longitude du dépôt
#' @param depot_lat Latitude du dépôt
#' @param decharge_lon Longitude de la décharge
#' @param decharge_lat Latitude de la décharge
#' @param dossier_cache Dossier de cache (par défaut "cache/")
#' @return Chemin du fichier sauvegardé
#' @export
sauvegarder_matrice_utilisateur <- function(res_dist, candidats, 
                                            depot_lon, depot_lat,
                                            decharge_lon, decharge_lat,
                                            dossier_cache = CACHE_DEFAUT) {
    
    # ── 1. Créer le dossier si nécessaire ─────────────────────
    if (!dir.exists(dossier_cache)) {
        dir.create(dossier_cache, recursive = TRUE)
        message("[Cache] 📁 Dossier créé : ", dossier_cache)
    }
    
    # ── 2. Générer la clé unique ──────────────────────────────
    cle <- cle_cache_utilisateur(
        candidats = candidats,
        depot_lon = depot_lon,
        depot_lat = depot_lat,
        decharge_lon = decharge_lon,
        decharge_lat = decharge_lat
    )
    
    fichier <- file.path(dossier_cache, paste0("matrice_utilisateur_", cle, ".Rds"))
    
    # ── 3. S'assurer que coords_snap existe ────────────────────
    if (is.null(res_dist$coords_snap)) {
        message("[Cache] ⚠️ coords_snap manquant, création à partir de coords_etendues")
        res_dist$coords_snap <- res_dist$coords_etendues
        res_dist$coords_snap$snap_lon <- res_dist$coords_snap$longitude
        res_dist$coords_snap$snap_lat <- res_dist$coords_snap$latitude
        res_dist$coords_snap$node_id <- paste0("node_", res_dist$coords_snap$index)
    }
    
    # ── 4. S'assurer que le graphe est sauvegardé ──────────────
    if (is.null(res_dist$graphe)) {
        message("[Cache] ⚠️ graphe manquant, tentative de chargement...")
        if (graphe_en_cache()) {
            res_dist$graphe <- charger_graphe_cache()
            message("[Cache] ✅ Graphe chargé et ajouté à res_dist")
        } else {
            warning("[Cache] ❌ Impossible de sauvegarder le graphe")
        }
    }
    
    # ── 5. Ajouter des métadonnées ─────────────────────────────
    attr(res_dist, "cache_info") <- list(
        date_creation = Sys.time(),
        nb_candidats = nrow(candidats),
        depot = c(depot_lon, depot_lat),
        decharge = c(decharge_lon, decharge_lat),
        source = res_dist$source %||% "unknown",
        has_coords_snap = "coords_snap" %in% names(res_dist),
        has_graphe = "graphe" %in% names(res_dist),
        n_coords_snap = ifelse("coords_snap" %in% names(res_dist), 
                               nrow(res_dist$coords_snap), 0)
    )
    
    # ── 6. Sauvegarder ─────────────────────────────────────────
    saveRDS(res_dist, file = fichier)
    
    # ── 7. Résumé ──────────────────────────────────────────────
    message("[Cache] 💾 Matrice sauvegardée : ", fichier)
    message("[Cache]   - Clé : ", substr(cle, 1, 12), "...")
    message("[Cache]   - coords_snap : ", 
            ifelse("coords_snap" %in% names(res_dist), 
                   paste0("✅ ", nrow(res_dist$coords_snap), " lignes"), 
                   "❌"))
    message("[Cache]   - graphe : ", 
            ifelse("graphe" %in% names(res_dist), "✅", "❌"))
    message("[Cache]   - Source : ", res_dist$source %||% "inconnue")
    
    return(fichier)
}

#' Charger une matrice utilisateur depuis le cache
#' 
#' @param candidats Data.frame des points candidats
#' @param depot_lon Longitude du dépôt
#' @param depot_lat Latitude du dépôt
#' @param decharge_lon Longitude de la décharge
#' @param decharge_lat Latitude de la décharge
#' @param dossier_cache Dossier de cache (par défaut "cache/")
#' @param verifier_complet Vérifier que coords_snap et graphe sont présents
#' @return La matrice chargée ou NULL si non trouvée ou incomplète
#' @export
charger_matrice_utilisateur <- function(candidats, depot_lon, depot_lat,
                                        decharge_lon, decharge_lat,
                                        dossier_cache = CACHE_DEFAUT,
                                        verifier_complet = TRUE) {
    
    # ── 1. Vérifier que le dossier existe ──────────────────────
    if (!dir.exists(dossier_cache)) {
        message("[Cache] 📁 Dossier cache inexistant : ", dossier_cache)
        return(NULL)
    }
    
    # ── 2. Générer la clé unique ──────────────────────────────
    cle <- tryCatch({
        cle_cache_utilisateur(
            candidats = candidats,
            depot_lon = depot_lon,
            depot_lat = depot_lat,
            decharge_lon = decharge_lon,
            decharge_lat = decharge_lat
        )
    }, error = function(e) {
        message("[Cache] ❌ Erreur génération clé : ", e$message)
        return(NULL)
    })
    
    if (is.null(cle)) return(NULL)
    
    fichier <- file.path(dossier_cache, paste0("matrice_utilisateur_", cle, ".Rds"))
    
    # ── 3. Vérifier que le fichier existe ──────────────────────
    if (!file.exists(fichier)) {
        message("[Cache] ℹ️ Aucun cache trouvé pour cette configuration")
        return(NULL)
    }
    
    # ── 4. Charger le fichier ──────────────────────────────────
    res_dist <- readRDS(fichier)
    
    # ── 5. Vérifier l'intégrité ────────────────────────────────
    if (verifier_complet) {
        # Vérifier coords_snap
        if (is.null(res_dist$coords_snap)) {
            message("[Cache] ⚠️ coords_snap manquant, tentative de réparation...")
            if (!is.null(res_dist$coords_etendues)) {
                res_dist$coords_snap <- res_dist$coords_etendues
                res_dist$coords_snap$snap_lon <- res_dist$coords_snap$longitude
                res_dist$coords_snap$snap_lat <- res_dist$coords_snap$latitude
                res_dist$coords_snap$node_id <- paste0("node_", res_dist$coords_snap$index)
                message("[Cache] ✅ coords_snap réparé")
            } else {
                message("[Cache] ❌ Impossible de réparer coords_snap")
                return(NULL)
            }
        }
        
        # Vérifier le graphe
        if (is.null(res_dist$graphe)) {
            message("[Cache] ⚠️ graphe manquant, tentative de chargement...")
            if (graphe_en_cache()) {
                res_dist$graphe <- charger_graphe_cache()
                message("[Cache] ✅ graphe chargé")
            } else {
                message("[Cache] ❌ Impossible de charger le graphe")
                return(NULL)
            }
        }
    }
    
    # ── 6. Résumé ──────────────────────────────────────────────
    message("[Cache] 📂 Matrice utilisateur chargée : ", basename(fichier))
    message("[Cache]   - Clé : ", substr(cle, 1, 12), "...")
    message("[Cache]   - coords_snap : ", 
            ifelse("coords_snap" %in% names(res_dist), 
                   paste0("✅ ", nrow(res_dist$coords_snap), " lignes"), 
                   "❌"))
    message("[Cache]   - graphe : ", 
            ifelse("graphe" %in% names(res_dist), "✅", "❌"))
    message("[Cache]   - Source : ", res_dist$source %||% "inconnue")
    
    # Ajouter un attribut pour indiquer la source
    attr(res_dist, "source") <- "cache_utilisateur"
    
    return(res_dist)
}

#' Vider le cache utilisateur
#' 
#' @param dossier_cache Dossier de cache (par défaut "cache/")
#' @param pattern Motif de recherche (par défaut "matrice_utilisateur_.*\\.Rds$")
#' @param ancien_que Supprimer les fichiers plus anciens que (en jours)
#' @return Nombre de fichiers supprimés
#' @export
vider_cache_utilisateur <- function(dossier_cache = CACHE_DEFAUT, 
                                    pattern = "matrice_utilisateur_.*\\.Rds$",
                                    ancien_que = NULL) {
    
    if (!dir.exists(dossier_cache)) {
        message("[Cache] 📁 Dossier cache inexistant : ", dossier_cache)
        return(0)
    }
    
    fichiers <- list.files(dossier_cache, pattern = pattern, full.names = TRUE)
    
    if (length(fichiers) == 0) {
        message("[Cache] ℹ️ Aucun fichier à supprimer")
        return(0)
    }
    
    # Supprimer les fichiers trop anciens
    if (!is.null(ancien_que) && is.numeric(ancien_que)) {
        date_limite <- Sys.time() - (ancien_que * 24 * 60 * 60)
        info_fichiers <- file.info(fichiers)
        fichiers <- fichiers[info_fichiers$mtime < date_limite]
        
        if (length(fichiers) == 0) {
            message("[Cache] ℹ️ Aucun fichier plus ancien que ", ancien_que, " jours")
            return(0)
        }
    }
    
    # Supprimer les fichiers
    file.remove(fichiers)
    message("[Cache] 🗑️ ", length(fichiers), " fichiers supprimés du cache")
    
    # Afficher les fichiers supprimés
    if (length(fichiers) <= 10) {
        for (f in fichiers) {
            message("[Cache]   - ", basename(f))
        }
    } else {
        message("[Cache]   - ", paste(head(basename(fichiers), 5), collapse = ", "), " ...")
    }
    
    return(length(fichiers))
}

#' Nettoyer le cache (supprimer les doublons)
#' 
#' @param dossier_cache Dossier de cache (par défaut "cache/")
#' @return Nombre de fichiers supprimés
#' @export
nettoyer_cache <- function(dossier_cache = CACHE_DEFAUT) {
    if (!dir.exists(dossier_cache)) {
        message("[Cache] 📁 Dossier cache inexistant : ", dossier_cache)
        return(0)
    }
    
    fichiers <- list.files(dossier_cache, 
                           pattern = "matrice_utilisateur_.*\\.Rds$", 
                           full.names = TRUE)
    
    if (length(fichiers) == 0) {
        message("[Cache] ℹ️ Aucun fichier à nettoyer")
        return(0)
    }
    
    message("[Cache] 🧹 Nettoyage des doublons...")
    
    fichiers_supprimes <- 0
    contenus <- list()
    
    for (f in fichiers) {
        tryCatch({
            obj <- readRDS(f)
            # Créer un hash du contenu (matrice + métadonnées)
            hash_contenu <- digest::digest(
                list(
                    dim = dim(obj$d_matrix),
                    values = round(obj$d_matrix, 6),
                    n = obj$n,
                    source = obj$source %||% "unknown"
                ),
                algo = "md5"
            )
            
            if (hash_contenu %in% names(contenus)) {
                # Doublon trouvé
                file.remove(f)
                fichiers_supprimes <- fichiers_supprimes + 1
                message("[Cache] 🗑️ Doublon supprimé : ", basename(f))
            } else {
                contenus[[hash_contenu]] <- f
            }
        }, error = function(e) {
            # Ignorer les fichiers corrompus
            message("[Cache] ⚠️ Fichier corrompu : ", basename(f))
        })
    }
    
    message("[Cache] 🧹 Nettoyage terminé : ", fichiers_supprimes, " doublons supprimés")
    return(fichiers_supprimes)
}

#' Vider TOUT le cache
#' 
#' @param dossier_cache Dossier de cache (par défaut "cache/")
#' @param dossier_data Dossier des données (par défaut "data/")
#' @param tout Supprimer aussi le graphe ? (défaut: FALSE)
#' @return Nombre de fichiers supprimés
#' @export
vider_cache_complet <- function(dossier_cache = CACHE_DEFAUT, 
                                dossier_data = DATA_DEFAUT,
                                tout = FALSE) {
    total <- 0
    
    # Vider le cache utilisateur
    total <- total + vider_cache_utilisateur(dossier_cache)
    
    # Supprimer le graphe si demandé
    if (tout) {
        # Dans data/
        fichier_graphe <- file.path(dossier_data, "graphe_dakar.Rds")
        if (file.exists(fichier_graphe)) {
            file.remove(fichier_graphe)
            message("[Cache] 🗑️ Graphe supprimé : ", fichier_graphe)
            total <- total + 1
        }
        
        # Dans le dossier courant
        if (file.exists("graphe_dakar.Rds")) {
            file.remove("graphe_dakar.Rds")
            message("[Cache] 🗑️ Graphe supprimé : graphe_dakar.Rds")
            total <- total + 1
        }
        
        # Dans utils/
        if (file.exists("utils/graphe_dakar.Rds")) {
            file.remove("utils/graphe_dakar.Rds")
            message("[Cache] 🗑️ Graphe supprimé : utils/graphe_dakar.Rds")
            total <- total + 1
        }
    }
    
    message("[Cache] 🗑️ Total : ", total, " fichiers supprimés")
    return(total)
}

#' Obtenir des informations sur le cache
#' 
#' @param dossier_cache Dossier de cache (par défaut "cache/")
#' @return Une liste avec les informations du cache
#' @export
info_cache <- function(dossier_cache = CACHE_DEFAUT) {
    info <- list(
        dossier = dossier_cache,
        existe = dir.exists(dossier_cache),
        fichiers = list(),
        taille = 0,
        nb_fichiers = 0,
        details = list(),
        complet = FALSE
    )
    
    if (info$existe) {
        fichiers <- list.files(dossier_cache, pattern = "matrice_utilisateur_.*\\.Rds$", 
                               full.names = TRUE)
        info$nb_fichiers <- length(fichiers)
        info$taille <- sum(file.size(fichiers)) / 1024 / 1024  # en Mo
        info$fichiers <- basename(fichiers)
        
        # Analyser chaque fichier
        for (f in fichiers) {
            tryCatch({
                obj <- readRDS(f)
                info$details[[basename(f)]] <- list(
                    source = obj$source %||% "inconnue",
                    n = obj$n,
                    dim = paste(dim(obj$d_matrix), collapse = "x"),
                    has_coords_snap = "coords_snap" %in% names(obj),
                    has_graphe = "graphe" %in% names(obj),
                    cache_info = attr(obj, "cache_info")
                )
            }, error = function(e) {
                info$details[[basename(f)]] <- list(error = e$message)
            })
        }
        
        # Vérifier si au moins un cache est complet
        complet <- sapply(info$details, function(x) {
            !is.null(x$has_coords_snap) && x$has_coords_snap &&
                !is.null(x$has_graphe) && x$has_graphe
        })
        info$complet <- any(complet)
    }
    
    return(info)
}

#' Vérifier l'intégrité d'un fichier de cache
#' 
#' @param fichier Chemin du fichier de cache
#' @return TRUE si le fichier est valide, FALSE sinon
#' @export
verifier_cache_fichier <- function(fichier) {
    if (!file.exists(fichier)) {
        message("[Cache] ❌ Fichier inexistant : ", fichier)
        return(FALSE)
    }
    
    tryCatch({
        obj <- readRDS(fichier)
        
        # Vérifier les composants requis
        requis <- c("d_matrix", "coords_etendues", "n", "idx_depot", "idx_decharge")
        manquants <- setdiff(requis, names(obj))
        if (length(manquants) > 0) {
            message("[Cache] ❌ Composants manquants : ", paste(manquants, collapse = ", "))
            return(FALSE)
        }
        
        # Vérifier que d_matrix est une matrice
        if (!is.matrix(obj$d_matrix)) {
            message("[Cache] ❌ d_matrix n'est pas une matrice")
            return(FALSE)
        }
        
        # Vérifier coords_snap (optionnel mais recommandé)
        if (is.null(obj$coords_snap)) {
            message("[Cache] ⚠️ coords_snap manquant (recommandé)")
        }
        
        # Vérifier graphe (optionnel mais recommandé)
        if (is.null(obj$graphe)) {
            message("[Cache] ⚠️ graphe manquant (recommandé)")
        }
        
        message("[Cache] ✅ Fichier valide : ", basename(fichier))
        message("[Cache]   - Dimensions : ", nrow(obj$d_matrix), "x", ncol(obj$d_matrix))
        message("[Cache]   - Source : ", obj$source %||% "inconnue")
        message("[Cache]   - coords_snap : ", ifelse("coords_snap" %in% names(obj), "✅", "❌"))
        message("[Cache]   - graphe : ", ifelse("graphe" %in% names(obj), "✅", "❌"))
        
        return(TRUE)
        
    }, error = function(e) {
        message("[Cache] ❌ Erreur : ", e$message)
        return(FALSE)
    })
}

# Opérateur %||% pour les valeurs NULL
`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}