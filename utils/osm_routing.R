# ============================================================
#  utils/osm_routing.R
#  Routage sur le réseau routier OSM via dodgr
#  
#  ⚠️ C'est la méthode PRINCIPALE pour les distances réelles
# ============================================================

library(dodgr)
library(osmdata)
library(sf)
library(dplyr)

# Cache du graphe OSM avec fichier persistant
charger_graphe_osm <- function(
        bbox      = c(-17.55, 14.60, -17.30, 14.80),
        type_voie = "motorcar",
        forcer    = FALSE) {
    
    if (!dir.exists("utils")) dir.create("utils")
    
    cle_cache <- paste0(type_voie, "_", paste(round(bbox, 3), collapse = "_"))
    fichier_cache <- paste0("utils/graphe_osm_", cle_cache, ".Rds")
    
    if (!forcer && file.exists(fichier_cache)) {
        message("[OSM] Graphe chargé depuis le cache : ", fichier_cache)
        return(readRDS(fichier_cache))
    }
    
    message("[OSM] Téléchargement du réseau routier OSM (Dakar)...")
    
    bbox_etroit <- c(-17.50, 14.63, -17.32, 14.77)
    
    osm_raw <- opq(bbox = bbox_etroit) %>%
        add_osm_feature(key = "highway") %>%
        osmdata_sf()
    
    message("[OSM] Construction du graphe dodgr...")
    graphe <- weight_streetnet(
        osm_raw,
        wt_profile = "motorcar",
        type_col   = "highway"
    )
    
    # Nettoyer
    graphe <- graphe %>% filter(!is.na(d) & d > 0 & d < 50000)
    
    saveRDS(graphe, file = fichier_cache)
    message("[OSM] Graphe sauvegardé : ", fichier_cache, " (", nrow(graphe), " arêtes)")
    
    graphe
}

#' Snap des points sur le réseau routier (version CORRIGÉE)
snap_candidats_reseau <- function(graphe, coords_ext) {
    
    # Récupérer les sommets du graphe
    verts <- dodgr_vertices(graphe)
    
    # Extraire les coordonnées des points à snapper
    pts <- as.matrix(coords_ext[, c("longitude", "latitude")])
    
    # Pour chaque point, trouver le nœud le plus proche
    noeuds_snap <- apply(pts, 1, function(pt) {
        # Distance euclidienne au carré (plus rapide)
        dists <- (verts$x - pt[1])^2 + (verts$y - pt[2])^2
        idx_min <- which.min(dists)
        verts$id[idx_min]
    })
    
    # Récupérer les coordonnées des nœuds snappés
    snap_info <- verts[match(noeuds_snap, verts$id), c("id", "x", "y")]
    names(snap_info) <- c("node_id", "snap_lon", "snap_lat")
    
    # Combiner avec les coordonnées originales
    result <- cbind(coords_ext, snap_info)
    
    # Afficher les points snappés pour vérification
    message("[OSM] Points projetés :")
    for (i in 1:nrow(result)) {
        message("  ", result$label[i], " : ",
                round(result$longitude[i], 5), ",", round(result$latitude[i], 5),
                " → ",
                round(result$snap_lon[i], 5), ",", round(result$snap_lat[i], 5))
    }
    
    return(result)
}

#' Matrice de distances ROUTIÈRES OSM (version CORRIGÉE)
matrice_distances_osm <- function(candidats,
                                  depot_lon, depot_lat,
                                  decharge_lon, decharge_lat,
                                  graphe = NULL,
                                  bbox   = c(-17.55, 14.60, -17.30, 14.80)) {
    
    if (is.null(graphe)) {
        graphe <- charger_graphe_osm(bbox = bbox)
    }
    
    n <- nrow(candidats)
    
    # Construire les nœuds étendus
    coords_etendues <- data.frame(
        index     = c(1L, seq(2L, n + 1L), n + 2L),
        role      = c("depot", rep("candidat", n), "decharge"),
        id_local  = c(0L, seq_len(n), n + 1L),
        longitude = c(depot_lon,       candidats$longitude, decharge_lon),
        latitude  = c(depot_lat,       candidats$latitude,  decharge_lat),
        label     = c("Dépôt (O)", paste0("Point_", seq_len(n)), "Décharge (S)")
    )
    
    message("[OSM] 📍 Projection des points sur le réseau routier...")
    coords_snap <- snap_candidats_reseau(graphe, coords_etendues)
    
    message("[OSM] 🗺️ Calcul de la matrice de distances ROUTIÈRES (", n + 2, " x ", n + 2, ")...")
    
    pts_from <- as.matrix(coords_snap[, c("snap_lon", "snap_lat")])
    pts_to <- pts_from
    
    # Calculer les distances sur le réseau
    d_metres <- tryCatch({
        dodgr_dists(graphe, from = pts_from, to = pts_to)
    }, error = function(e) {
        message("[OSM] ⚠️ Erreur dodgr_dists : ", e$message)
        message("[OSM] 🔄 Fallback vers Haversine")
        geosphere::distm(pts_from, pts_to, fun = geosphere::distHaversine)
    })
    
    # Remplacer les NA par Haversine (plus réaliste)
    nb_na <- sum(is.na(d_metres))
    if (nb_na > 0) {
        message("[OSM] ⚠️ ", nb_na, " paires de nœuds non connectées")
        d_metres_hav <- geosphere::distm(pts_from, pts_to, fun = geosphere::distHaversine)
        d_metres[is.na(d_metres)] <- d_metres_hav[is.na(d_metres)]
    }
    
    d_matrix <- d_metres / 1000  # Conversion en km
    diag(d_matrix) <- 0
    
    # Afficher un résumé
    message("[OSM] ✅ Matrice OSM prête.")
    message("[OSM] Distance min : ", round(min(d_matrix[d_matrix > 0]), 3), " km")
    message("[OSM] Distance max : ", round(max(d_matrix), 3), " km")
    
    # Afficher la matrice pour vérification
    message("[OSM] Aperçu de la matrice (km) :")
    print(round(d_matrix, 2))
    
    list(
        d_matrix        = d_matrix,
        coords_etendues = coords_etendues,
        coords_snap     = coords_snap,
        n               = n,
        idx_depot       = 1L,
        idx_decharge    = as.integer(n + 2),
        graphe          = graphe,
        source          = "osm"
    )
}

#' ═══════════════════════════════════════════════════════════════
#'  Géométries des tournées OSM - VERSION ULTIME ROBUSTE
#' ═══════════════════════════════════════════════════════════════
geometries_tournees_osm <- function(sol_z, res_dist) {
    
    graphe <- res_dist$graphe
    coords_snap <- res_dist$coords_snap
    
    message("[OSM] geometries_tournees_osm appelé")
    message("[OSM] Nombre d'arcs : ", nrow(sol_z))
    
    if (is.null(graphe) || is.null(coords_snap) || nrow(sol_z) == 0) {
        message("[OSM] ⚠️ Données manquantes")
        return(NULL)
    }
    
    resultats <- list()
    verts <- dodgr_vertices(graphe)
    
    for (k in seq_len(nrow(sol_z))) {
        v_id <- sol_z$v[k]
        i <- sol_z$i[k]
        j <- sol_z$j[k]
        
        # Récupérer les points snappés
        pt_i <- coords_snap[coords_snap$index == i, ]
        pt_j <- coords_snap[coords_snap$index == j, ]
        
        # Si les points ne sont pas snappés, utiliser les coordonnées originales
        if (nrow(pt_i) == 0 || is.na(pt_i$snap_lon[1])) {
            pt_i <- data.frame(
                index = i,
                snap_lon = coords_snap$longitude[coords_snap$index == i],
                snap_lat = coords_snap$latitude[coords_snap$index == i],
                label = coords_snap$label[coords_snap$index == i]
            )
        }
        if (nrow(pt_j) == 0 || is.na(pt_j$snap_lon[1])) {
            pt_j <- data.frame(
                index = j,
                snap_lon = coords_snap$longitude[coords_snap$index == j],
                snap_lat = coords_snap$latitude[coords_snap$index == j],
                label = coords_snap$label[coords_snap$index == j]
            )
        }
        
        from_pts <- as.matrix(pt_i[, c("snap_lon", "snap_lat")])
        to_pts <- as.matrix(pt_j[, c("snap_lon", "snap_lat")])
        
        message("[OSM] Arc ", k, " : véhicule ", v_id, " ", i, " → ", j)
        message("[OSM]   from : ", paste(round(from_pts, 5), collapse = ", "))
        message("[OSM]   to   : ", paste(round(to_pts, 5), collapse = ", "))
        
        # Essayer plusieurs méthodes pour obtenir le chemin
        geom <- NULL
        
        # Méthode 1 : dodgr_paths direct
        chemin <- tryCatch({
            dodgr_paths(graphe, from = from_pts, to = to_pts)
        }, error = function(e) {
            NULL
        })
        
        if (!is.null(chemin) && length(chemin[[1]][[1]]) > 1) {
            ids_path <- chemin[[1]][[1]]
            coords_path <- verts[match(ids_path, verts$id), c("x", "y")]
            coords_path <- coords_path[!is.na(coords_path$x), ]
            
            if (nrow(coords_path) > 1) {
                geom <- st_linestring(as.matrix(coords_path))
                message("[OSM] ✅ Chemin trouvé (", nrow(coords_path), " points)")
            }
        }
        
        # Méthode 2 : Si échec, essayer de trouver un chemin entre les nœuds les plus proches
        if (is.null(geom)) {
            message("[OSM] 🔄 Tentative avec les nœuds les plus proches")
            
            # Trouver les nœuds les plus proches des points
            idx_from <- which.min((verts$x - from_pts[1])^2 + (verts$y - from_pts[2])^2)
            idx_to <- which.min((verts$x - to_pts[1])^2 + (verts$y - to_pts[2])^2)
            
            from_pts2 <- as.matrix(verts[idx_from, c("x", "y")])
            to_pts2 <- as.matrix(verts[idx_to, c("x", "y")])
            
            chemin2 <- tryCatch({
                dodgr_paths(graphe, from = from_pts2, to = to_pts2)
            }, error = function(e) {
                NULL
            })
            
            if (!is.null(chemin2) && length(chemin2[[1]][[1]]) > 1) {
                ids_path <- chemin2[[1]][[1]]
                coords_path <- verts[match(ids_path, verts$id), c("x", "y")]
                coords_path <- coords_path[!is.na(coords_path$x), ]
                
                if (nrow(coords_path) > 1) {
                    geom <- st_linestring(as.matrix(coords_path))
                    message("[OSM] ✅ Chemin trouvé avec nœuds les plus proches (", nrow(coords_path), " points)")
                }
            }
        }
        
        # Fallback final : ligne droite
        if (is.null(geom)) {
            geom <- st_linestring(matrix(
                c(pt_i$snap_lon[1], pt_i$snap_lat[1],
                  pt_j$snap_lon[1], pt_j$snap_lat[1]),
                ncol = 2, byrow = TRUE
            ))
            message("[OSM] 📏 Fallback : ligne droite")
        }
        
        resultats[[k]] <- list(
            v = v_id,
            i = i,
            j = j,
            label_i = pt_i$label[1],
            label_j = pt_j$label[1],
            distance_km = round(res_dist$d_matrix[i, j], 3),
            geometry = geom
        )
    }
    
    message("[OSM] geometries_tournees_osm terminé : ", length(resultats), " géométries")
    return(resultats)
}

#' Tracer les tournées OSM sur une carte Leaflet
tracer_tournees_osm <- function(m, geometries, couleurs) {
    
    COLORS <- c("#378ADD", "#D85A30", "#1D9E75", "#7F77DD", "#BA7517",
                "#E41A1C", "#FF7F00", "#4DAF4A", "#984EA3", "#A65628")
    
    for (geo in geometries) {
        if (is.null(geo)) next
        
        v_id <- geo$v
        couleur <- COLORS[(v_id - 1) %% length(COLORS) + 1]
        
        tryCatch({
            coords <- st_coordinates(geo$geometry)
            
            if (nrow(coords) > 1) {
                m <- m %>% addPolylines(
                    lng = coords[, 1],
                    lat = coords[, 2],
                    color = couleur,
                    weight = 4,
                    opacity = 0.9,
                    label = paste0("Véhicule ", v_id, " : ",
                                   geo$label_i, " → ", geo$label_j,
                                   " (", geo$distance_km, " km)")
                )
            }
        }, error = function(e) {
            message("[OSM] ⚠️ Erreur tracé : ", e$message)
        })
    }
    m
}