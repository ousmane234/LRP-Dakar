# ============================================================
#  utils/osm_routing.R
#  Routage sur le réseau routier OSM via dodgr
#
#  Remplace les distances Haversine (lignes droites) par des
#  distances et itinéraires réels sur les routes de Dakar.
#
#  Pourquoi dodgr et pas osrm ?
#    - dodgr est offline : pas de quota, pas de dépendance réseau
#    - Le graphe est construit une fois et réutilisé (cache)
#    - dodgr_paths() retourne les géométries pour tracer les routes
#
#  Fonctions principales :
#    charger_graphe_osm()        : télécharge et prépare le graphe routier
#    matrice_distances_osm()     : remplace distm() Haversine
#    geometries_tournees_osm()   : retourne les LINESTRING pour Leaflet
#    snap_candidats_reseau()     : projette les centroïdes sur le réseau
# ============================================================

library(dodgr)
library(osmdata)
library(sf)
library(dplyr)

# ------------------------------------------------------------
#  CACHE du graphe — évite de retélécharger à chaque run
# ------------------------------------------------------------
############.graphe_cache <- new.env(parent = emptyenv())


# ------------------------------------------------------------
#  charger_graphe_osm()
#  Télécharge le réseau routier OSM pour Dakar et construit
#  le graphe pondéré dodgr.
#
#  Arguments :
#    bbox      : vecteur c(xmin, ymin, xmax, ymax) en WGS84
#                Par défaut : emprise de la presqu'île de Dakar
#    type_voie : profil de pondération dodgr ("motorcar" ou "foot")
#    forcer    : si TRUE, ignore le cache et retélécharge
#
#  Retourne : graphe dodgr (data.frame pondéré)
# ------------------------------------------------------------

# charger_graphe_osm <- function(
#         bbox      = c(-17.55, 14.60, -17.30, 14.80),   # Dakar presqu'île
#         type_voie = "motorcar",
#         forcer    = FALSE) {
#     
#     cle_cache <- paste0(type_voie, "_", paste(round(bbox, 3), collapse = "_"))
#     
#     if (!forcer && exists(cle_cache, envir = .graphe_cache)) {
#         message("[OSM] Graphe chargé depuis le cache.")
#         return(get(cle_cache, envir = .graphe_cache))
#     }
#     
#     message("[OSM] Téléchargement du réseau routier OSM...")
#     
#     # Requête OSM : routes carrossables dans la bbox
#     osm_raw <- opq(bbox = bbox) %>%
#         add_osm_feature(key = "highway") %>%
#         osmdata_sf()
#     
#     # Construction du graphe dodgr
#     graphe <- weight_streetnet(
#         osm_raw,
#         wt_profile = type_voie,
#         type_col   = "highway"
#     )
#     
#     # Mise en cache
#     assign(cle_cache, graphe, envir = .graphe_cache)
#     message("[OSM] Graphe construit : ", nrow(graphe), " arêtes.")
#     
#     graphe
# }
# ------------------------------------------------------------
#  charger_graphe_osm() avec CACHE PHYSIQUE PERSISTANT (.Rds)
# ------------------------------------------------------------
charger_graphe_osm <- function(
        bbox      = c(-17.55, 14.60, -17.30, 14.80),   # Dakar presqu'île
        type_voie = "motorcar",
        forcer    = FALSE) {
    
    # 1. Définir un nom de fichier unique basé sur les paramètres
    cle_cache  <- paste0(type_voie, "_", paste(round(bbox, 3), collapse = "_"))
    fichier_cache <- paste0("utils/graphe_osm_", cle_cache, ".Rds")
    
    # 2. Si le fichier existe déjà sur le disque et qu'on ne force pas, on le lit directement
    if (!forcer && file.exists(fichier_cache)) {
        message("[OSM] Graphe trouvé sur le disque. Chargement instantané...")
        return(readRDS(fichier_cache))
    }
    
    # 3. Sinon, on procède au téléchargement lourd (une seule fois)
    message("[OSM] Fichier cache non trouvé. Téléchargement du réseau routier OSM (Dakar)...")
    
    osm_raw <- osmdata::opq(bbox = bbox) %>%
        osmdata::add_osm_feature(key = "highway") %>%
        osmdata::osmdata_sf()
    
    message("[OSM] Construction et pondération du graphe dodgr...")
    graphe <- dodgr::weight_streetnet(
        osm_raw,
        wt_profile = type_voie,
        type_col   = "highway"
    )
    
    # 4. On sauvegarde le résultat sur le disque pour les prochaines fois
    message("[OSM] Sauvegarde du graphe sur le disque : ", fichier_cache)
    saveRDS(graphe, file = fichier_cache)
    
    graphe
}



# ------------------------------------------------------------
#  snap_candidats_reseau()
#  Projette les points (centroïdes, dépôt, décharge) sur le
#  nœud du graphe routier le plus proche.
#
#  Arguments :
#    graphe     : graphe dodgr
#    coords_ext : data.frame retourné par construire_matrice_distances()
#                 → coords_etendues (longitude, latitude, index, label)
#
#  Retourne : data.frame avec colonnes supplémentaires
#             snap_lon, snap_lat, node_id (id du nœud dodgr snappé)
# ------------------------------------------------------------

snap_candidats_reseau <- function(graphe, coords_ext) {
    pts_mat <- as.matrix(coords_ext[, c("longitude", "latitude")])
    
    # dodgr_nearest_vertices : trouve le nœud le plus proche pour chaque point
    noeuds_snap <- dodgr_nearest_vertices(graphe, pts_mat)
    
    # Récupérer les coordonnées des nœuds snappés
    verts <- dodgr_vertices(graphe)
    
    snap_info <- verts[match(noeuds_snap, verts$id), c("id", "x", "y")]
    names(snap_info) <- c("node_id", "snap_lon", "snap_lat")
    
    cbind(coords_ext, snap_info)
}


# ------------------------------------------------------------
#  matrice_distances_osm()
#  Calcule la matrice de distances routières (km) entre tous
#  les nœuds du graphe étendu via dodgr.
#
#  Remplace construire_matrice_distances() (qui utilisait Haversine).
#  La structure retournée est identique pour compatibilité avec
#  resoudre_lrp().
#
#  Arguments :
#    candidats    : data.frame des points candidats (lon, lat)
#    depot_lon/lat, decharge_lon/lat : coordonnées dépôt et décharge
#    graphe       : graphe dodgr (si NULL, appelle charger_graphe_osm())
#    bbox         : bbox pour charger_graphe_osm si graphe=NULL
#
#  Retourne : même structure que construire_matrice_distances()
#    + $graphe            : le graphe dodgr (pour réutilisation)
#    + $coords_snap       : coords snappées sur le réseau
# ------------------------------------------------------------

matrice_distances_osm <- function(candidats,
                                  depot_lon, depot_lat,
                                  decharge_lon, decharge_lat,
                                  graphe = NULL,
                                  bbox   = c(-17.55, 14.60, -17.30, 14.80)) {
    
    # ── Graphe ────────────────────────────────────────────────
    if (is.null(graphe)) {
        graphe <- charger_graphe_osm(bbox = bbox)
    }
    
    n <- nrow(candidats)
    
    # ── Construire le tableau des nœuds étendus (comme distances.R) ──
    coords_etendues <- data.frame(
        index     = c(1L, seq(2L, n + 1L), n + 2L),
        role      = c("depot", rep("candidat", n), "decharge"),
        id_local  = c(0L, seq_len(n), n + 1L),
        longitude = c(depot_lon,       candidats$longitude, decharge_lon),
        latitude  = c(depot_lat,       candidats$latitude,  decharge_lat),
        label     = c("Dépôt (O)", paste0("Point_", seq_len(n)), "Décharge (S)")
    )
    
    # ── Snap sur le réseau routier ────────────────────────────
    message("[OSM] Projection des points sur le réseau routier...")
    coords_snap <- snap_candidats_reseau(graphe, coords_etendues)
    
    # ── Matrice de distances routières via dodgr ──────────────
    message("[OSM] Calcul de la matrice de distances routières (", n + 2, " x ", n + 2, ")...")
    
    pts_from <- as.matrix(coords_snap[, c("snap_lon", "snap_lat")])
    pts_to   <- pts_from   # matrice carrée : tous vers tous
    
    # dodgr_dists retourne une matrice en mètres
    d_metres <- dodgr_dists(graphe, from = pts_from, to = pts_to)
    
    # Remplacer NA (nœuds non connexes) par une grande valeur
    nb_na <- sum(is.na(d_metres))
    if (nb_na > 0) {
        warning("[OSM] ", nb_na, " paires de nœuds non connexes — remplacées par 9999 km.")
        d_metres[is.na(d_metres)] <- 9999 * 1000
    }
    
    d_matrix <- d_metres / 1000   # → km
    diag(d_matrix) <- 0
    
    message("[OSM] Matrice prête. Distance min : ", round(min(d_matrix[d_matrix > 0]), 2),
            " km | max : ", round(max(d_matrix), 2), " km")
    
    # ── Retour (même structure que construire_matrice_distances) ──
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


# ------------------------------------------------------------
#  geometries_tournees_osm()
#  Pour chaque arc actif z[v,i,j]=1, calcule le chemin routier
#  réel et retourne une liste de LINESTRING sf pour Leaflet.
#
#  Arguments :
#    sol_z       : data.frame des arcs actifs (v, i, j) — sol$sol_z
#    res_dist    : liste retournée par matrice_distances_osm()
#
#  Retourne : liste de data.frames avec
#             v, i, j, geometry (sf LINESTRING), distance_km
# ------------------------------------------------------------

geometries_tournees_osm <- function(sol_z, res_dist) {
    
    graphe      <- res_dist$graphe
    coords_snap <- res_dist$coords_snap
    
    resultats <- list()
    
    for (k in seq_len(nrow(sol_z))) {
        v_id <- sol_z$v[k]
        i    <- sol_z$i[k]
        j    <- sol_z$j[k]
        
        # Coordonnées snappées des deux nœuds
        pt_i <- coords_snap[coords_snap$index == i, c("snap_lon", "snap_lat")]
        pt_j <- coords_snap[coords_snap$index == j, c("snap_lon", "snap_lat")]
        
        # Chemin réel sur le graphe
        chemin <- tryCatch(
            dodgr_paths(graphe,
                        from = as.matrix(pt_i),
                        to   = as.matrix(pt_j)),
            error = function(e) NULL
        )
        
        if (!is.null(chemin) && length(chemin[[1]][[1]]) > 0) {
            # Récupérer les coordonnées des nœuds du chemin
            verts     <- dodgr_vertices(graphe)
            ids_path  <- chemin[[1]][[1]]
            coords_path <- verts[match(ids_path, verts$id), c("x", "y")]
            
            # Construire LINESTRING sf
            geom <- st_linestring(as.matrix(coords_path))
        } else {
            # Fallback : ligne droite si le chemin échoue
            geom <- st_linestring(matrix(
                c(pt_i$snap_lon, pt_i$snap_lat,
                  pt_j$snap_lon, pt_j$snap_lat),
                ncol = 2, byrow = TRUE
            ))
        }
        
        resultats[[k]] <- list(
            v           = v_id,
            i           = i,
            j           = j,
            label_i     = coords_snap$label[coords_snap$index == i],
            label_j     = coords_snap$label[coords_snap$index == j],
            distance_km = round(res_dist$d_matrix[i, j], 3),
            geometry    = geom
        )
    }
    
    resultats
}


# ------------------------------------------------------------
#  tracer_tournees_osm()
#  Ajoute les itinéraires OSM sur une carte Leaflet existante.
#  Remplace la boucle addPolylines dans server.R.
#
#  Arguments :
#    m           : objet leaflet existant
#    geometries  : liste retournée par geometries_tournees_osm()
#    couleurs    : vecteur de couleurs (une par véhicule)
#
#  Retourne : objet leaflet enrichi
# ------------------------------------------------------------

tracer_tournees_osm <- function(m, geometries, couleurs) {
    
    COLORS <- c("#378ADD", "#D85A30", "#1D9E75", "#7F77DD", "#BA7517",
                "#E41A1C", "#FF7F00", "#4DAF4A", "#984EA3", "#A65628")
    
    for (geo in geometries) {
        v_id    <- geo$v
        couleur <- COLORS[(v_id - 1) %% length(COLORS) + 1]
        coords  <- st_coordinates(geo$geometry)
        
        m <- m %>% addPolylines(
            lng     = coords[, 1],
            lat     = coords[, 2],
            color   = couleur,
            weight  = 3,
            opacity = 0.85,
            label   = paste0("Véhicule ", v_id, " : ",
                             geo$label_i, " → ", geo$label_j,
                             " (", geo$distance_km, " km)")
        )
    }
    m
}