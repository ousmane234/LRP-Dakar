# ============================================================
#  modules/chargement_matrice.R
#  Chargement de la matrice de distances depuis un fichier CSV
#  (généré par scripts/calcul_matrice_dakar.R)
# ============================================================

#' Charger la matrice de distances depuis un fichier CSV
#' 
#' @param fichier_csv Chemin vers le fichier CSV de la matrice
#' @param fichier_coords Chemin vers le fichier CSV des coordonnées (optionnel)
#' @param depot_lon Longitude du dépôt (optionnel, sinon utilise les labels)
#' @param depot_lat Latitude du dépôt (optionnel, sinon utilise les labels)
#' @param decharge_lon Longitude de la décharge (optionnel, sinon utilise les labels)
#' @param decharge_lat Latitude de la décharge (optionnel, sinon utilise les labels)
#' @param separateur Séparateur de colonnes (défaut: ",")
#' @return Liste avec la même structure que res_dist
#' @export
charger_matrice_csv <- function(fichier_csv, 
                                fichier_coords = NULL,
                                depot_lon = NULL, 
                                depot_lat = NULL,
                                decharge_lon = NULL, 
                                decharge_lat = NULL,
                                separateur = ",") {
    
    # ── 1. Vérifier que le fichier existe ──────────────────────
    if (!file.exists(fichier_csv)) {
        stop("❌ Fichier non trouvé : ", fichier_csv)
    }
    
    # ── 2. Charger la matrice ──────────────────────────────────
    df <- tryCatch({
        read.csv(fichier_csv, row.names = 1, check.names = FALSE, sep = separateur)
    }, error = function(e) {
        # Essayer avec un autre séparateur
        read.csv(fichier_csv, row.names = 1, check.names = FALSE, sep = ";")
    })
    
    if (nrow(df) != ncol(df)) {
        stop("❌ La matrice n'est pas carrée (", nrow(df), "x", ncol(df), ")")
    }
    
    d_matrix <- as.matrix(df)
    labels <- rownames(d_matrix)
    
    if (is.null(labels) || length(labels) == 0) {
        stop("❌ La matrice n'a pas de noms de lignes")
    }
    
    n <- nrow(d_matrix) - 2  # -2 pour dépôt et décharge
    
    message("[Chargement] ✅ Matrice chargée : ", nrow(d_matrix), " x ", ncol(d_matrix))
    message("[Chargement]   - Points candidats : ", n)
    
    # ── 3. Identifier dépôt et décharge dans les labels ──────
    idx_depot <- grep("Dépôt|Depot|depot|O", labels)
    if (length(idx_depot) == 0) {
        warning("[Chargement] ⚠️ 'Dépôt' non trouvé, utilisation de la première ligne")
        idx_depot <- 1
    } else {
        idx_depot <- idx_depot[1]
    }
    
    idx_decharge <- grep("Décharge|Decharge|decharge|S", labels)
    if (length(idx_decharge) == 0) {
        warning("[Chargement] ⚠️ 'Décharge' non trouvé, utilisation de la dernière ligne")
        idx_decharge <- nrow(d_matrix)
    } else {
        idx_decharge <- idx_decharge[1]
    }
    
    message("[Chargement]   - Dépôt : ligne ", idx_depot, " ('", labels[idx_depot], "')")
    message("[Chargement]   - Décharge : ligne ", idx_decharge, " ('", labels[idx_decharge], "')")
    
    # ── 4. Charger ou créer les coordonnées ────────────────────
    if (!is.null(fichier_coords) && file.exists(fichier_coords)) {
        coords <- read.csv(fichier_coords)
        message("[Chargement] 📍 Coordonnées chargées : ", nrow(coords), " points")
        
        if (nrow(coords) != n + 2) {
            warning("[Chargement] ⚠️ Incohérence : matrice ", n + 2, 
                    " nœuds, coordonnées ", nrow(coords))
        }
    } else {
        # Créer les coordonnées à partir des labels et des paramètres
        coords <- data.frame(
            index = seq_len(n + 2),
            role = c("depot", rep("candidat", n), "decharge"),
            id_local = c(0, seq_len(n), n + 1),
            longitude = NA_real_,
            latitude = NA_real_,
            label = labels,
            stringsAsFactors = FALSE
        )
        
        # ✅ CORRECTION : Utiliser les coordonnées de l'UI si fournies
        if (!is.null(depot_lon) && !is.null(depot_lat)) {
            coords$longitude[idx_depot] <- depot_lon
            coords$latitude[idx_depot] <- depot_lat
            message("[Chargement] 📍 Dépôt : coordonnées utilisateur (", depot_lon, ", ", depot_lat, ")")
        } else {
            # Essayer d'extraire des labels
            for (i in seq_len(nrow(coords))) {
                label <- coords$label[i]
                coord_match <- regexec("([0-9.-]+)[,;\\s]+([0-9.-]+)", label)
                match_result <- regmatches(label, coord_match)
                if (length(match_result) > 0 && length(match_result[[1]]) >= 3) {
                    coords$longitude[i] <- as.numeric(match_result[[1]][2])
                    coords$latitude[i] <- as.numeric(match_result[[1]][3])
                }
            }
            message("[Chargement] 📍 Coordonnées extraites des labels")
        }
        
        if (!is.null(decharge_lon) && !is.null(decharge_lat)) {
            coords$longitude[idx_decharge] <- decharge_lon
            coords$latitude[idx_decharge] <- decharge_lat
            message("[Chargement] 📍 Décharge : coordonnées utilisateur (", decharge_lon, ", ", decharge_lat, ")")
        }
        
        # Pour les points candidats sans coordonnées, utiliser des valeurs par défaut
        # mais avertir l'utilisateur
        if (any(is.na(coords$longitude) | is.na(coords$latitude))) {
            warning("[Chargement] ⚠️ Certaines coordonnées sont manquantes")
            # Remplacer les NA par 0 (fallback)
            coords$longitude[is.na(coords$longitude)] <- 0
            coords$latitude[is.na(coords$latitude)] <- 0
        }
    }
    
    # ── 5. Retour ──────────────────────────────────────────────
    list(
        d_matrix = d_matrix,
        coords_etendues = coords,
        n = n,
        idx_depot = as.integer(idx_depot),
        idx_decharge = as.integer(idx_decharge),
        source = "csv"
    )
}

#' Vérifier si un fichier est une matrice de distances valide
#' 
#' @param fichier_csv Chemin vers le fichier CSV
#' @param separateur Séparateur de colonnes (défaut: ",")
#' @return TRUE si valide, FALSE sinon
#' @export
verifier_matrice_csv <- function(fichier_csv, separateur = ",") {
    tryCatch({
        if (!file.exists(fichier_csv)) {
            message("[Vérification] ❌ Fichier non trouvé : ", fichier_csv)
            return(FALSE)
        }
        
        df <- read.csv(fichier_csv, row.names = 1, check.names = FALSE, sep = separateur)
        
        # Vérifier que c'est une matrice carrée
        if (nrow(df) < 3) {
            message("[Vérification] ❌ Matrice trop petite (", nrow(df), " lignes)")
            return(FALSE)
        }
        if (!all(rownames(df) == colnames(df))) {
            message("[Vérification] ❌ Labels lignes/colonnes différents")
            return(FALSE)
        }
        
        # Vérifier que les valeurs sont numériques
        if (!is.numeric(as.matrix(df))) {
            message("[Vérification] ❌ La matrice contient des valeurs non-numériques")
            return(FALSE)
        }
        
        message("[Vérification] ✅ Matrice valide : ", nrow(df), " x ", ncol(df))
        return(TRUE)
        
    }, error = function(e) {
        message("[Vérification] ❌ Erreur : ", e$message)
        return(FALSE)
    })
}

#' Exporter la matrice au format CSV
#' 
#' @param res_dist Liste de matrice (comme retournée par charger_matrice_csv)
#' @param fichier_sortie Chemin de sortie
#' @param separateur Séparateur de colonnes (défaut: ",")
#' @return Chemin du fichier créé
#' @export
exporter_matrice_csv <- function(res_dist, fichier_sortie, separateur = ",") {
    if (is.null(res_dist$d_matrix)) {
        stop("❌ res_dist$d_matrix est NULL")
    }
    
    df <- as.data.frame(res_dist$d_matrix)
    
    if (!is.null(res_dist$coords_etendues$label)) {
        rownames(df) <- res_dist$coords_etendues$label
        colnames(df) <- res_dist$coords_etendues$label
    }
    
    write.csv(df, file = fichier_sortie, row.names = TRUE)
    message("[Export] ✅ Matrice exportée : ", fichier_sortie)
    return(fichier_sortie)
}