# ============================================================
#  modules/chargement_matrice.R
#  Chargement de la matrice de distances depuis un fichier CSV
#  (généré par scripts/calcul_matrice_dakar.R)
# ============================================================

#' Charger la matrice de distances depuis un fichier CSV
#' 
#' @param fichier_csv Chemin vers le fichier CSV de la matrice
#' @param fichier_coords Chemin vers le fichier CSV des coordonnées (optionnel)
#' @return Liste avec la même structure que res_dist
charger_matrice_csv <- function(fichier_csv, fichier_coords = NULL) {
    
    # ── 1. Charger la matrice ──────────────────────────────────
    df <- read.csv(fichier_csv, row.names = 1, check.names = FALSE)
    d_matrix <- as.matrix(df)
    
    # Extraire les labels
    labels <- rownames(d_matrix)
    n <- length(labels) - 2  # -2 pour dépôt et décharge
    
    message("[Chargement] Matrice chargée : ", nrow(d_matrix), " x ", ncol(d_matrix))
    
    # ── 2. Charger les coordonnées ─────────────────────────────
    if (!is.null(fichier_coords) && file.exists(fichier_coords)) {
        coords <- read.csv(fichier_coords)
        message("[Chargement] Coordonnées chargées : ", nrow(coords), " points")
    } else {
        # Créer des coordonnées à partir des labels
        coords <- data.frame(
            index = 1:(n + 2),
            role = c("depot", rep("candidat", n), "decharge"),
            id_local = c(0, seq_len(n), n + 1),
            longitude = rep(0, n + 2),
            latitude = rep(0, n + 2),
            label = labels,
            stringsAsFactors = FALSE
        )
        message("[Chargement] ⚠️ Coordonnées non fournies, création factice")
    }
    
    # ── 3. Vérifier la cohérence ───────────────────────────────
    if (nrow(coords) != n + 2) {
        warning("[Chargement] Incohérence : matrice ", n + 2, 
                " nœuds, coordonnées ", nrow(coords))
    }
    
    # ── 4. Retour ──────────────────────────────────────────────
    list(
        d_matrix = d_matrix,
        coords_etendues = coords,
        n = n,
        idx_depot = 1L,
        idx_decharge = as.integer(n + 2),
        source = "csv"
    )
}

#' Vérifier si un fichier est une matrice de distances valide
#' 
#' @param fichier_csv Chemin vers le fichier CSV
#' @return TRUE si valide, FALSE sinon
verifier_matrice_csv <- function(fichier_csv) {
    tryCatch({
        df <- read.csv(fichier_csv, row.names = 1, check.names = FALSE)
        
        # Vérifier que c'est une matrice carrée
        if (nrow(df) < 3) {
            message("[Vérification] ❌ Matrice trop petite")
            return(FALSE)
        }
        if (!all(rownames(df) == colnames(df))) {
            message("[Vérification] ❌ Labels lignes/colonnes différents")
            return(FALSE)
        }
        
        # Vérifier qu'il y a "Dépôt" et "Décharge"
        labels <- rownames(df)
        if (!any(grepl("Dépôt", labels))) {
            message("[Vérification] ❌ 'Dépôt' manquant")
            return(FALSE)
        }
        if (!any(grepl("Décharge", labels))) {
            message("[Vérification] ❌ 'Décharge' manquant")
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
exporter_matrice_csv <- function(res_dist, fichier_sortie) {
    df <- as.data.frame(res_dist$d_matrix)
    rownames(df) <- res_dist$coords_etendues$label
    colnames(df) <- res_dist$coords_etendues$label
    write.csv(df, file = fichier_sortie)
    message("[Export] ✅ Matrice exportée : ", fichier_sortie)
}