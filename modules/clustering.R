# ============================================================
#  modules/clustering.R
#  K-means contraint pour le placement des points de collecte
# ============================================================

library(dplyr)

#' Détection du coude par distance à la ligne
#' 
#' @param k_vals Vecteur des valeurs de k
#' @param inerties Vecteur des inerties correspondantes
#' @return La valeur de k optimale
#' @export
detecter_coude <- function(k_vals, inerties) {
    # Vérifications de base
    if (length(k_vals) < 2) {
        return(k_vals[1])
    }
    if (length(k_vals) < 3) {
        # Avec seulement 2 points, prendre celui avec la plus forte pente
        pente <- (inerties[2] - inerties[1]) / (k_vals[2] - k_vals[1])
        return(ifelse(pente > 0, k_vals[1], k_vals[2]))
    }
    
    # Méthode de la distance à la ligne entre le premier et le dernier point
    x1 <- k_vals[1];              y1 <- inerties[1]
    x2 <- k_vals[length(k_vals)]; y2 <- inerties[length(inerties)]
    
    # Éviter la division par zéro
    denom <- sqrt((y2 - y1)^2 + (x2 - x1)^2)
    if (denom == 0) {
        return(k_vals[which.min(inerties)])
    }
    
    distances <- sapply(seq_along(k_vals), function(idx) {
        xp <- k_vals[idx]; yp <- inerties[idx]
        abs((y2 - y1) * xp - (x2 - x1) * yp + x2 * y1 - y2 * x1) / denom
    })
    
    k_vals[which.max(distances)]
}

#' Initialisation K-means++ améliorée
#' 
#' @param coords Matrice des coordonnées
#' @param k Nombre de clusters
#' @return Matrice des centres initiaux
#' @export
kmeans_plus_plus_init <- function(coords, k) {
    n <- nrow(coords)
    
    # Vérifier qu'il y a assez de points
    if (n < k) {
        warning("[Clustering] k (", k, ") > nombre de points (", n, "), réduction à ", n)
        k <- max(1, n)
    }
    
    if (k == 1) {
        return(colMeans(coords))
    }
    
    # Premier centre : un point aléatoire
    centres <- coords[sample(n, 1), , drop = FALSE]
    
    for (iter in seq_len(k - 1)) {
        # Distance au centre le plus proche pour chaque point
        dists <- apply(coords, 1, function(pt) {
            min(apply(centres, 1, function(c) sum((pt - c)^2)))
        })
        
        # Éviter les probabilités nulles
        if (sum(dists) == 0) {
            # Si toutes les distances sont nulles, prendre un point aléatoire
            nouveau <- coords[sample(n, 1), , drop = FALSE]
        } else {
            prob <- dists / sum(dists)
            # S'assurer qu'aucune probabilité n'est NA
            prob[is.na(prob)] <- 0
            if (sum(prob) == 0) prob <- rep(1/n, n)
            nouveau <- coords[sample(n, 1, prob = prob), , drop = FALSE]
        }
        
        centres <- rbind(centres, nouveau)
    }
    
    centres
}

#' Splitter un cluster surchargé en 2 sous-clusters
#' 
#' @param menages_cluster Data.frame des ménages d'un cluster
#' @param cap_point Capacité maximale du point
#' @return Vecteur des clusters (1 ou 2)
#' @export
splitter_cluster_surcharge <- function(menages_cluster, cap_point) {
    # Vérifier qu'il y a assez de points
    if (nrow(menages_cluster) < 2) {
        return(rep(1, nrow(menages_cluster)))
    }
    
    coords <- as.matrix(menages_cluster[, c("longitude", "latitude")])
    
    # Si tous les points sont identiques, les garder ensemble
    if (var(coords[,1]) == 0 && var(coords[,2]) == 0) {
        return(rep(1, nrow(menages_cluster)))
    }
    
    # K-means avec 2 clusters
    km2 <- tryCatch({
        kmeans(coords, centers = 2, nstart = 5)
    }, error = function(e) {
        # Fallback : diviser aléatoirement
        message("[Clustering] ⚠️ Erreur kmeans split, division aléatoire")
        return(list(cluster = sample(1:2, nrow(menages_cluster), replace = TRUE)))
    })
    
    km2$cluster
}

#' Clustering contraint avec K-means
#' 
#' @param menages Data.frame des ménages (doit contenir longitude, latitude, poids_dechets)
#' @param k_min Nombre minimum de clusters
#' @param k_max Nombre maximum de clusters
#' @param n_init Nombre d'initialisations K-means
#' @param cap_point Capacité maximale d'un point de collecte
#' @param appliquer_contrainte Appliquer la contrainte de charge ?
#' @return Liste contenant candidats, menages, inerties, etc.
#' @export
clustering_contraint <- function(menages,
                                 k_min,
                                 k_max,
                                 n_init   = 10,
                                 cap_point,
                                 appliquer_contrainte = TRUE) {
    
    # ── Vérifications ──────────────────────────────────────────
    if (is.null(menages) || nrow(menages) == 0) {
        stop("❌ Aucun ménage à clusteriser")
    }
    
    cols_requises <- c("longitude", "latitude", "poids_dechets")
    manquantes <- setdiff(cols_requises, names(menages))
    if (length(manquantes) > 0) {
        stop("❌ Colonnes manquantes : ", paste(manquantes, collapse = ", "))
    }
    
    # Ajuster k_min et k_max
    n_points <- nrow(menages)
    if (k_min < 1) k_min <- 1
    if (k_max > n_points) k_max <- n_points
    if (k_min > k_max) k_min <- k_max
    
    coords <- as.matrix(menages[, c("longitude", "latitude")])
    k_vals <- k_min:k_max
    
    if (length(k_vals) == 0) {
        stop("❌ Intervalle k vide : k_min=", k_min, ", k_max=", k_max)
    }
    
    message("[Clustering] Début : ", n_points, " ménages, k de ", k_min, " à ", k_max)
    
    # ── Calcul des inerties ────────────────────────────────────
    inerties <- sapply(k_vals, function(k) {
        init <- kmeans_plus_plus_init(coords, k)
        km <- kmeans(coords, centers = init, nstart = n_init)
        km$tot.withinss
    })
    
    df_inerties <- data.frame(k = k_vals, inertie = inerties)
    
    # ── Sélection du k optimal ────────────────────────────────
    k_optimal <- detecter_coude(k_vals, inerties)
    message("[Clustering] k optimal : ", k_optimal)
    
    # ── Clustering final ──────────────────────────────────────
    init_opt <- kmeans_plus_plus_init(coords, k_optimal)
    km_final <- kmeans(coords, centers = init_opt, nstart = n_init)
    
    menages$cluster <- km_final$cluster
    
    # ── Gestion des surcharges ────────────────────────────────
    charges <- tapply(menages$poids_dechets, menages$cluster, sum)
    clusters_surcharges <- as.integer(names(charges[charges > cap_point]))
    
    if (appliquer_contrainte && length(clusters_surcharges) > 0) {
        message("[Clustering] ⚠️ ", length(clusters_surcharges), " cluster(s) surchargé(s)")
        
        nouveau_cluster <- menages$cluster
        offset <- max(menages$cluster)
        
        for (cl_id in clusters_surcharges) {
            idx_cl <- which(menages$cluster == cl_id)
            if (length(idx_cl) < 2) {
                message("[Clustering] ⚠️ Cluster ", cl_id, " trop petit pour être splité")
                next
            }
            
            sous_assign <- splitter_cluster_surcharge(
                menages[idx_cl, c("longitude", "latitude", "poids_dechets")],
                cap_point
            )
            
            if (length(unique(sous_assign)) > 1) {
                offset <- offset + 1
                nouveau_cluster[idx_cl[sous_assign == 2]] <- offset
                message("[Clustering] ✅ Cluster ", cl_id, " splité en 2")
            }
        }
        
        # Mettre à jour les clusters
        menages$cluster <- as.integer(factor(nouveau_cluster))
        k_final <- max(menages$cluster)
        
        # Recalculer les charges
        charges <- tapply(menages$poids_dechets, menages$cluster, sum)
        clusters_surcharges_post <- as.integer(names(charges[charges > cap_point]))
        
        if (length(clusters_surcharges_post) > 0) {
            message("[Clustering] ⚠️ ", length(clusters_surcharges_post), 
                    " cluster(s) restent surchargé(s)")
        }
        
    } else {
        k_final <- k_optimal
        clusters_surcharges_post <- clusters_surcharges
    }
    
    # ── Calcul des centres ─────────────────────────────────────
    clusters_uniq <- sort(unique(menages$cluster))
    centres <- do.call(rbind, lapply(clusters_uniq, function(cl) {
        idx <- which(menages$cluster == cl)
        data.frame(
            longitude = mean(menages$longitude[idx], na.rm = TRUE),
            latitude  = mean(menages$latitude[idx], na.rm = TRUE)
        )
    }))
    
    # Créer le data.frame des candidats dans le bon ordre
    charges_sorted <- charges[as.character(clusters_uniq)]
    candidats <- data.frame(
        id        = seq_along(clusters_uniq),
        longitude = centres$longitude,
        latitude  = centres$latitude,
        cluster   = clusters_uniq,
        charge    = as.numeric(charges_sorted),
        valide    = as.numeric(charges_sorted) <= cap_point,
        stringsAsFactors = FALSE
    )
    
    message("[Clustering] ✅ Terminé : ", nrow(candidats), " points, k final = ", k_final)
    
    list(
        candidats            = candidats,
        menages              = menages,
        inerties             = df_inerties,
        k_optimal            = k_optimal,
        k_final              = k_final,
        clusters_surcharges  = clusters_surcharges_post
    )
}