# ============================================================
#  modules/clustering.R
#  K-means contraint pour le placement des points de collecte
# ============================================================

library(dplyr)

detecter_coude <- function(k_vals, inerties) {
    if (length(k_vals) < 3) return(k_vals[1])
    
    x1 <- k_vals[1];              y1 <- inerties[1]
    x2 <- k_vals[length(k_vals)]; y2 <- inerties[length(inerties)]
    
    distances <- sapply(seq_along(k_vals), function(idx) {
        xp <- k_vals[idx]; yp <- inerties[idx]
        abs((y2 - y1) * xp - (x2 - x1) * yp + x2 * y1 - y2 * x1) /
            sqrt((y2 - y1)^2 + (x2 - x1)^2)
    })
    
    k_vals[which.max(distances)]
}

kmeans_plus_plus_init <- function(coords, k) {
    n <- nrow(coords)
    centres <- coords[sample(n, 1), , drop = FALSE]
    
    for (iter in seq_len(k - 1)) {
        dists <- apply(coords, 1, function(pt) {
            min(apply(centres, 1, function(c) sum((pt - c)^2)))
        })
        prob <- dists / sum(dists)
        nouveau <- coords[sample(n, 1, prob = prob), , drop = FALSE]
        centres <- rbind(centres, nouveau)
    }
    centres
}

splitter_cluster_surcharge <- function(menages_cluster, cap_point) {
    coords <- as.matrix(menages_cluster[, c("longitude", "latitude")])
    km2 <- kmeans(coords, centers = 2, nstart = 5)
    km2$cluster
}

clustering_contraint <- function(menages,
                                 k_min,
                                 k_max,
                                 n_init   = 10,
                                 cap_point,
                                 appliquer_contrainte = TRUE) {
    
    coords <- as.matrix(menages[, c("longitude", "latitude")])
    k_vals <- k_min:k_max
    
    inerties <- sapply(k_vals, function(k) {
        init <- kmeans_plus_plus_init(coords, k)
        km <- kmeans(coords, centers = init, nstart = n_init)
        km$tot.withinss
    })
    df_inerties <- data.frame(k = k_vals, inertie = inerties)
    
    k_optimal <- detecter_coude(k_vals, inerties)
    
    init_opt <- kmeans_plus_plus_init(coords, k_optimal)
    km_final <- kmeans(coords, centers = init_opt, nstart = n_init)
    
    menages$cluster <- km_final$cluster
    charges <- tapply(menages$poids_dechets, menages$cluster, sum)
    clusters_surcharges <- as.integer(names(charges[charges > cap_point]))
    
    if (appliquer_contrainte && length(clusters_surcharges) > 0) {
        nouveau_cluster <- menages$cluster
        offset <- max(menages$cluster)
        
        for (cl_id in clusters_surcharges) {
            idx_cl <- which(menages$cluster == cl_id)
            if (length(idx_cl) < 2) next
            
            sous_assign <- splitter_cluster_surcharge(
                menages[idx_cl, c("longitude", "latitude", "poids_dechets")],
                cap_point
            )
            offset <- offset + 1
            nouveau_cluster[idx_cl[sous_assign == 2]] <- offset
        }
        
        menages$cluster <- as.integer(factor(nouveau_cluster))
        k_final <- max(menages$cluster)
        charges <- tapply(menages$poids_dechets, menages$cluster, sum)
        clusters_surcharges_post <- as.integer(names(charges[charges > cap_point]))
    } else {
        k_final <- k_optimal
        clusters_surcharges_post <- clusters_surcharges
    }
    
    centres <- do.call(rbind, lapply(sort(unique(menages$cluster)), function(cl) {
        idx <- which(menages$cluster == cl)
        data.frame(
            longitude = mean(menages$longitude[idx]),
            latitude  = mean(menages$latitude[idx])
        )
    }))
    
    candidats <- data.frame(
        id        = seq_len(nrow(centres)),
        longitude = centres$longitude,
        latitude  = centres$latitude,
        cluster   = seq_len(nrow(centres)),
        charge    = as.numeric(charges[order(as.integer(names(charges)))]),
        valide    = as.numeric(charges[order(as.integer(names(charges)))]) <= cap_point
    )
    
    list(
        candidats            = candidats,
        menages              = menages,
        inerties             = df_inerties,
        k_optimal            = k_optimal,
        k_final              = k_final,
        clusters_surcharges  = clusters_surcharges_post
    )
}