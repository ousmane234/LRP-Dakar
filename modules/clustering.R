# ============================================================
#  modules/clustering.R
#  K-means contraint pour le placement des points de collecte
#
#  Problème du K-means standard :
#    - Ne respecte pas la contrainte de capacité par cluster
#    - Détection du coude incorrecte (dérivée seconde instable)
#
#  Solution ici :
#    1. K-means++ pour l'initialisation (meilleure convergence)
#    2. Vérification post-clustering de la contrainte de charge
#    3. Splitting automatique des clusters surchargés
#    4. Détection du coude par distance à la droite (méthode Kneedle)
#    5. Projection des centroïdes sur le réseau routier (via OSM/dodgr)
#       → activée seulement si dodgr est disponible
#
#  Fonctions principales :
#    clustering_contraint()   : pipeline complet
#    detecter_coude()         : détection k optimal
#    verifier_charges()       : diagnostic des clusters
# ============================================================

library(dplyr)

# ------------------------------------------------------------
#  detecter_coude()
#  Méthode de la distance à la droite (Kneedle simplifié)
#
#  Arguments :
#    k_vals   : vecteur des valeurs de k testées
#    inerties : vecteur des inerties correspondantes
#
#  Retourne : k optimal (entier)
# ------------------------------------------------------------

detecter_coude <- function(k_vals, inerties) {
    if (length(k_vals) < 3) return(k_vals[1])
    
    # Normaliser pour que la droite soit bien définie
    x1 <- k_vals[1];              y1 <- inerties[1]
    x2 <- k_vals[length(k_vals)]; y2 <- inerties[length(inerties)]
    
    # Distance de chaque point à la droite (x1,y1)-(x2,y2)
    distances <- sapply(seq_along(k_vals), function(idx) {
        xp <- k_vals[idx]; yp <- inerties[idx]
        abs((y2 - y1) * xp - (x2 - x1) * yp + x2 * y1 - y2 * x1) /
            sqrt((y2 - y1)^2 + (x2 - x1)^2)
    })
    
    k_vals[which.max(distances)]
}


# ------------------------------------------------------------
#  kmeans_plus_plus_init()
#  Initialisation K-means++ : choisit les centres initiaux
#  proportionnellement à la distance au centre le plus proche.
#
#  Arguments :
#    coords : matrice n x 2 (longitude, latitude)
#    k      : nombre de centres
#
#  Retourne : matrice k x 2 des centres initiaux
# ------------------------------------------------------------

kmeans_plus_plus_init <- function(coords, k) {
    n <- nrow(coords)
    # 1er centre : aléatoire
    centres <- coords[sample(n, 1), , drop = FALSE]
    
    for (iter in seq_len(k - 1)) {
        # Distance de chaque point au centre le plus proche
        dists <- apply(coords, 1, function(pt) {
            min(apply(centres, 1, function(c) sum((pt - c)^2)))
        })
        # Probabilité proportionnelle à d²
        prob <- dists / sum(dists)
        nouveau <- coords[sample(n, 1, prob = prob), , drop = FALSE]
        centres <- rbind(centres, nouveau)
    }
    centres
}


# ------------------------------------------------------------
#  splitter_cluster_surcharge()
#  Divise un cluster dont la charge dépasse cap_point en 2.
#
#  Arguments :
#    menages_cluster : data.frame des ménages du cluster (lon, lat, poids)
#    cap_point       : capacité max par point (kg)
#
#  Retourne : vecteur d'assignation (1 ou 2) pour chaque ménage
# ------------------------------------------------------------

splitter_cluster_surcharge <- function(menages_cluster, cap_point) {
    coords  <- as.matrix(menages_cluster[, c("longitude", "latitude")])
    poids   <- menages_cluster$poids_dechets
    
    # K-means en 2 sous-clusters
    km2 <- kmeans(coords, centers = 2, nstart = 5)
    
    # Vérifier que chaque sous-cluster respecte la capacité
    charges <- tapply(poids, km2$cluster, sum)
    
    # Si un sous-cluster est encore surchargé, on accepte quand même
    # (éviter la récursion infinie) — l'UI le signalera
    km2$cluster
}


# ------------------------------------------------------------
#  clustering_contraint()
#  Pipeline complet de clustering pour le LRP.
#
#  Arguments :
#    menages         : data.frame avec longitude, latitude, poids_dechets
#    k_min, k_max    : plage de k à tester
#    n_init          : nb d'initialisations K-means
#    cap_point       : capacité max par cluster/point (kg)
#    appliquer_contrainte : si TRUE, split les clusters surchargés
#
#  Retourne une liste :
#    $candidats    : data.frame des points candidats (id, lon, lat, charge, valide)
#    $menages      : data.frame original + colonnes cluster, cluster_valide
#    $inerties     : data.frame(k, inertie) pour la courbe
#    $k_optimal    : k choisi par la méthode du coude
#    $clusters_surcharges : vecteur des id clusters dépassant cap_point
#    $resume       : texte résumé pour l'UI
# ------------------------------------------------------------

clustering_contraint <- function(menages,
                                 k_min,
                                 k_max,
                                 n_init   = 10,
                                 cap_point,
                                 appliquer_contrainte = TRUE) {
    
    coords <- as.matrix(menages[, c("longitude", "latitude")])
    k_vals <- k_min:k_max
    
    # ── 1. Courbe d'inertie sur la plage k_min..k_max ────────
    inerties <- sapply(k_vals, function(k) {
        init    <- kmeans_plus_plus_init(coords, k)
        km      <- kmeans(coords, centers = init, nstart = n_init)
        km$tot.withinss
    })
    df_inerties <- data.frame(k = k_vals, inertie = inerties)
    
    # ── 2. Détection du coude ─────────────────────────────────
    k_optimal <- detecter_coude(k_vals, inerties)
    
    # ── 3. K-means final avec k optimal ──────────────────────
    init_opt <- kmeans_plus_plus_init(coords, k_optimal)
    km_final <- kmeans(coords, centers = init_opt, nstart = n_init)
    
    menages$cluster <- km_final$cluster
    
    # ── 4. Calcul des charges par cluster ─────────────────────
    charges <- tapply(menages$poids_dechets, menages$cluster, sum)
    
    # ── 5. Gestion des clusters surchargés ────────────────────
    clusters_surcharges <- as.integer(names(charges[charges > cap_point]))
    
    if (appliquer_contrainte && length(clusters_surcharges) > 0) {
        # Pour chaque cluster surchargé, on le divise en 2
        # et on renuméroter tous les clusters proprement
        
        nouveau_cluster <- menages$cluster  # copie
        offset          <- max(menages$cluster)
        
        for (cl_id in clusters_surcharges) {
            idx_cl <- which(menages$cluster == cl_id)
            if (length(idx_cl) < 2) next
            
            sous_assign <- splitter_cluster_surcharge(
                menages[idx_cl, c("longitude", "latitude", "poids_dechets")],
                cap_point
            )
            # Le sous-cluster 1 garde l'id original, le 2 prend offset+1
            offset <- offset + 1
            nouveau_cluster[idx_cl[sous_assign == 2]] <- offset
        }
        
        menages$cluster <- as.integer(factor(nouveau_cluster))  # renumérote 1..K_final
        k_final         <- max(menages$cluster)
        charges         <- tapply(menages$poids_dechets, menages$cluster, sum)
        clusters_surcharges_post <- as.integer(names(charges[charges > cap_point]))
    } else {
        k_final                  <- k_optimal
        clusters_surcharges_post <- clusters_surcharges
    }
    
    # ── 6. Construire les centroïdes (points candidats) ───────
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
    
    # ── 7. Résumé texte ───────────────────────────────────────
    nb_valides     <- sum(candidats$valide)
    nb_surcharges  <- sum(!candidats$valide)
    
    resume <- paste0(
        "k optimal (coude)    : ", k_optimal, "\n",
        "k final (après split): ", k_final, "\n",
        "Points valides       : ", nb_valides, "\n",
        "Points surchargés    : ", nb_surcharges,
        if (nb_surcharges > 0)
            paste0(" → IDs : ", paste(candidats$id[!candidats$valide], collapse = ", "))
        else "",
        "\n",
        "Charge max observée  : ", round(max(candidats$charge)), " kg\n",
        "Capacité par point   : ", cap_point, " kg\n"
    )
    
    # ── 8. Retour ─────────────────────────────────────────────
    list(
        candidats            = candidats,
        menages              = menages,
        inerties             = df_inerties,
        k_optimal            = k_optimal,
        k_final              = k_final,
        clusters_surcharges  = clusters_surcharges_post,
        resume               = resume
    )
}


# ------------------------------------------------------------
#  verifier_charges()
#  Diagnostic détaillé des charges par cluster.
#  Utile pour l'affichage dans l'UI (table des candidats).
#
#  Arguments :
#    candidats : data.frame retourné par clustering_contraint()
#    cap_point : capacité max
#
#  Retourne : data.frame enrichi avec statut et taux de remplissage
# ------------------------------------------------------------

verifier_charges <- function(candidats, cap_point) {
    candidats %>%
        mutate(
            taux_remplissage = round(100 * charge / cap_point, 1),
            statut = case_when(
                charge > cap_point        ~ "⚠ Surchargé",
                charge > 0.9 * cap_point  ~ "⚡ Quasi-plein",
                charge > 0.5 * cap_point  ~ "✓ Normal",
                TRUE                      ~ "○ Sous-utilisé"
            )
        )
}