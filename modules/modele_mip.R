# ============================================================
#  modules/modele_mip.R
#  Formulation et résolution du modèle LRP (Location Routing Problem)
#
#  Corrections par rapport au code original :
#    1. Dépôt (O) et décharge (S) intégrés comme nœuds réels
#    2. Contraintes de départ/arrivée (C1, C2) ajoutées
#    3. Conservation du flux sur les nœuds intermédiaires (C3)
#    4. Matrice de distances étendue (n+2 nœuds)
#    5. Distance par véhicule calculée en post-traitement
#
#  Convention d'index dans le modèle MIP :
#    - Les variables z[v, i, j] utilisent des index 1..(n+2)
#      où 1 = dépôt, 2..(n+1) = candidats, (n+2) = décharge
#    - Les variables y[j] et x[i,j] gardent l'index 1..n (candidats)
#
#  Fonction principale : resoudre_lrp()
# ============================================================

library(ompr)
library(ompr.roi)
library(ROI.plugin.glpk)


# ------------------------------------------------------------
#  resoudre_lrp()
#
#  Arguments :
#    res_dist  : liste retournée par construire_matrice_distances()
#    menages   : data.frame avec colonne poids_dechets
#    nb_vehicules : nombre de véhicules |V|
#    cap_camion   : capacité Qv (kg)
#    cap_point    : capacité Cj par point (kg) — valeur uniforme
#    d_max        : distance max ménage→point (mètres)
#    n_max        : budget nb points ouverts
#    solveur      : "glpk", "cbc" ou "gurobi"
#    time_limit   : limite en secondes
#
#  Retourne une liste :
#    $statut          : statut du solveur
#    $cout_total      : distance totale optimisée (km)
#    $points_ouverts  : vecteur des index candidats ouverts (1..n)
#    $sol_y           : data.frame complet des y[j]
#    $sol_z           : data.frame des arcs actifs (z=1) avec v, i, j
#    $sol_x           : data.frame des affectations ménages
#    $distances_par_vehicule : vecteur nommé (distance en km par véhicule)
#    $couverture      : % ménages affectés
#    $tournees_detail : data.frame lisible pour la table UI
# ------------------------------------------------------------

resoudre_lrp <- function(res_dist,
                         menages,
                         nb_vehicules,
                         cap_camion,
                         cap_point,
                         d_max,
                         n_max,
                         solveur    = "glpk",
                         time_limit = 120) {
    
    # ── Index de commodité ──────────────────────────────────────
    n      <- res_dist$n                  # nb points candidats
    m      <- nrow(menages)               # nb ménages
    V      <- nb_vehicules
    Qv     <- cap_camion
    Cj     <- rep(cap_point, n)           # capacité uniforme par point
    W      <- menages$poids_dechets       # poids de chaque ménage
    
    O      <- res_dist$idx_depot          # index dépôt    = 1
    S      <- res_dist$idx_decharge       # index décharge = n+2
    N_ext  <- n + 2                       # taille du graphe étendu
    d      <- res_dist$d_matrix           # matrice (N_ext x N_ext)
    
    # Index des candidats dans le graphe étendu : 2..(n+1)
    cand_ext <- seq(2, n + 1)
    
    # ── Masque d_max : xij=0 si distance ménage→candidat > d_max ──
    # (distances ménages→candidats en mètres via Haversine)
    coords_cand <- res_dist$coords_etendues[cand_ext, c("longitude","latitude")]
    coords_men  <- menages[, c("longitude","latitude")]
    
    d_men_cand <- geosphere::distm(
        as.matrix(coords_men),
        as.matrix(coords_cand),
        fun = geosphere::distHaversine
    )   # matrice m x n en mètres
    
    autorise_xij <- (d_men_cand <= d_max)   # TRUE si affectation permise
    
    # ── Modèle MIP ──────────────────────────────────────────────
    model <- MIPModel() %>%
        
        # ── Variables ──────────────────────────────────────────
        # y[j]      : point candidat j ouvert  (j dans 1..n)
        add_variable(y[j],    j = 1:n,                    type = "binary") %>%
        # x[i,j]   : ménage i affecté au candidat j
        add_variable(x[i, j], i = 1:m, j = 1:n,          type = "binary") %>%
        # z[v,i,j] : véhicule v emprunte arc i→j dans graphe étendu (1..N_ext)
        add_variable(z[v, i, j], v = 1:V, i = 1:N_ext, j = 1:N_ext, type = "binary") %>%
        # u[v,j]   : charge cumulée du véhicule v quand il quitte j (MTZ)
        add_variable(u[v, j], v = 1:V, j = 1:N_ext,      type = "continuous", lb = 0) %>%
        
        # ── Fonction objectif : distance totale des arcs empruntés ──
        set_objective(
            sum_expr(d[i, j] * z[v, i, j], v = 1:V, i = 1:N_ext, j = 1:N_ext),
            sense = "min"
        ) %>%
        
        # ── C1 : chaque véhicule part du dépôt exactement une fois ──
        add_constraint(
            sum_expr(z[v, O, j], j = cand_ext) == 1,
            v = 1:V
        ) %>%
        
        # ── C2 : chaque véhicule arrive à la décharge exactement une fois ──
        add_constraint(
            sum_expr(z[v, i, S], i = cand_ext) == 1,
            v = 1:V
        ) %>%
        
        # ── C3 : conservation du flux sur chaque nœud candidat ──
        # Ce qui entre = ce qui sort (pour chaque véhicule sur chaque candidat)
        add_constraint(
            sum_expr(z[v, i, j], i = 1:N_ext) == sum_expr(z[v, j, k], k = 1:N_ext),
            v = 1:V, j = cand_ext
        ) %>%
        
        # ── C4 : chaque point ouvert est visité par exactement un véhicule ──
        add_constraint(
            sum_expr(z[v, i, j], v = 1:V, i = 1:N_ext) == y[j - 1],
            j = cand_ext     # j-1 car y est indexé 1..n, j dans graphe étendu = 2..(n+1)
        ) %>%
        
        # ── C5 : pas d'arc d'un nœud vers lui-même ──
        add_constraint(z[v, i, i] == 0, v = 1:V, i = 1:N_ext) %>%
        
        # ── C6 : pas d'arc depuis la décharge ──
        add_constraint(z[v, S, j] == 0, v = 1:V, j = 1:N_ext) %>%
        
        # ── C7 : pas d'arc vers le dépôt ──
        add_constraint(z[v, i, O] == 0, v = 1:V, i = 1:N_ext) %>%
        
        # ── C8 : chaque ménage est affecté à exactement un point ──
        add_constraint(sum_expr(x[i, j], j = 1:n) == 1, i = 1:m) %>%
        
        # ── C9 : affectation seulement si le point est ouvert ──
        add_constraint(x[i, j] <= y[j], i = 1:m, j = 1:n) %>%
        
        # ── C10 : respect de la distance max ménage→point ──
        add_constraint(x[i, j] <= as.integer(autorise_xij[i, j]), i = 1:m, j = 1:n) %>%
        
        # ── C11 : capacité du point de collecte ──
        add_constraint(
            sum_expr(W[i] * x[i, j], i = 1:m) <= Cj[j] * y[j],
            j = 1:n
        ) %>%
        
        # ── C12 : nombre max de points ouverts (contrainte budget) ──
        add_constraint(sum_expr(y[j], j = 1:n) <= n_max) %>%
        
        # ── C13 : MTZ — élimination des sous-tours ──
        # Charge cumulée : si le véhicule v va de i vers j, u augmente de Cj
        add_constraint(
            u[v, j] >= u[v, i] + Cj[j - 1] * z[v, i, j] - Qv * (1 - z[v, i, j]),
            v = 1:V,
            i = cand_ext,
            j = cand_ext
        ) %>%
        
        # ── C14 : borne supérieure de la charge ──
        add_constraint(u[v, j] <= Qv, v = 1:V, j = 1:N_ext) %>%
        
        # ── C15 : charge initiale au dépôt = 0 ──
        add_constraint(u[v, O] == 0, v = 1:V)
    
    # ── Résolution ──────────────────────────────────────────────
    result <- solve_model(
        model,
        with_ROI(
            solver  = solveur,
            control = list(tm_limit = time_limit * 1000,
                           verbose  = TRUE)
        )
    )
    
    # ── Extraction des solutions ─────────────────────────────────
    statut   <- solver_status(result)
    cout_obj <- objective_value(result)
    
    sol_y <- get_solution(result, y[j])
    sol_x <- get_solution(result, x[i, j])
    sol_z <- get_solution(result, z[v, i, j])
    
    points_ouverts <- sol_y$j[sol_y$value > 0.5]
    
    arcs_actifs <- sol_z[sol_z$value > 0.5, ]   # uniquement z=1
    
    # ── Distance par véhicule ────────────────────────────────────
    distances_par_vehicule <- sapply(1:V, function(v_id) {
        arcs_v <- arcs_actifs[arcs_actifs$v == v_id, ]
        if (nrow(arcs_v) == 0) return(0)
        sum(sapply(seq_len(nrow(arcs_v)), function(k) {
            d[arcs_v$i[k], arcs_v$j[k]]
        }))
    })
    names(distances_par_vehicule) <- paste0("Vehicule_", 1:V)
    
    # ── Couverture ───────────────────────────────────────────────
    menages_couverts <- sum(rowSums(sol_x[sol_x$value > 0.5, c("value")]) > 0,
                            na.rm = TRUE)
    # Plus simple :
    affectes <- unique(sol_x$i[sol_x$value > 0.5])
    couverture_pct <- round(100 * length(affectes) / m, 1)
    
    # ── Table lisible pour l'UI ──────────────────────────────────
    tournees_detail <- arcs_actifs %>%
        dplyr::mutate(
            distance_km = mapply(function(ii, jj) round(d[ii, jj], 3), i, j),
            noeud_depart = dplyr::case_when(
                i == O ~ "Dépôt",
                i == S ~ "Décharge",
                TRUE   ~ paste0("Point_", i - 1)   # i-1 car dépôt occupe index 1
            ),
            noeud_arrivee = dplyr::case_when(
                j == O ~ "Dépôt",
                j == S ~ "Décharge",
                TRUE   ~ paste0("Point_", j - 1)
            )
        ) %>%
        dplyr::select(Vehicule = v, De = noeud_depart, Vers = noeud_arrivee, Distance_km = distance_km) %>%
        dplyr::arrange(Vehicule, De)
    
    # ── Retour ───────────────────────────────────────────────────
    list(
        statut                  = statut,
        cout_total              = round(cout_obj, 3),
        points_ouverts          = points_ouverts,
        sol_y                   = sol_y,
        sol_z                   = arcs_actifs,
        sol_x                   = sol_x[sol_x$value > 0.5, ],
        distances_par_vehicule  = round(distances_par_vehicule, 3),
        couverture              = couverture_pct,
        tournees_detail         = tournees_detail,
        nb_vehicules            = V,
        nb_points               = length(points_ouverts)
    )
}


# ------------------------------------------------------------
#  resumer_solution()
#  Affiche un résumé console de la solution (debug / log solveur)
# ------------------------------------------------------------

resumer_solution <- function(sol) {
    cat("=== Solution LRP ===\n")
    cat("Statut           :", sol$statut, "\n")
    cat("Distance totale  :", sol$cout_total, "km\n")
    cat("Points ouverts   :", sol$nb_points, "\n")
    cat("Véhicules        :", sol$nb_vehicules, "\n")
    cat("Couverture       :", sol$couverture, "%\n\n")
    cat("Distance par véhicule :\n")
    for (nm in names(sol$distances_par_vehicule)) {
        cat(" ", nm, ":", sol$distances_par_vehicule[nm], "km\n")
    }
    cat("\nTournées :\n")
    print(sol$tournees_detail)
}