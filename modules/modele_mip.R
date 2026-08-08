# ============================================================
#  modules/modele_mip.R - VERSION AVEC OSM
#  Modele LRP avec distances OSM
#  ------------------------------------------------------------
#  Optimisation lexicographique en deux etapes :
#    Etape A : maximiser la couverture reelle (menages a <= Dmax
#              de leur point affecte), sous budget Nmax et
#              capacite. Modele leger (pas de routage).
#    Etape B : minimiser la distance des tournees, en garantissant
#              une couverture au moins egale a (couverture_max -
#              tol) trouvee en etape A. C'est le modele complet
#              (localisation + affectation + routage).
#
#  Ce choix remplace une ancienne version a somme ponderee
#  (poids_desserte) : dans un MILP (variables entieres), une somme
#  ponderee (a) demande un poids arbitraire a calibrer, (b) peut
#  ne jamais atteindre certains optima de compromis a cause de la
#  non-convexite du probleme entier. L'approche lexicographique /
#  epsilon-contrainte evite les deux ecueils : aucun poids a
#  deviner, et la couverture maximale est garantie par construction,
#  jamais sacrifiee pour quelques km de tournee en moins.
# ============================================================

library(ompr)
library(ompr.roi)
library(ROI.plugin.glpk)


# ------------------------------------------------------------
#  calculer_distances_menages_candidats()
#  Calcule la matrice des distances menage -> candidat (metres),
#  via le graphe OSM si disponible, sinon Haversine en repli.
# ------------------------------------------------------------
calculer_distances_menages_candidats <- function(res_dist, menages, cand_ext) {
    
    if (!is.null(res_dist$graphe) && !is.null(res_dist$coords_snap)) {
        message("[MIP] Utilisation du graphe OSM pour les distances menage-candidat")
        
        coords_men  <- as.matrix(menages[, c("longitude", "latitude")])
        coords_cand <- as.matrix(res_dist$coords_snap[cand_ext, c("snap_lon", "snap_lat")])
        
        d_men_cand <- tryCatch({
            dodgr::dodgr_dists(res_dist$graphe, from = coords_men, to = coords_cand)
        }, error = function(e) {
            message("[MIP] Erreur OSM, fallback Haversine : ", e$message)
            geosphere::distm(coords_men, coords_cand, fun = geosphere::distHaversine)
        })
        
        if (any(is.na(d_men_cand))) {
            message("[MIP] Des distances OSM sont NA, fallback Haversine")
            d_men_cand_hav <- geosphere::distm(coords_men, coords_cand,
                                               fun = geosphere::distHaversine)
            d_men_cand[is.na(d_men_cand)] <- d_men_cand_hav[is.na(d_men_cand)]
        }
        
    } else {
        message("[MIP] Fallback Haversine pour les distances menage-candidat")
        coords_cand <- as.matrix(res_dist$coords_etendues[cand_ext, c("longitude", "latitude")])
        coords_men  <- as.matrix(menages[, c("longitude", "latitude")])
        d_men_cand  <- geosphere::distm(coords_men, coords_cand, fun = geosphere::distHaversine)
    }
    
    d_men_cand
}


# ------------------------------------------------------------
#  resoudre_couverture_max()  --  ETAPE A
#  Modele leger (x, y, c uniquement -- pas de routage) qui calcule
#  la couverture maximale atteignable sous Nmax et la capacite.
#  Sert de plafond de reference pour l'etape B, et peut aussi etre
#  appele seul pour un diagnostic rapide (ex. balayage de Nmax).
#
#  Retourne une liste : couverture_max_n, couverture_max_pct,
#  y_ref, x_ref, statut.
# ------------------------------------------------------------
resoudre_couverture_max <- function(d_men_cand, W, Cj, n_max, d_max,
                                    solveur = "glpk", time_limit = 60) {
    
    m <- nrow(d_men_cand)
    n <- ncol(d_men_cand)
    
    model_a <- MIPModel() %>%
        add_variable(y[j],    j = 1:n,           type = "binary") %>%
        add_variable(x[i, j], i = 1:m, j = 1:n,  type = "binary") %>%
        add_variable(c[i],    i = 1:m,           type = "binary") %>%
        
        set_objective(sum_expr(c[i], i = 1:m), sense = "max") %>%
        
        add_constraint(sum_expr(x[i, j], j = 1:n) == 1, i = 1:m) %>%
        add_constraint(x[i, j] <= y[j], i = 1:m, j = 1:n) %>%
        add_constraint(sum_expr(y[j], j = 1:n) <= n_max) %>%
        add_constraint(
            sum_expr(W[i] * x[i, j], i = 1:m) <= Cj[j] * y[j],
            j = 1:n
        ) %>%
        add_constraint(
            c[i] <= sum_expr(x[i, j], j = 1:n, d_men_cand[i, j] <= d_max),
            i = 1:m
        )
    
    result_a <- tryCatch({
        solve_model(model_a, with_ROI(solver = solveur,
                                      control = list(tm_limit = time_limit * 1000, verbose = FALSE)))
    }, error = function(e) {
        message("[MIP][EtapeA] Erreur solveur : ", e$message)
        NULL
    })
    
    if (is.null(result_a)) {
        return(list(statut = "error", couverture_max_n = 0,
                    couverture_max_pct = 0, y_ref = NULL, x_ref = NULL))
    }
    
    sol_c <- get_solution(result_a, c[i])
    couverture_max_n <- sum(sol_c$value > 0.5)
    
    list(
        statut              = solver_status(result_a),
        couverture_max_n    = couverture_max_n,
        couverture_max_pct  = round(100 * couverture_max_n / m, 1),
        y_ref               = get_solution(result_a, y[j]),
        x_ref               = get_solution(result_a, x[i, j])
    )
}


#' Resoudre le LRP avec distances OSM (optimisation lexicographique)
#'
#' @param res_dist Resultat de matrice_distances_osm() (doit contenir graphe)
#' @param menages Data.frame des menages
#' @param nb_vehicules Nombre de vehicules disponibles (usage optionnel,
#'   pas tous forcement utilises)
#' @param cap_camion Capacite des camions
#' @param cap_point Capacite des points de collecte
#' @param d_max Distance max menage -> point consideree "couverte" (m)
#' @param n_max Nombre max de points ouverts (budget)
#' @param solveur Solveur a utiliser
#' @param time_limit Limite de temps pour l'etape B (routage complet)
#' @param time_limit_couverture Limite de temps pour l'etape A (couverture,
#'   modele plus leger -- defaut 60s)
#' @param tol_couverture Tolerance (en nombre de menages) sous la couverture
#'   maximale trouvee en etape A. 0 = on exige exactement la meilleure
#'   couverture possible. >0 = on accepte de sacrifier jusqu'a tol_couverture
#'   menages de couverture si cela reduit significativement le cout de
#'   tournee.
#' @return Solution du modele
#' @export
resoudre_lrp <- function(res_dist,
                         menages,
                         nb_vehicules,
                         cap_camion,
                         cap_point,
                         d_max,
                         n_max,
                         solveur               = "glpk",
                         time_limit            = 120,
                         time_limit_couverture = 60,
                         tol_couverture        = 0) {
    
    # -- Verifications --------------------------------------------
    if (is.null(res_dist)) stop("res_dist est NULL")
    if (is.null(menages) || nrow(menages) == 0) stop("Aucun menage fourni")
    
    n <- res_dist$n
    m <- nrow(menages)
    V <- nb_vehicules
    Qv <- cap_camion
    Cj <- rep(cap_point, n)
    W <- menages$poids_dechets
    
    O <- res_dist$idx_depot
    S <- res_dist$idx_decharge
    N_ext <- n + 2
    d <- res_dist$d_matrix
    cand_ext <- seq(2, n + 1)
    
    message("[MIP] Debut de l'optimisation (lexicographique, 2 etapes)")
    message("[MIP]   - Menages : ", m)
    message("[MIP]   - Points candidats : ", n)
    message("[MIP]   - Vehicules disponibles : ", V)
    message("[MIP]   - Capacite point : ", cap_point, " kg")
    message("[MIP]   - Capacite camion : ", cap_camion, " kg")
    message("[MIP]   - Distance de couverture (Dmax) : ", d_max, " m")
    message("[MIP]   - Nmax : ", n_max)
    message("[MIP]   - Tolerance couverture (etape B) : ", tol_couverture, " menage(s)")
    
    # -- Distances menages -> candidats -----------------------------
    message("[MIP] Calcul des distances d'affectation avec OSM...")
    d_men_cand <- calculer_distances_menages_candidats(res_dist, menages, cand_ext)
    
    couvrable_theorique <- (d_men_cand <= d_max)
    pct_couvrable_theorique <- 100 * sum(rowSums(couvrable_theorique) > 0) / m
    message("[MIP]   - Menages couvrables en theorie (<=", d_max,
            "m, tous candidats confondus, meme fermes) : ",
            round(pct_couvrable_theorique, 1), "%")
    
    # -- ETAPE A : couverture maximale atteignable -------------------
    message("[MIP] Etape A -- calcul de la couverture maximale sous Nmax...")
    res_a <- resoudre_couverture_max(d_men_cand, W, Cj, n_max, d_max,
                                     solveur = solveur,
                                     time_limit = time_limit_couverture)
    
    if (res_a$statut == "error") {
        return(list(statut = "error",
                    message = "Erreur lors de l'etape A (couverture maximale)"))
    }
    
    couverture_max_n   <- res_a$couverture_max_n
    couverture_max_pct <- res_a$couverture_max_pct
    message("[MIP]   Couverture maximale atteignable avec Nmax=", n_max, " : ",
            couverture_max_n, " / ", m, " menages (", couverture_max_pct, "%)")
    
    seuil_couverture <- max(couverture_max_n - tol_couverture, 0)
    
    # -- ETAPE B : routage optimal sous couverture garantie -----------
    message("[MIP] Etape B -- construction du modele complet (localisation + affectation + routage)...")
    
    model <- MIPModel() %>%
        add_variable(y[j],    j = 1:n,                              type = "binary") %>%
        add_variable(x[i, j], i = 1:m, j = 1:n,                    type = "binary") %>%
        add_variable(c[i],    i = 1:m,                              type = "binary") %>%
        add_variable(z[v, i, j], v = 1:V, i = 1:N_ext, j = 1:N_ext, type = "binary") %>%
        add_variable(u[v, j], v = 1:V, j = 1:N_ext,                 type = "continuous", lb = 0) %>%
        
        set_objective(
            sum_expr(d[i, j] * z[v, i, j], v = 1:V, i = 1:N_ext, j = 1:N_ext),
            sense = "min"
        ) %>%
        
        add_constraint(sum_expr(z[v, O, j], j = cand_ext) <= 1, v = 1:V) %>%
        add_constraint(sum_expr(z[v, i, S], i = cand_ext) <= 1, v = 1:V) %>%
        add_constraint(
            sum_expr(z[v, O, j], j = cand_ext) == sum_expr(z[v, i, S], i = cand_ext),
            v = 1:V
        ) %>%
        add_constraint(
            sum_expr(z[v, i, j], i = 1:N_ext) == sum_expr(z[v, j, k], k = 1:N_ext),
            v = 1:V, j = cand_ext
        ) %>%
        add_constraint(
            sum_expr(z[v, i, j], v = 1:V, i = 1:N_ext) == y[j - 1],
            j = cand_ext
        ) %>%
        add_constraint(z[v, i, i] == 0, v = 1:V, i = 1:N_ext) %>%
        add_constraint(z[v, S, j] == 0, v = 1:V, j = 1:N_ext) %>%
        add_constraint(z[v, i, O] == 0, v = 1:V, i = 1:N_ext) %>%
        
        add_constraint(sum_expr(x[i, j], j = 1:n) == 1, i = 1:m) %>%
        add_constraint(x[i, j] <= y[j], i = 1:m, j = 1:n) %>%
        
        add_constraint(
            sum_expr(W[i] * x[i, j], i = 1:m) <= Cj[j] * y[j],
            j = 1:n
        ) %>%
        add_constraint(sum_expr(y[j], j = 1:n) <= n_max) %>%
        
        add_constraint(
            c[i] <= sum_expr(x[i, j], j = 1:n, d_men_cand[i, j] <= d_max),
            i = 1:m
        ) %>%
        add_constraint(sum_expr(c[i], i = 1:m) >= seuil_couverture) %>%
        
        add_constraint(
            u[v, j] >= u[v, i] + Cj[j - 1] * z[v, i, j] - Qv * (1 - z[v, i, j]),
            v = 1:V, i = cand_ext, j = cand_ext
        ) %>%
        add_constraint(u[v, j] <= Qv, v = 1:V, j = 1:N_ext) %>%
        add_constraint(u[v, O] == 0, v = 1:V)
    
    # -- Resolution -----------------------------------------------
    message("[MIP] Resolution (etape B) avec ", solveur, "...")
    
    result <- tryCatch({
        solve_model(
            model,
            with_ROI(solver = solveur,
                     control = list(tm_limit = time_limit * 1000, verbose = TRUE))
        )
    }, error = function(e) {
        message("[MIP] Erreur solveur : ", e$message)
        NULL
    })
    
    if (is.null(result)) {
        return(list(statut = "error", message = "Erreur lors de la resolution (etape B)"))
    }
    
    statut   <- solver_status(result)
    cout_obj <- objective_value(result)
    
    message("[MIP] Resolution terminee - Statut : ", statut)
    message("[MIP]   - Distance totale des tournees : ", round(cout_obj, 2), " m")
    
    sol_y <- get_solution(result, y[j])
    sol_x <- get_solution(result, x[i, j])
    sol_z <- get_solution(result, z[v, i, j])
    sol_c <- get_solution(result, c[i])
    
    points_ouverts <- sol_y$j[sol_y$value > 0.5]
    arcs_actifs    <- sol_z[sol_z$value > 0.5, ]
    
    vehicules_utilises    <- sort(unique(arcs_actifs$v[arcs_actifs$i == O]))
    nb_vehicules_utilises <- length(vehicules_utilises)
    message("[MIP]   - Vehicules utilises : ", nb_vehicules_utilises, " / ", V,
            " | Points ouverts : ", length(points_ouverts))
    
    distances_par_vehicule <- sapply(1:V, function(v_id) {
        arcs_v <- arcs_actifs[arcs_actifs$v == v_id, ]
        if (nrow(arcs_v) == 0) return(0)
        sum(sapply(seq_len(nrow(arcs_v)), function(k) d[arcs_v$i[k], arcs_v$j[k]]))
    })
    names(distances_par_vehicule) <- paste0("Vehicule_", 1:V)
    
    couverture_pct <- round(100 * sum(sol_c$value > 0.5) / m, 1)
    
    sol_x_actifs <- sol_x[sol_x$value > 0.5, ]
    sol_x_actifs$distance_reelle_m <- mapply(
        function(ii, jj) d_men_cand[ii, jj], sol_x_actifs$i, sol_x_actifs$j
    )
    sol_x_actifs$couvert <- sol_x_actifs$distance_reelle_m <= d_max
    
    menages_non_couverts <- sol_x_actifs$i[!sol_x_actifs$couvert]
    distance_moy_non_couverts <- if (length(menages_non_couverts) > 0) {
        round(mean(sol_x_actifs$distance_reelle_m[!sol_x_actifs$couvert]), 1)
    } else 0
    
    message("[MIP]   - Couverture reelle obtenue : ", couverture_pct,
            "% (cible etape A : ", couverture_max_pct, "%, tolerance : ",
            tol_couverture, " menage(s))")
    message("[MIP]   - Menages non couverts : ", length(menages_non_couverts),
            " (distance moyenne : ", distance_moy_non_couverts, " m)")
    
    if (couverture_max_pct < 80) {
        warning("[MIP] La couverture maximale ATTEIGNABLE sous Nmax=", n_max,
                " est deja sous 80% (", couverture_max_pct,
                "%) -- ce n'est pas un defaut du solveur, c'est un budget ",
                "insuffisant. Augmenter Nmax pour ameliorer la couverture.")
    }
    
    tournees_detail <- arcs_actifs %>%
        dplyr::mutate(
            distance_km = mapply(function(ii, jj) round(d[ii, jj], 3), i, j),
            noeud_depart = dplyr::case_when(
                i == O ~ "Depot", i == S ~ "Decharge", TRUE ~ paste0("Point_", i - 1)
            ),
            noeud_arrivee = dplyr::case_when(
                j == O ~ "Depot", j == S ~ "Decharge", TRUE ~ paste0("Point_", j - 1)
            )
        ) %>%
        dplyr::select(Vehicule = v, De = noeud_depart, Vers = noeud_arrivee, Distance_km = distance_km) %>%
        dplyr::arrange(Vehicule, De)
    
    list(
        statut                     = statut,
        cout_total                 = round(cout_obj, 3),
        points_ouverts             = points_ouverts,
        sol_y                      = sol_y,
        sol_z                      = arcs_actifs,
        sol_x                      = sol_x_actifs,
        distances_par_vehicule     = round(distances_par_vehicule, 3),
        couverture                 = couverture_pct,
        couverture_max_atteignable = couverture_max_pct,
        tol_couverture_utilisee    = tol_couverture,
        tournees_detail            = tournees_detail,
        nb_vehicules               = V,
        nb_vehicules_utilises      = nb_vehicules_utilises,
        nb_points                  = length(points_ouverts),
        distances_affectation      = d_men_cand,
        pct_couvrable_theorique    = pct_couvrable_theorique,
        menages_non_couverts       = menages_non_couverts,
        nb_non_couverts            = length(menages_non_couverts),
        distance_moy_non_couverts  = distance_moy_non_couverts,
        sol_x_distances            = sol_x_actifs
    )
}
