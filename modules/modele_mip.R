# ============================================================
#  modules/modele_mip.R
#  Formulation et résolution du modèle LRP
# ============================================================

library(ompr)
library(ompr.roi)
library(ROI.plugin.glpk)

resoudre_lrp <- function(res_dist,
                         menages,
                         nb_vehicules,
                         cap_camion,
                         cap_point,
                         d_max,
                         n_max,
                         solveur    = "glpk",
                         time_limit = 120) {
    
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
    
    coords_cand <- res_dist$coords_etendues[cand_ext, c("longitude","latitude")]
    coords_men <- menages[, c("longitude","latitude")]
    
    d_men_cand <- geosphere::distm(
        as.matrix(coords_men),
        as.matrix(coords_cand),
        fun = geosphere::distHaversine
    )
    autorise_xij <- (d_men_cand <= d_max)
    
    model <- MIPModel() %>%
        add_variable(y[j],    j = 1:n,                    type = "binary") %>%
        add_variable(x[i, j], i = 1:m, j = 1:n,          type = "binary") %>%
        add_variable(z[v, i, j], v = 1:V, i = 1:N_ext, j = 1:N_ext, type = "binary") %>%
        add_variable(u[v, j], v = 1:V, j = 1:N_ext,      type = "continuous", lb = 0) %>%
        
        set_objective(
            sum_expr(d[i, j] * z[v, i, j], v = 1:V, i = 1:N_ext, j = 1:N_ext),
            sense = "min"
        ) %>%
        
        add_constraint(sum_expr(z[v, O, j], j = cand_ext) == 1, v = 1:V) %>%
        add_constraint(sum_expr(z[v, i, S], i = cand_ext) == 1, v = 1:V) %>%
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
        add_constraint(x[i, j] <= as.integer(autorise_xij[i, j]), i = 1:m, j = 1:n) %>%
        add_constraint(
            sum_expr(W[i] * x[i, j], i = 1:m) <= Cj[j] * y[j],
            j = 1:n
        ) %>%
        add_constraint(sum_expr(y[j], j = 1:n) <= n_max) %>%
        add_constraint(
            u[v, j] >= u[v, i] + Cj[j - 1] * z[v, i, j] - Qv * (1 - z[v, i, j]),
            v = 1:V, i = cand_ext, j = cand_ext
        ) %>%
        add_constraint(u[v, j] <= Qv, v = 1:V, j = 1:N_ext) %>%
        add_constraint(u[v, O] == 0, v = 1:V)
    
    result <- solve_model(
        model,
        with_ROI(
            solver = solveur,
            control = list(tm_limit = time_limit * 1000, verbose = TRUE)
        )
    )
    
    statut <- solver_status(result)
    cout_obj <- objective_value(result)
    
    sol_y <- get_solution(result, y[j])
    sol_x <- get_solution(result, x[i, j])
    sol_z <- get_solution(result, z[v, i, j])
    
    points_ouverts <- sol_y$j[sol_y$value > 0.5]
    arcs_actifs <- sol_z[sol_z$value > 0.5, ]
    
    distances_par_vehicule <- sapply(1:V, function(v_id) {
        arcs_v <- arcs_actifs[arcs_actifs$v == v_id, ]
        if (nrow(arcs_v) == 0) return(0)
        sum(sapply(seq_len(nrow(arcs_v)), function(k) {
            d[arcs_v$i[k], arcs_v$j[k]]
        }))
    })
    names(distances_par_vehicule) <- paste0("Vehicule_", 1:V)
    
    affectes <- unique(sol_x$i[sol_x$value > 0.5])
    couverture_pct <- round(100 * length(affectes) / m, 1)
    
    tournees_detail <- arcs_actifs %>%
        dplyr::mutate(
            distance_km = mapply(function(ii, jj) round(d[ii, jj], 3), i, j),
            noeud_depart = dplyr::case_when(
                i == O ~ "Dépôt",
                i == S ~ "Décharge",
                TRUE ~ paste0("Point_", i - 1)
            ),
            noeud_arrivee = dplyr::case_when(
                j == O ~ "Dépôt",
                j == S ~ "Décharge",
                TRUE ~ paste0("Point_", j - 1)
            )
        ) %>%
        dplyr::select(Vehicule = v, De = noeud_depart, Vers = noeud_arrivee, Distance_km = distance_km) %>%
        dplyr::arrange(Vehicule, De)
    
    list(
        statut = statut,
        cout_total = round(cout_obj, 3),
        points_ouverts = points_ouverts,
        sol_y = sol_y,
        sol_z = arcs_actifs,
        sol_x = sol_x[sol_x$value > 0.5, ],
        distances_par_vehicule = round(distances_par_vehicule, 3),
        couverture = couverture_pct,
        tournees_detail = tournees_detail,
        nb_vehicules = V,
        nb_points = length(points_ouverts)
    )
}