library(sf)
library(dodgr)
library(dplyr)

# 1. Charger le fichier
routes <- st_read("data/routes_dakar.geojson")

# 2. Garder seulement les lignes
routes_lines <- routes[st_geometry_type(routes) %in% c("LINESTRING", "MULTILINESTRING"), ]

# 3. Voir les noms des colonnes pour identifier l'ID
names(routes_lines)

# 4. Convertir en graphe en spécifiant la colonne ID
# Utiliser "id" comme colonne d'identifiant
graphe <- weight_streetnet(
    routes_lines, 
    wt_profile = "motorcar",
    id_col = "id"  # Spécifier la colonne ID
)

# 5. Nettoyer
graphe <- graphe %>% filter(!is.na(d) & d > 0 & d < 50000)

# 6. Sauvegarder
saveRDS(graphe, file = "data/graphe_dakar.Rds")

message("✅ Arêtes : ", nrow(graphe))
message("📍 Nœuds : ", nrow(dodgr_vertices(graphe)))