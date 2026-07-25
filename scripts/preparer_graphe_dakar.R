# ============================================================
#  scripts/preparer_graphe_dakar.R
#  Préparation du graphe OSM pour Dakar
#  À exécuter UNE SEULE FOIS pour créer data/graphe_dakar.Rds
# ============================================================

# ── Packages ──────────────────────────────────────────────────
library(osmdata)
library(dodgr)
library(dplyr)
library(sf)

# ── Paramètres ──────────────────────────────────────────────
# Bbox de Dakar (étendue)
bbox_dakar <- c(-17.55, 14.60, -17.30, 14.80)
fichier_sortie <- "data/graphe_dakar.Rds"

# ── 1. Télécharger les données OSM ──────────────────────────
message("🌍 Téléchargement du réseau routier OSM pour Dakar...")
message("   Bbox : ", paste(bbox_dakar, collapse = ", "))

osm_data <- opq(bbox = bbox_dakar) %>%
    add_osm_feature(key = "highway") %>%
    osmdata_sf()

message("✅ Données OSM téléchargées")
message("   Points : ", nrow(osm_data$osm_points))
message("   Lignes : ", nrow(osm_data$osm_lines))
message("   Polygones : ", nrow(osm_data$osm_polygons))

# ── 2. Construction du graphe dodgr ─────────────────────────
message("🔧 Construction du graphe dodgr...")

graphe <- weight_streetnet(
    osm_data,
    wt_profile = "motorcar",
    type_col = "highway"
)

message("✅ Graphe construit")
message("   Arêtes : ", nrow(graphe))
message("   Colonnes : ", paste(names(graphe), collapse = ", "))

# ── 3. Nettoyage du graphe ──────────────────────────────────
message("🧹 Nettoyage du graphe...")

# Supprimer les arêtes avec distance NULL ou négative
graphe <- graphe %>% 
    filter(!is.na(d) & d > 0)

# Supprimer les arêtes trop longues (> 50 km = erreur)
graphe <- graphe %>% 
    filter(d < 50000)

message("✅ Nettoyage terminé")
message("   Arêtes restantes : ", nrow(graphe))

# ── 4. Ajout des sommets ─────────────────────────────────────
message("🔢 Ajout des sommets...")

# Extraire les sommets
vertices <- dodgr_vertices(graphe)
message("   Sommets : ", nrow(vertices))

# ── 5. Sauvegarde ────────────────────────────────────────────
message("💾 Sauvegarde du graphe...")

# Créer le dossier si nécessaire
if (!dir.exists("data")) {
    dir.create("data")
}

saveRDS(graphe, file = fichier_sortie)
message("✅ Graphe sauvegardé : ", fichier_sortie)
message("   Taille : ", round(file.size(fichier_sortie) / 1024 / 1024, 2), " Mo")

# ── 6. Vérification ──────────────────────────────────────────
message("🔍 Vérification...")
verif <- readRDS(fichier_sortie)
message("   Arêtes : ", nrow(verif))
message("   Sommets : ", nrow(dodgr_vertices(verif)))

message("\n✅ TERMINÉ !")
message("Pour utiliser ce graphe dans l'application :")
message("  1. Vérifiez que le fichier ", fichier_sortie, " existe")
message("  2. Lancez l'application Shiny")
message("  3. Le graphe sera automatiquement chargé")