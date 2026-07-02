# =============================================================================
# Code S1: Network Resilience Framework - Lake Chamo, Ethiopia
# Authors: Tadesse H, Atalaya K, Tesfaye G, Tesfaye G
# =============================================================================

# 1. Load packages
library(igraph)
library(terra)
library(sf)
library(ggplot2)
library(tidyterra)
library(rstac)
library(dplyr)
library(tidyr)
library(patchwork)
library(RColorBrewer)
library(ggspatial)
library(Hmisc)

# 2. Set working directory
setwd("C:/Users/pc/Desktop/New folder")

# 3. Load data
sites <- data.frame(
  Site = c("Segene", "Alfacho", "Open", "Sego", "Kulfo"),
  Longitude = c(37.626375, 37.595067, 37.562747, 37.504869, 37.560331),
  Latitude = c(5.853239, 5.794622, 5.833644, 5.836714, 5.929206),
  Resilience_Index = c(0.027, 1.000, 0.000, 0.718, 0.393)
)

# 4. Convert sites to sf
sites_sf <- st_as_sf(sites, coords = c("Longitude", "Latitude"), crs = 4326)

# 5. Network construction function
construct_network <- function(data, threshold = 0.5) {
  cor_matrix <- rcorr(as.matrix(data), type = "spearman")
  adj_matrix <- abs(cor_matrix$r) > threshold
  diag(adj_matrix) <- FALSE
  network <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
  metrics <- data.frame(
    Nodes = vcount(network),
    Edges = ecount(network),
    Density = edge_density(network),
    Avg_Degree = mean(degree(network)),
    Clustering = transitivity(network, type = "global"),
    Modularity = modularity(cluster_louvain(network))
  )
  return(list(network = network, metrics = metrics, cor_matrix = cor_matrix))
}

# 6. Keystone species function
identify_keystone <- function(network) {
  degree_cent <- degree(network)
  betweenness_cent <- betweenness(network)
  eigenvector_cent <- eigen_centrality(network)$vector
  degree_norm <- (degree_cent - min(degree_cent)) / (max(degree_cent) - min(degree_cent))
  betweenness_norm <- (betweenness_cent - min(betweenness_cent)) / (max(betweenness_cent) - min(betweenness_cent))
  eigenvector_norm <- (eigenvector_cent - min(eigenvector_cent)) / (max(eigenvector_cent) - min(eigenvector_cent))
  keystone_score <- (degree_norm + betweenness_norm + eigenvector_norm) / 3
  results <- data.frame(
    Variable = V(network)$name,
    Degree = degree_cent,
    Betweenness = betweenness_cent,
    Eigenvector = eigenvector_cent,
    Keystone_Score = keystone_score
  )
  return(results[order(results$Keystone_Score, decreasing = TRUE), ])
}

# 7. Eco-energy function
calculate_eco_energy <- function(data, cor_matrix, threshold = 0.5) {
  data_std <- scale(data)
  adj_matrix <- abs(cor_matrix$r) > threshold
  diag(adj_matrix) <- FALSE
  X <- as.matrix(data_std)
  A <- matrix(as.numeric(adj_matrix), nrow = nrow(adj_matrix), ncol = ncol(adj_matrix))
  energy <- -0.5 * t(X) %*% A %*% X
  return(sum(diag(energy)))
}

# 8. Satellite processing function
process_sentinel2 <- function(start_date, end_date, bbox, sites_utm) {
  s2 <- stac("https://earth-search.aws.element84.com/v1")
  date_str <- paste0(start_date, "T00:00:00Z/", end_date, "T23:59:59Z")
  results <- s2 |> stac_search(collections = "sentinel-2-l2a", bbox = bbox, datetime = date_str, limit = 20) |> post_request()
  if (length(results$features) == 0) {
    start <- as.Date(start_date) - 7
    end <- as.Date(end_date) + 7
    date_str <- paste0(start, "T00:00:00Z/", end, "T23:59:59Z")
    results <- s2 |> stac_search(collections = "sentinel-2-l2a", bbox = bbox, datetime = date_str, limit = 20) |> post_request()
  }
  if (length(results$features) == 0) return(rep(NA, nrow(sites_utm)))
  cloud_cover <- sapply(results$features, function(f) ifelse(!is.null(f$properties$`eo:cloud_cover`), f$properties$`eo:cloud_cover`, 100))
  best_idx <- which.min(cloud_cover)
  best_img <- results$features[[best_idx]]
  assets <- best_img$assets
  temp_red <- tempfile(fileext = ".tif")
  temp_nir <- tempfile(fileext = ".tif")
  download.file(assets$red$href, temp_red, mode = "wb", quiet = TRUE)
  download.file(assets$nir$href, temp_nir, mode = "wb", quiet = TRUE)
  band4 <- rast(temp_red)
  band8 <- rast(temp_nir)
  ndvi <- (band8 - band4) / (band8 + band4)
  ndvi_vals <- terra::extract(ndvi, sites_utm)[,2]
  file.remove(temp_red)
  file.remove(temp_nir)
  return(ndvi_vals)
}

# 9. Load NDVI and extract values
ndvi <- rast("lake_chamo_ndvi.tif")
sites_utm <- st_transform(sites_sf, crs(ndvi))
ndvi_values <- terra::extract(ndvi, sites_utm)
sites$NDVI <- ndvi_values[,2]

# 10. Process seasonal NDVI
seasons <- list(Dry = c("2024-01-15", "2024-01-15"), Short_Rainy = c("2024-02-15", "2024-02-15"), Transition = c("2024-04-15", "2024-04-15"), Long_Rainy = c("2023-08-21", "2023-08-21"))
bbox <- c(37.4, 5.7, 37.8, 6.0)
all_ndvi <- data.frame(Site = sites$Site)
for (season_name in names(seasons)) {
  dates <- seasons[[season_name]]
  all_ndvi[[paste0(season_name, "_NDVI")]] <- process_sentinel2(dates[1], dates[2], bbox, sites_utm)
}
all_ndvi$Resilience_Index <- sites$Resilience_Index

# 11. Statistics
ndvi_columns <- grep("_NDVI$", names(all_ndvi), value = TRUE)
all_ndvi$Avg_NDVI <- rowMeans(all_ndvi[, ndvi_columns], na.rm = TRUE)
cor_test <- cor.test(all_ndvi$Resilience_Index, all_ndvi$Avg_NDVI, use = "complete.obs")
lm_model <- lm(Avg_NDVI ~ Resilience_Index, data = all_ndvi)

# 12. Figures
algae_colors <- colorRampPalette(c("#08306B", "#08519C", "#2171B5", "#4292C6", "#6BAED6", "#9ECAE1", "#C6DBEF", "#EFF3FF", "#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45", "#006D2C"))(100)

p_validation <- ggplot(all_ndvi, aes(x = Resilience_Index, y = Avg_NDVI, label = Site)) + geom_point(size = 6, color = "darkblue") + geom_text(vjust = -0.8, hjust = 0.5, size = 5, fontface = "bold") + geom_smooth(method = "lm", se = TRUE, color = "red", alpha = 0.2) + labs(title = "Validation: Network Resilience vs Satellite NDVI", subtitle = "Lake Chamo, Ethiopia", x = "Resilience Index", y = "Average NDVI", caption = paste0("r = ", round(cor_test$estimate, 3), " | p = ", round(cor_test$p.value, 4), " | R² = ", round(summary(lm_model)$r.squared, 3))) + theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16), plot.subtitle = element_text(hjust = 0.5, size = 12), axis.title = element_text(size = 14), axis.text = element_text(size = 12), panel.grid.minor = element_blank())
ggsave("validation_plot.png", p_validation, width = 10, height = 8, dpi = 300)

time_series_long <- all_ndvi %>% select(-Avg_NDVI) %>% pivot_longer(cols = matches("_NDVI$"), names_to = "Season", values_to = "NDVI") %>% mutate(Season = factor(Season, levels = c("Dry_NDVI", "Short_Rainy_NDVI", "Transition_NDVI", "Long_Rainy_NDVI"), labels = c("Dry", "Short Rainy", "Transition", "Long Rainy")))
p_seasonal <- ggplot(time_series_long, aes(x = Season, y = NDVI, color = Site, group = Site)) + geom_line(size = 1.2) + geom_point(size = 4) + geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) + labs(title = "Seasonal NDVI Patterns at Lake Chamo Sites", x = "Season", y = "NDVI", color = "Site") + theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16), plot.subtitle = element_text(hjust = 0.5, size = 12), axis.title = element_text(size = 14), legend.position = "bottom")
ggsave("seasonal_ndvi_patterns.png", p_seasonal, width = 12, height = 8, dpi = 300)

p_combined <- (p_validation + p_seasonal) + plot_annotation(title = "Supplementary Figure S1: Satellite Validation of Network Resilience", subtitle = "Lake Chamo, Ethiopia", theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18), plot.subtitle = element_text(hjust = 0.5, size = 14)))
ggsave("supplementary_figure_s1.png", p_combined, width = 14, height = 10, dpi = 300)

# 13. Save outputs
summary_table <- all_ndvi %>% mutate(Interpretation = case_when(Avg_NDVI > 0 ~ "High algae", Avg_NDVI > -0.2 ~ "Moderate algae", TRUE ~ "Low algae")) %>% select(Site, Resilience_Index, Dry_NDVI, Short_Rainy_NDVI, Transition_NDVI, Long_Rainy_NDVI, Avg_NDVI, Interpretation)
write.csv(summary_table, "supplementary_table_s1.csv", row.names = FALSE)
write.csv(all_ndvi, "lake_chamo_satellite_validation_data.csv", row.names = FALSE)

print("✅ Analysis complete")