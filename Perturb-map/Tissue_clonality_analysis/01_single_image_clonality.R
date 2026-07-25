# Load required libraries
library(ggplot2)
library(dplyr)
library(viridis)
library(igraph)
library(dbscan)
library(meanShiftR)
library(data.table)
library(sf)
library(tidyr)

# Load helper functions
source("../R/utils_spatial_clonality.R")

### Step 1: Load and Prepare Data ###

# Load image data
img <- read.csv("../data/example.csv", header = TRUE, stringsAsFactors = FALSE)

# Add a unique cell ID
img <- img %>%
  mutate(cell_id = as.character(1:nrow(img)))

# Filter for cells with valid gene.filtered
img_filtered <- img %>%
  filter(!is.na(gene.filtered))

### Step 2: Fraction of Non-Group Neighbors ###

# Create kNN graph
test.graph <- knn_from_coord(
  matrix_to_use = img_filtered %>% select(x, y) %>% as.matrix(),
  cell_names = img_filtered$cell_id,
  k = 50,
  distance_thresh = 75
)

# Calculate fraction of non-group neighbors
graph_id_vec <- img_filtered$gene.filtered
names(graph_id_vec) <- img_filtered$cell_id
graph_id_vec <- graph_id_vec[V(test.graph)$name]

ngdf <- calc_nonGroupMembers(test.graph, graph_id_vec)

# Merge non-group neighbor data with filtered image data
img_filtered <- left_join(img_filtered, ngdf, by = c("cell_id" = "id"))

# Plot tissue clonality
g1 <- ggplot(img_filtered, aes(x = x, y = y)) +
  geom_point(
    data = img_filtered %>% filter(!is.na(gene.filtered)),
    aes(color = frac_non_group),
    size = 1,
    alpha = 0.4
  ) +
  theme_minimal() +
  scale_color_viridis(name = "% Non-Matching Neighbors", option = "A") +
  ggtitle("Tissue Clonality via Fraction of Non-Group Neighbors")

print(g1)
ggsave("../results/plot.pdf", g1, width = 9, height = 4, dpi = 300)

### Step 3: Focal Point Distributions ###

# Identify dense areas using mean shift clustering
h <- 250

ms_results <- meanShift(
  as.matrix(img_filtered %>% select(x, y)),
  as.matrix(img_filtered %>% select(x, y)),
  algorithm = "LINEAR",
  kernelType = "NORMAL",
  bandwidth = c(h, h),
  alpha = 0,
  iterations = 1000
)

img_filtered <- img_filtered %>%
  mutate(focal_point = ms_results$assignment[, 1])

# Filter clusters with >= 50 cells
cluster_sizes <- table(img_filtered$focal_point)
filtered_focal_points <- names(cluster_sizes[cluster_sizes >= 50])

img_filtered <- img_filtered %>%
  filter(focal_point %in% filtered_focal_points)

# Compute centroids
fp_centroid <- img_filtered %>%
  group_by(focal_point) %>%
  summarise(
    x_centroid = mean(x),
    y_centroid = mean(y),
    .groups = "drop"
  )

# Define a window size for square polygons around each centroid
window_scale <- 200

# Generate square polygons around each focal point centroid
fp_centroid_window <- lapply(1:nrow(fp_centroid), function(i) {
  centroid <- fp_centroid[i, c("x_centroid", "y_centroid")]
  
  matrix_coords <- matrix(
    c(
      centroid$x_centroid + window_scale, centroid$y_centroid + window_scale,
      centroid$x_centroid + window_scale, centroid$y_centroid - window_scale,
      centroid$x_centroid - window_scale, centroid$y_centroid - window_scale,
      centroid$x_centroid - window_scale, centroid$y_centroid + window_scale,
      centroid$x_centroid + window_scale, centroid$y_centroid + window_scale
    ),
    ncol = 2,
    byrow = TRUE
  )
  
  st_polygon(list(matrix_coords))
})

names(fp_centroid_window) <- fp_centroid$focal_point

# Ensure focal_point is a factor
img_filtered <- img_filtered %>%
  mutate(focal_point = as.factor(focal_point))

# Visualize density and focal points
g2 <- ggplot(img_filtered, aes(x = x, y = y)) +
  geom_density_2d_filled(aes(fill = after_stat(level)), alpha = 0.5) +
  geom_point(aes(color = focal_point), size = 0.5, alpha = 0.7) +
  geom_text(
    data = fp_centroid,
    aes(x = x_centroid, y = y_centroid, label = focal_point),
    color = "white",
    size = 6,
    fontface = "bold"
  ) +
  scale_fill_manual(
    name = "Focal Points",
    values = c(
      "#fff5f0", "#fee0d2", "#fcbba1", "#fc9272",
               "#fb6a4a", "#ef3b2c", "#cb181d", "#a50f15",
               "#67000d", "#67000d", "#67000d", "#67000d",
               "#67000d", "#67000d", "#67000d", "#67000d"
    )
  ) +
  scale_color_manual(
    name = "Focal Points",
    values = rep("grey40", 18)
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank()
  ) +
  coord_fixed() +
  ggtitle("Focal Point Kernel Density with Polygons")

# Add polygons for focal points
for (i in seq_along(fp_centroid_window)) {
  bound.coord <- st_coordinates(fp_centroid_window[[i]]) %>%
    as.data.frame() %>%
    dplyr::rename(x = X, y = Y)
  
  g2 <- g2 + geom_polygon(
    data = bound.coord,
    aes(x = x, y = y),
    color = "white",
    fill = NA,
    linewidth = 1
  )
}

print(g2)
ggsave("../results/plot2.pdf", g2, width = 8, height = 5, dpi = 300)

### Step 4: Pro-Code Frequencies in Focal Points ###

filtered_focal_points <- unique(img_filtered$focal_point)

fp_cells <- lapply(filtered_focal_points, function(fp) {
  img_filtered %>%
    filter(focal_point == fp) %>%
    count(gene.filtered) %>%
    mutate(
      freq = n / sum(n),
      focal_point = as.character(fp)
    )
})

fp_freq <- bind_rows(fp_cells) %>%
  mutate(
    focal_point = factor(focal_point, levels = as.character(filtered_focal_points))
  )

write.csv(fp_freq, "../results/focal_point_freq.csv", row.names = FALSE)

# Stacked bar plot of pro-code frequencies
bp <- ggplot(fp_freq, aes(x = focal_point, y = freq, fill = gene.filtered)) +
  geom_bar(stat = "identity", color = "black") +
  theme_classic() +
  scale_fill_viridis(discrete = TRUE, na.value = "grey", option = "A") +
  ggtitle("Pro-Code Frequencies at Focal Points")

print(bp)
ggsave("../results/plot3.pdf", bp, width = 7, height = 4, dpi = 300)
