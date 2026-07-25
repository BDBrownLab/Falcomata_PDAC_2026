# Load required libraries
library(ggplot2)
library(dplyr)
library(viridis)
library(igraph)
library(dbscan)
library(meanShiftR)
library(data.table)

# Load helper functions
source("../R/utils_spatial_clonality.R")

# Load metadata
metadata <- read.csv("../metadata/sample_metadata.csv", stringsAsFactors = FALSE)

# Define function to process a single image
process_image <- function(image_path, sample_id) {
  # Load image data
  img <- read.csv(image_path, header = TRUE, stringsAsFactors = FALSE)
  
  # Add a unique cell ID
  img <- img %>%
    mutate(cell_id = as.character(1:nrow(img)))
  
  # Filter for cells with valid gene.filtered
  img_filtered <- img %>%
    filter(!is.na(gene.filtered))
  
  # kNN graph creation
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
  
  # Mean shift clustering for focal points
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
  
  # Calculate Shannon index and evenness for each focal point
  shannon_data <- img_filtered %>%
    group_by(focal_point) %>%
    count(gene.filtered) %>%
    mutate(freq = n / sum(n)) %>%
    summarise(
      shannon_index = -sum(freq * log(freq), na.rm = TRUE),
      evenness = shannon_index / log(n_distinct(gene.filtered)),
      .groups = "drop"
    )
  
  # Add sample_id
  shannon_data <- shannon_data %>%
    mutate(sample_id = sample_id)
  
  return(shannon_data)
}

# Apply the processing function to all images
image_folder <- "../../results/"
all_images <- list.files(image_folder, full.names = TRUE, pattern = "\\.csv$")

# Extract sample IDs from file names
sample_ids <- basename(all_images) %>% gsub("\\.csv$", "", .)

# Process each image
shannon_results <- lapply(seq_along(all_images), function(i) {
  process_image(all_images[i], sample_ids[i])
}) %>%
  bind_rows()

# Join Shannon index data with metadata
shannon_metadata <- left_join(shannon_results, metadata, by = c("sample_id" = "sample_id"))
write.csv(shannon_metadata, "../results/sample_shannon.csv", row.names = FALSE)

# Calculate median Shannon index per mouse
mouse_shannon <- shannon_metadata %>%
  group_by(mouse_day) %>%
  summarise(
    median_shannon = median(shannon_index, na.rm = TRUE),
    median_evenness = median(evenness, na.rm = TRUE),
    .groups = "drop"
  )

# Extract day from mouse_day
mouse_shannon <- mouse_shannon %>%
  mutate(day = sub(".*_(D\\d+)_.*", "\\1", mouse_day))

print(mouse_shannon)

write.csv(mouse_shannon, "../results/mouse_shannon.csv", row.names = FALSE)

# Compare median Shannon diversity by time point
p1 <- ggplot(mouse_shannon, aes(x = day, y = median_shannon, fill = day)) +
  geom_boxplot(outlier.color = "red", outlier.size = 1) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_minimal() +
  scale_fill_viridis(discrete = TRUE, name = "Time Point") +
  labs(
    title = "Comparison of Median Shannon Diversity Index by Time Point",
    x = "Time Point",
    y = "Median Shannon Index"
  )

print(p1)

p2 <- ggplot(mouse_shannon, aes(x = day, y = median_evenness, fill = day)) +
  geom_boxplot(outlier.color = "red", outlier.size = 1) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_minimal() +
  scale_fill_viridis(discrete = TRUE, name = "Time Point") +
  labs(
    title = "Comparison of Median Evenness by Time Point",
    x = "Time Point",
    y = "Median Evenness"
  )

print(p2)
