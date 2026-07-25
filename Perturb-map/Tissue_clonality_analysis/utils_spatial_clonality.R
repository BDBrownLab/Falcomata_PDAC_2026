# Load required libraries
library(dplyr)
library(igraph)
library(dbscan)
library(data.table)

# Function to create kNN graph from spatial coordinates
knn_from_coord <- function(matrix_to_use, cell_names, k, distance_thresh = NULL) {
  names(cell_names) <- cell_names
  
  nn_network <- dbscan::kNN(
    x = matrix_to_use,
    k = k,
    sort = TRUE,
    search = "kdtree"
  )
  
  nn_network_dt <- data.table::data.table(
    from = rep(1:nrow(nn_network$id), k),
    to = as.vector(nn_network$id),
    weight = 1 / (1 + as.vector(nn_network$dist)),
    distance = as.vector(nn_network$dist)
  )
  
  nn_network_dt[, from_cell_ID := cell_names[from]]
  nn_network_dt[, to_cell_ID := cell_names[to]]
  
  if (!is.null(distance_thresh)) {
    nn_network_dt <- nn_network_dt[distance < distance_thresh]
  }
  
  all_index <- unique(c(nn_network_dt$from_cell_ID, nn_network_dt$to_cell_ID))
  
  igraph::graph_from_data_frame(
    nn_network_dt[, .(from_cell_ID, to_cell_ID, weight, distance)],
    directed = TRUE,
    vertices = all_index
  )
}

# Function to calculate the fraction of non-group neighbors
calc_nonGroupMembers <- function(in_graph, graph_id_vec) {
  pc_groups <- split(graph_id_vec, graph_id_vec)
  
  non_group_neighbors <- sapply(names(graph_id_vec), function(x) {
    setdiff(
      igraph::neighbors(in_graph, x, mode = "out") %>% names(),
      pc_groups[[graph_id_vec[x]]] %>% names()
    ) %>% length()
  })
  
  deg_out <- degree(in_graph, names(graph_id_vec), mode = "out")
  frac_non_group <- ifelse(deg_out > 0, non_group_neighbors / deg_out, NA_real_)
  
  data.frame(
    id = names(graph_id_vec),
    non_group_neighbors = non_group_neighbors,
    deg = deg_out,
    frac_non_group = frac_non_group
  )
}