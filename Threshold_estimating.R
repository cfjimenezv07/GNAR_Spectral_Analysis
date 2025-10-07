# =========================
# === Threshld_estimating.R ===
# =========================
library(igraph)

# =========================
# === Helper Functions ===
# =========================

# Create adjacency matrix given a threshold
create_adj_matrix <- function(fevd_mat, threshold) {
  K <- nrow(fevd_mat)
  adj_mat <- matrix(0, K, K)
  for (i in 1:K) {
    for (j in 1:K) {
      if (i != j && (fevd_mat[i, j] >= threshold || fevd_mat[j, i] >= threshold)) {
        adj_mat[i, j] <- 1
        adj_mat[j, i] <- 1
      }
    }
  }
  diag(adj_mat) <- 0
  return(adj_mat)
}

# Compute network metrics for a range of thresholds
compute_network_metrics <- function(fevd_mat, thresholds) {
  K <- nrow(fevd_mat)
  num_components <- numeric(length(thresholds))
  largest_comp_size <- numeric(length(thresholds))
  connected_flag <- logical(length(thresholds))
  
  for (idx in seq_along(thresholds)) {
    adj_mat <- create_adj_matrix(fevd_mat, thresholds[idx])
    g <- graph_from_adjacency_matrix(adj_mat, mode = "undirected", diag = FALSE)
    comps <- components(g)
    
    num_components[idx] <- comps$no
    largest_comp_size[idx] <- max(comps$csize)
    connected_flag[idx] <- (comps$no == 1)
  }
  
  max_connected_threshold <- max(thresholds[connected_flag])
  
  list(
    thresholds = thresholds,
    num_components = num_components,
    largest_comp_size = largest_comp_size,
    connected_flag = connected_flag,
    max_connected_threshold = max_connected_threshold
  )
}

# Plot network metrics vs threshold
plot_network_metrics <- function(metrics, title_suffix = "") {
  # Plot Number of Connected Components
  plot(metrics$thresholds, metrics$num_components, type = "s",
       xlab = "Threshold", ylab = "Number of Connected Components",
       main = paste("Network Connectivity vs Threshold", title_suffix),
       ylim = c(1, max(metrics$num_components)))
  abline(v = metrics$max_connected_threshold, col = "red", lty = 2)
  text(metrics$max_connected_threshold, max(metrics$num_components),
       labels = paste0("Critical Threshold = ", round(metrics$max_connected_threshold, 4)),
       pos = 4, col = "red")
  
  # Plot Largest Component Size
  K <- max(metrics$largest_comp_size)
  plot(metrics$thresholds, metrics$largest_comp_size, type = "s",
       xlab = "Threshold", ylab = "Size of Largest Connected Component",
       main = paste("Largest Component Size vs Threshold", title_suffix),
       ylim = c(1, K))
  abline(v = metrics$max_connected_threshold, col = "red", lty = 2)
  text(metrics$max_connected_threshold, K,
       labels = paste0("Critical Threshold = ", round(metrics$max_connected_threshold, 4)),
       pos = 4, col = "red")
}

# =========================
# === 1. Full Network ===
# =========================
thresholds <- seq(0, max(fevd), length.out = 100)
metrics_full <- compute_network_metrics(fevd, thresholds)
cat("Largest threshold where the network is still connected:", metrics_full$max_connected_threshold, "\n")
plot_network_metrics(metrics_full)

# =========================
# === 2. Filtered Network (no .jp/.cn) ===
# =========================
thresholds_nojp_cn <- seq(0, max(fevd_nojp_cn), length.out = 100)
metrics_nojp_cn <- compute_network_metrics(fevd_nojp_cn, thresholds_nojp_cn)
cat("Largest threshold where the network (no .jp/.cn) is still connected:", metrics_nojp_cn$max_connected_threshold, "\n")
plot_network_metrics(metrics_nojp_cn, title_suffix = "(excluding .jp/.cn)")
