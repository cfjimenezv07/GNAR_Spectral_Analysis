################################################################################
#                          GNAR SPECTRUM ANALYSIS PIPELINE
# ==============================================================================
# Title:     Soft Thresholding of Spectral Matrices on Networks
# Author:    Cristian Felipe Jimenez Varon
# Purpose:   (1) Apply soft thresholding to estimated spectral matrices
#            (2) Compute GNAR spectrum, coherence, and partial coherence RMSEs
# Date:      October 2025
################################################################################

# ==============================================================================
#                                SETUP
# ==============================================================================

# --- Define Directories ---
dir.pc  <- "G:/My Drive/2024/York/"
dir.mac <- "~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/2024/York/"
dir     <- dir.mac  # change to dir.pc if running on PC

setwd(paste0(dir, "Rcodes"))
dir.r <- paste0(dir, "Results/")

# --- Load Required Libraries ---
library(igraph)
library(abind)
library(GNAR)
library(fields)
library(signal)
library(qfa)

# --- Define Method Names ---
method_names <- c(
  "GNAR_est", "VAR_est", "VAR_par_A", "VAR_par_A1",
  "SDM_est_A", "SDM_est_A1", "GNAR_np"
)

# ==============================================================================
#                         UTILITY FUNCTIONS
# ==============================================================================

# --- Soft Thresholding for Complex Numbers ---
soft_thresh <- function(z, rho) {
  magnitude <- Mod(z)
  angle <- Arg(z)
  thresholded_mag <- pmax(magnitude - rho, 0)
  thresholded_mag * exp(1i * angle)
}

# --- Generate a Connected GNAR Network ---
generate_gnar_network <- function(n_nodes, p = 0.2) {
  g <- igraph::sample_gnp(n = n_nodes, p = p, directed = FALSE)
  while (!igraph::is_connected(g)) {
    g <- igraph::sample_gnp(n = n_nodes, p = p, directed = FALSE)
  }
  adj_matrix <- as.matrix(igraph::as_adjacency_matrix(g))
  GNAR::as.GNARnet(adj_matrix)
}

# --- Compute Pairwise Graph Distances ---
compute_graph_distances <- function(adj_matrix) {
  G <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
  igraph::distances(G)
}

# --- Invert a 3D Array of Spectral Matrices (f ¡ú S) ---
invert_spectrum <- function(spectrum) {
  array(apply(spectrum, 3, function(fk) solve(fk)), dim = dim(spectrum))
}

# --- Compute Threshold Values per Frequency ---
compute_thresholds_per_freq <- function(S, graph_dists, r_star) {
  n_freq <- dim(S)[3]
  xi_list <- vector("list", r_star + 2)
  xi_list[[1]] <- rep(Inf, n_freq)  # For r = 0
  
  for (r in 1:r_star) {
    dist_mask <- graph_dists %in% c(2 * r - 1, 2 * r)
    dist_indices <- arrayInd(which(dist_mask), dim(graph_dists))
    
    xi_r <- sapply(1:n_freq, function(k) {
      vals <- abs(S[,,k])[dist_indices]
      if (length(vals) == 0) NA else min(vals)
    })
    xi_list[[r + 1]] <- xi_r
  }
  
  xi_list[[r_star + 2]] <- rep(0, n_freq)
  xi_list
}

# --- Apply Soft-Thresholding to Inverse Spectrum ---
threshold_inverse_spectrum <- function(S, xi_list, r_star) {
  n <- dim(S)[1]
  n_freq <- dim(S)[3]
  result <- vector("list", r_star)
  names(result) <- paste0("r=", 1:r_star)
  
  for (r in 1:r_star) {
    S_r <- array(0 + 0i, dim = c(n, n, n_freq))
    for (k in 1:n_freq) {
      Sk <- S[,,k]
      diag_adjusted <- Sk + diag(xi_list[[r + 1]][k] * exp(1i * Arg(diag(Sk))))
      Sk_thresh <- matrix(0 + 0i, n, n)
      for (i in 1:n) {
        for (j in 1:n) {
          Sk_thresh[i, j] <- soft_thresh(diag_adjusted[i, j], xi_list[[r + 1]][k])
        }
      }
      S_r[,,k] <- Sk_thresh
    }
    result[[r]] <- S_r
  }
  result
}

# --- Invert Thresholded Spectra (S_r ¡ú f_r) ---
invert_thresholded_spectra <- function(thresholded_spectra) {
  lapply(thresholded_spectra, function(S_r) {
    array(apply(S_r, 3, function(mat) solve(mat)), dim = dim(S_r))
  })
}

# --- Process a Single Spectrum ---
process_one_spectrum <- function(spectrum, adj_matrix, r_star = 3) {
  graph_dists <- compute_graph_distances(adj_matrix)
  S <- invert_spectrum(spectrum)
  xi_list <- compute_thresholds_per_freq(S, graph_dists, r_star)
  thresholded_spectra <- threshold_inverse_spectrum(S, xi_list, r_star)
  f_r_array <- invert_thresholded_spectra(thresholded_spectra)
  list(S_r = thresholded_spectra, f_r = f_r_array)
}

# ==============================================================================
#                      PROCESS ESTIMATED SPECTRA
# ==============================================================================

process_all_methods <- function(method_list, adj_matrix,
                                method_names = NULL,
                                T_values = c(100, 200, 500, 1000),
                                dir.r = dir.r,
                                verbose = TRUE) {
  
  M <- length(method_list)
  if (is.null(method_names)) method_names <- paste0("method_", seq_len(M))
  if (length(method_names) != M) stop("Mismatch in method_names and method_list length")
  
  dir.create(dir.r, showWarnings = FALSE, recursive = TRUE)
  result <- setNames(vector("list", M), method_names)
  
  for (m in seq_len(M)) {
    method_name <- method_names[m]
    if (verbose) cat("\n>>> Processing method:", method_name, "\n")
    method_data <- method_list[[m]]
    
    result[[m]] <- lapply(seq_along(method_data), function(i) {
      T_value <- T_values[i]
      lapply(seq_along(method_data[[i]]), function(j) {
        r_star_j <- if (j == 1) 1 else if (j == 2) 2 else 3
        
        lapply(seq_along(method_data[[i]][[j]]), function(k) {
          if (verbose)
            cat(sprintf("Method: %s ~ T: %d ~ Model: %d ~ Sim: %d ~ r* = %d\n",
                        method_name, T_value, j, k, r_star_j))
          spectrum <- method_data[[i]][[j]][[k]]
          process_one_spectrum(spectrum, adj_matrix, r_star = r_star_j)
        })
      })
    })
    
    saveRDS(result[[m]], file = file.path(dir.r, paste0("sft_spectrum_", method_name, ".rds")))
  }
  return(result)
}

# ==============================================================================
#                       PROCESS TRUE SPECTRUM
# ==============================================================================

process_true_spectrum <- function(true_spec_list, adj_matrix,
                                  T_values = c(100, 200, 500, 1000),
                                  dir.r = "./results/true_spectra",
                                  verbose = TRUE) {
  
  stopifnot(length(true_spec_list) == length(T_values))
  dir.create(dir.r, showWarnings = FALSE, recursive = TRUE)
  
  result <- vector("list", length(T_values))
  names(result) <- paste0("T=", T_values)
  
  for (i in seq_along(true_spec_list)) {
    T_val <- T_values[i]
    model_list <- true_spec_list[[i]]
    
    if (verbose) cat("\n>>> Processing T =", T_val, "\n")
    
    result[[i]] <- lapply(seq_along(model_list), function(j) {
      r_star <- if (j == 1) 1 else if (j == 2) 2 else 3
      if (verbose) cat(sprintf("Model: %d (r* = %d)\n", j, r_star))
      
      spectrum <- model_list[[j]]
      process_one_spectrum(spectrum, adj_matrix, r_star = r_star)
    })
  }
  
  saveRDS(result, file = file.path(dir.r, "true_sft_spectrum.rds"))
  return(result)
}

# ==============================================================================
#                              MAIN EXECUTION
# ==============================================================================

set.seed(42)
gnar_network <- generate_gnar_network(n_nodes = 10, p = 0.05)
adj_matrix <- as.matrix(gnar_network)

estimated_gnar_spec_10 <- readRDS(file.path(dir, "Results/estimated_gnar_spec_10.rds"))
method_list <- list(estimated_gnar_spec_10)
T_values <- c(100, 200, 500, 1000)

result_est <- process_all_methods(
  method_list = method_list,
  adj_matrix = adj_matrix,
  method_names = c("GNAR_est"),
  T_values = T_values,
  dir.r = dir.r,
  verbose = TRUE
)

True_gnar_spec <- readRDS(file.path(dir, "Results/True_gnar_spec.rds"))
True <- process_true_spectrum(
  true_spec_list = True_gnar_spec,
  adj_matrix = adj_matrix,
  T_values = T_values,
  dir.r = dir.r,
  verbose = TRUE
)

# ==============================================================================
#              RMSE EVALUATION FOR GNAR SPECTRUM, COHERENCE, AND PCOH
# ==============================================================================

library(abind)

# --- Load Processed Spectra ---
sft_spec_gnar_est <- readRDS(file.path(dir.r, "sft_spectrum_GNAR_est.rds"))
True_gnar_spec <- readRDS(file.path(dir, "Results/True_gnar_spec.rds"))

# --- Parameters ---
T_values      <- c(100, 200, 500, 1000)
n_datasets    <- length(T_values)
n_simulations <- 500
n_models      <- length(True_gnar_spec[[1]])
r_star        <- 3
n_nodes       <- dim(True_gnar_spec[[1]][[1]])[1]

# --- Initialize Arrays ---
rmse_array_common      <- array(NA, dim = c(n_models, n_datasets, r_star))
coh_rmse_array_common  <- array(NA, dim = c(n_models, n_datasets, r_star))
pcoh_rmse_array_common <- array(NA, dim = c(n_models, n_datasets, r_star))

# --- Helper Function ---
get_fixed_r <- function(model_idx) {
  if (model_idx == 1) return(1)
  if (model_idx == 2) return(2)
  return(3)
}

# --- Main Loop ---
for (dataset_idx in seq_len(n_datasets)) {
  cat("Dataset T =", T_values[dataset_idx], "\n")
  for (model_idx in seq_len(n_models)) {
    cat("  Model", model_idx, "\n")
    fixed_r <- get_fixed_r(model_idx)
    
    true_spec     <- True_gnar_spec[[dataset_idx]][[model_idx]]
    n_freq        <- dim(true_spec)[3]
    true_inv_spec <- array(apply(true_spec, 3, solve), dim = dim(true_spec))
    
    for (r_level in seq_len(fixed_r)) {
      cat("    r =", r_level, "\n")
      rmse_sim <- coh_sim <- pcoh_sim <- numeric(n_simulations)
      
      for (sim_idx in seq_len(n_simulations)) {
        est_spectrum <- sft_spec_gnar_est[[dataset_idx]][[model_idx]][[sim_idx]][[2]][[r_level]]
        est_inv_spec <- sft_spec_gnar_est[[dataset_idx]][[model_idx]][[sim_idx]][[1]][[r_level]]
        
        # --- RMSE (Spectrum) ---
        rmse_sim[sim_idx] <- mean(sapply(1:n_freq, function(f)
          sqrt(mean(Mod(est_spectrum[,,f] - true_spec[,,f])^2))
        ))
        
        # --- RMSE (Coherence) ---
        coh_sim[sim_idx] <- mean(sapply(1:n_freq, function(f) {
          coh_true <- coh_est <- matrix(1, n_nodes, n_nodes)
          for (i in 1:n_nodes) for (j in 1:n_nodes)
            if (i != j) {
              Sii_t <- Re(true_spec[i, i, f]); Sjj_t <- Re(true_spec[j, j, f])
              Sii_e <- Re(est_spectrum[i, i, f]); Sjj_e <- Re(est_spectrum[j, j, f])
              coh_true[i, j] <- Mod(true_spec[i, j, f])^2 / (Sii_t * Sjj_t)
              coh_est[i, j]  <- Mod(est_spectrum[i, j, f])^2 / (Sii_e * Sjj_e)
            }
          sqrt(mean((coh_est - coh_true)^2))
        }))
        
        # --- RMSE (Partial Coherence) ---
        pcoh_sim[sim_idx] <- mean(sapply(1:n_freq, function(f) {
          G_true <- true_inv_spec[,,f]
          G_est  <- est_inv_spec[,,f]
          D_true <- diag(1 / sqrt(abs(Re(diag(G_true)))))
          D_est  <- diag(1 / sqrt(abs(Re(diag(G_est)))))
          Gamma_true <- -D_true %*% G_true %*% D_true
          Gamma_est  <- -D_est  %*% G_est  %*% D_est
          sqrt(mean((Mod(Gamma_est)^2 - Mod(Gamma_true)^2)^2))
        }), na.rm = TRUE)
      }
      
      # --- Average over Simulations ---
      rmse_array_common[model_idx, dataset_idx, r_level]      <- mean(rmse_sim, na.rm = TRUE)
      coh_rmse_array_common[model_idx, dataset_idx, r_level]  <- mean(coh_sim, na.rm = TRUE)
      pcoh_rmse_array_common[model_idx, dataset_idx, r_level] <- mean(pcoh_sim, na.rm = TRUE)
    }
  }
}

# ==============================================================================
#                     FORMAT RESULTS INTO TABLES
# ==============================================================================

reshape_rmse_array <- function(rmse_array, n_models, n_datasets, r_star, T_values) {
  rmse_matrix <- matrix(NA, nrow = n_models, ncol = n_datasets * r_star)
  col_names <- character(n_datasets * r_star)
  col_idx <- 1
  for (d in seq_len(n_datasets)) {
    for (r in seq_len(r_star)) {
      rmse_matrix[, col_idx] <- round(rmse_array[, d, r] * 100, 2)
      col_names[col_idx] <- paste0("T", T_values[d], "_r", r)
      col_idx <- col_idx + 1
    }
  }
  colnames(rmse_matrix) <- col_names
  rownames(rmse_matrix) <- paste0("Model_", seq_len(n_models))
  as.data.frame(rmse_matrix)
}

rmse_table_common      <- reshape_rmse_array(rmse_array_common,      n_models, n_datasets, r_star, T_values)
coh_rmse_table_common  <- reshape_rmse_array(coh_rmse_array_common,  n_models, n_datasets, r_star, T_values)
pcoh_rmse_table_common <- reshape_rmse_array(pcoh_rmse_array_common, n_models, n_datasets, r_star, T_values)

cat("=== RMSE Table (Spectrum) ===\n"); print(rmse_table_common)
cat("\n=== RMSE Table (Coherence) ===\n"); print(coh_rmse_table_common)
cat("\n=== RMSE Table (Partial Coherence) ===\n"); print(pcoh_rmse_table_common)

# ==============================================================================
#                             SAVE RESULTS
# ==============================================================================

saveRDS(
  list(
    rmse_array_common,
    coh_rmse_array_common,
    pcoh_rmse_array_common
  ),
  paste0(dir.r, "sft_rmse_results_common_true_GNAR.rds")
)

################################################################################
#                                    END
################################################################################
