################################################################################
# Simulation 3: GNAR Spectrum Estimation under Model Misspecification
# Network size: 10 nodes
# Author: Cristian Felipe Jimenez Varon
# Description:
#   This script evaluates the GNAR spectral estimator under model misspecification
#   for a 10-node network. For each simulated dataset, the GNAR model is fitted
#   with incorrect order specifications, and estimation accuracy is evaluated via
#   RMSE metrics for the spectral density, coherence, and partial coherence.
################################################################################

# ==============================================================================
# 1. Setup
# ==============================================================================

dir.pc  <- "G:/My Drive/2024/York/"
dir.mac <- "~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/2024/York/"
dir     <- dir.pc  # Set directory depending on system

setwd(paste0(dir, "Rcodes"))
dir.r <- paste0(dir, "Results/")

# Load required libraries
library(GNAR)
library(QZ)
library(parallel)
library(igraph)

# ==============================================================================
# 2. GNAR Network Definition
# ==============================================================================

generate_gnar_network <- function(n_nodes, p = 0.2) {
  g <- igraph::sample_gnp(n = n_nodes, p = p, directed = FALSE)
  while (!igraph::is_connected(g)) {
    g <- igraph::sample_gnp(n = n_nodes, p = p, directed = FALSE)
  }
  adj_matrix <- as.matrix(igraph::as_adjacency_matrix(g))
  GNAR::as.GNARnet(adj_matrix)
}

set.seed(42)
num_simulations <- 500
n_nodes <- 10
n_models <- 5
max_alpha <- 4
max_beta  <- 5

gnar_network <- generate_gnar_network(n_nodes = n_nodes, p = 0.05)

# ==============================================================================
# 3. Load Auxiliary Functions and Data
# ==============================================================================

source("Aux_GNAR_spec.R")
# This is obtained in Simulation_10N_N.R
Datasets <- readRDS(file.path(dir.r, "Datasets_10.rds"))

# ==============================================================================
# 4. GNAR Spectrum Estimation under Misspecification
# ==============================================================================

Miss_est_gnar_spec_10 <- vector("list", length(Datasets))

for (ii in seq_along(Datasets)) {
  cat(sprintf("Dataset %d/%d\n", ii, length(Datasets)))
  data_list <- Datasets[[ii]]
  dataset_results <- vector("list", n_models)
  
  for (i in seq_len(n_models)) {
    cat(sprintf("  Model %d/%d in Dataset %d\n", i, n_models, ii))
    model_results <- vector("list", num_simulations)
    
    for (j in seq_len(num_simulations)) {
      cat(sprintf("    Simulation %d/%d | Model %d | Dataset %d\n", j, num_simulations, i, ii))
      model_results[[j]] <- GNARSpec_MM(
        data_list[[i]][[j]],
        net = gnar_network,
        max_alpha,
        max_beta,
        globalalpha = TRUE
      )
    }
    dataset_results[[i]] <- model_results
  }
  Miss_est_gnar_spec_10[[ii]] <- dataset_results
}

# ==============================================================================
# 5. RMSE Evaluation
# ==============================================================================

True_gnar_spec <- readRDS(file.path(dir.r, "True_gnar_spec.rds"))
True <- True_gnar_spec

n_datasets    <- 4
n_simulations <- 500
n_models      <- length(True[[1]])
n_nodes       <- dim(True[[1]][[1]])[1]
n_methods     <- 1
method_list   <- list(Miss_est_gnar_spec_10)

rmse_table      <- matrix(0, nrow = n_models, ncol = n_datasets)
coh_rmse_table  <- matrix(0, nrow = n_models, ncol = n_datasets)
pcoh_rmse_table <- matrix(0, nrow = n_models, ncol = n_datasets)

# ==============================================================================
# 6. Compute RMSEs for Spectrum, Coherence, and Partial Coherence
# ==============================================================================

for (dataset_idx in 1:n_datasets) {
  cat("Processing dataset", dataset_idx, "\n")
  
  for (model_idx in 1:n_models) {
    cat("  Model", model_idx, "\n")
    true_spectrum <- True[[dataset_idx]][[model_idx]]
    sel_freq <- dim(true_spectrum)[3]
    
    rmse_freq_sim  <- numeric(n_simulations)
    coh_freq_sim   <- numeric(n_simulations)
    pcoh_freq_sim  <- numeric(n_simulations)
    
    for (sim_idx in 1:n_simulations) {
      est_spectrum <- method_list[[1]][[dataset_idx]][[model_idx]][[sim_idx]]
      
      # --- Spectrum RMSE ---
      rmse_per_freq <- sapply(1:sel_freq, function(f)
        sqrt(mean(Mod(est_spectrum[,,f] - true_spectrum[,,f])^2))
      )
      rmse_freq_sim[sim_idx] <- mean(rmse_per_freq)
      
      # --- Coherence RMSE ---
      coh_per_freq <- numeric(sel_freq)
      for (f in 1:sel_freq) {
        coh_true <- matrix(1, n_nodes, n_nodes)
        coh_est  <- matrix(1, n_nodes, n_nodes)
        for (i in 1:n_nodes) {
          for (j in 1:n_nodes) {
            if (i != j) {
              Sii_true <- Mod(true_spectrum[i, i, f])
              Sjj_true <- Mod(true_spectrum[j, j, f])
              Sii_est  <- Mod(est_spectrum[i, i, f])
              Sjj_est  <- Mod(est_spectrum[j, j, f])
              coh_true[i, j] <- Mod(true_spectrum[i, j, f])^2 / (Sii_true * Sjj_true)
              coh_est[i, j]  <- Mod(est_spectrum[i, j, f])^2 / (Sii_est * Sjj_est)
            }
          }
        }
        coh_per_freq[f] <- sqrt(mean((coh_est - coh_true)^2))
      }
      coh_freq_sim[sim_idx] <- mean(coh_per_freq)
      
      # --- Partial Coherence RMSE ---
      pcoh_per_freq <- numeric(sel_freq)
      for (f in 1:sel_freq) {
        f_true <- matrix(true_spectrum[,,f], n_nodes, n_nodes)
        f_est  <- matrix(est_spectrum[,,f],  n_nodes, n_nodes)
        
        G_true <- tryCatch(solve(f_true), error = function(e) NULL)
        G_est  <- tryCatch(solve(f_est),  error = function(e) NULL)
        if (is.null(G_true) || is.null(G_est)) next
        
        D_true <- diag(1 / sqrt(abs(Re(diag(G_true)))))
        D_est  <- diag(1 / sqrt(abs(Re(diag(G_est)))))
        Gamma_true <- -D_true %*% G_true %*% D_true
        Gamma_est  <- -D_est  %*% G_est  %*% D_est
        pcoh_per_freq[f] <- sqrt(mean((Mod(Gamma_est)^2 - Mod(Gamma_true)^2)^2))
      }
      pcoh_freq_sim[sim_idx] <- mean(pcoh_per_freq, na.rm = TRUE)
    }
    
    rmse_table[model_idx, dataset_idx]      <- mean(rmse_freq_sim, na.rm = TRUE)
    coh_rmse_table[model_idx, dataset_idx]  <- mean(coh_freq_sim, na.rm = TRUE)
    pcoh_rmse_table[model_idx, dataset_idx] <- mean(pcoh_freq_sim, na.rm = TRUE)
  }
}

# ==============================================================================
# 7. Output
# ==============================================================================

resuts_rmse_MM_10 <- list(rmse_table, coh_rmse_table, pcoh_rmse_table)

print(round(resuts_rmse_MM_10[[1]] * 100, 2))  # Spectrum RMSE (%)
print(round(resuts_rmse_MM_10[[2]] * 100, 2))  # Coherence RMSE (%)
print(round(resuts_rmse_MM_10[[3]] * 100, 2))  # Partial Coherence RMSE (%)

################################################################################
# End of Simulation 3: Model Misspecification (10-Node Network)
################################################################################
