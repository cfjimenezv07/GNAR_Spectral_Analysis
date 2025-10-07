################################################################################
# Simulation Study: GNAR Spectrum Estimation (10 Nodes)
# Author: Cristian Felipe Jimenez Varon
# University of York, Department of Mathematics
# ------------------------------------------------------------------------------
# This script simulates GNAR processes, estimates parametric and nonparametric
# spectra, and computes RMSE metrics for various estimators.
################################################################################

# --- Setup --------------------------------------------------------------------

library(igraph)
library(GNAR)

# --- Network generation -------------------------------------------------------

generate_gnar_network <- function(n_nodes, p = 0.2) {
  g <- igraph::sample_gnp(n = n_nodes, p = p, directed = FALSE)
  while (!igraph::is_connected(g)) {
    g <- igraph::sample_gnp(n = n_nodes, p = p, directed = FALSE)
  }
  adj_matrix <- as.matrix(igraph::as_adjacency_matrix(g))
  GNAR::as.GNARnet(adj_matrix)
}

set.seed(42)
gnar_network <- generate_gnar_network(n_nodes = 10, p = 0.05)
plot(gnar_network)

# --- Simulation parameters ----------------------------------------------------

num_simulations <- 500
n_nodes <- 10
n_models <- 5

# --- True parameter settings --------------------------------------------------

params_list <- list(
  "GNAR(2,[1,1])" = list(
    alpha = list(rep(0.2, n_nodes), rep(0.2, n_nodes)),
    beta  = list(0.2, 0.1)
  ),
  "GNAR(2,[1,2])" = list(
    alpha = list(rep(0.1, n_nodes), rep(0.1, n_nodes)),
    beta  = list(c(0.075), c(0.05, 0.15))
  ),
  "GNAR(2,[2,3])" = list(
    alpha = list(rep(0.2, n_nodes), rep(0.1, n_nodes)),
    beta  = list(c(0.075, 0.05), c(0.05, 0.05, 0.1))
  ),
  "GNAR(3,[1,2,3])" = list(
    alpha = list(rep(0.1, n_nodes), rep(0.075, n_nodes), rep(0.05, n_nodes)),
    beta  = list(c(0.1), c(0.075, 0.075), c(0.05, 0.05, 0.05))
  ),
  "GNAR(3,[3,3,3])" = list(
    alpha = list(rep(0.15, n_nodes), rep(0.1, n_nodes), rep(0.05, n_nodes)),
    beta  = list(rep(0.05, 3), rep(0.05, 3), rep(0.05, 3))
  )
)

# --- Generate simulated datasets ---------------------------------------------

results_100  <- list()
results_200  <- list()
results_500  <- list()
results_1000 <- list()

for (model_name in names(params_list)) {
  alpha_params <- params_list[[model_name]]$alpha
  beta_params  <- params_list[[model_name]]$beta
  
  results_100[[model_name]]  <- replicate(num_simulations, GNARsim(100,  gnar_network, alpha_params, beta_params), simplify = FALSE)
  results_200[[model_name]]  <- replicate(num_simulations, GNARsim(200,  gnar_network, alpha_params, beta_params), simplify = FALSE)
  results_500[[model_name]]  <- replicate(num_simulations, GNARsim(500,  gnar_network, alpha_params, beta_params), simplify = FALSE)
  results_1000[[model_name]] <- replicate(num_simulations, GNARsim(1000, gnar_network, alpha_params, beta_params), simplify = FALSE)
}

Datasets_10 <- list(results_100, results_200, results_500, results_1000)
data_l <- c(100, 200, 500, 1000)

################################################################################
# Parametric GNAR Spectrum Estimation
################################################################################

source("Aux_GNAR_spec.R")

model_specs <- list(
  list(alphaOrder = 2, betaOrder = c(1, 1)),
  list(alphaOrder = 2, betaOrder = c(1, 2)),
  list(alphaOrder = 2, betaOrder = c(2, 3)),
  list(alphaOrder = 3, betaOrder = c(1, 2, 3)),
  list(alphaOrder = 3, betaOrder = c(3, 3, 3))
)

estimated_gnar_spec_10 <- vector("list", length(Datasets_10))

for (dataset_idx in seq_along(Datasets_10)) {
  data_list <- Datasets_10[[dataset_idx]]
  
  estimated_gnar_spec_10[[dataset_idx]] <- lapply(seq_along(model_specs), function(model_idx) {
    spec <- model_specs[[model_idx]]
    lapply(seq_len(num_simulations), function(sim_idx) {
      vts <- data_list[[model_idx]][[sim_idx]]
      GNARSpec(vts, gnar_network, spec$alphaOrder, spec$betaOrder, globalalpha = TRUE)
    })
  })
}

################################################################################
# Parametric VAR Spectrum Estimation
################################################################################

model_orders <- c(2, 2, 2, 3, 3)

all_VAR_spectra_10 <- lapply(seq_along(data_l), function(i) {
  lapply(seq_len(n_models), function(j) {
    alphaOrder <- model_orders[j]
    lapply(seq_len(num_simulations), function(k) {
      vts <- Datasets_10[[i]][[j]][[k]]
      tryCatch(compute_VAR_spectrum_single(vts, alphaOrder), error = function(e) NULL)
    })
  })
})

################################################################################
# Nonparametric Spectrum Estimation (Network Agnostic)
################################################################################

Np_GNAR_AN_10_v1 <- lapply(seq_along(Datasets_10), function(ii) {
  results <- Datasets_10[[ii]]
  leng <- data_l[ii]
  lapply(seq_along(results), function(model_idx) {
    lapply(results[[model_idx]], function(simulation_data) {
      nonparametric_gnar_spectrum(simulation_data, leng, regularize = FALSE)$spectrum
    })
  })
})

################################################################################
# True GNAR Spectrum Computation
################################################################################

residual_cov_matrix <- diag(n_nodes)
True_gnar_spec <- vector("list", length(Datasets_10))

for (ii in seq_along(Datasets_10)) {
  nf <- data_l[ii]
  ff <- seq(0, nf - 1) / nf
  freq <- ff[ff > 0 & ff < 0.5]
  
  True_gnar_spec[[ii]] <- lapply(params_list, function(params) {
    compute_gnar_spectrum(params, gnar_network, freq, sigma = residual_cov_matrix)
  })
}

################################################################################
# RMSE Computation for Spectral Estimators
################################################################################

method_list <- list(
  estimated_gnar_spec_10,
  all_VAR_spectra_10,
  Np_GNAR_AN_10_v1
)

n_datasets <- length(Datasets_10)
n_methods  <- length(method_list)
n_models   <- length(params_list)
n_simulations <- num_simulations
n_nodes <- dim(True_gnar_spec[[1]][[1]])[1]

rmse_table <- matrix(0, nrow = n_models, ncol = n_datasets * n_methods)

for (dataset_idx in seq_len(n_datasets)) {
  for (model_idx in seq_len(n_models)) {
    true_spectrum <- True_gnar_spec[[dataset_idx]][[model_idx]]
    sel_freq <- dim(true_spectrum)[3]
    
    for (method_idx in seq_len(n_methods)) {
      method_data <- method_list[[method_idx]]
      
      rmse_freq_sim <- sapply(seq_len(n_simulations), function(sim_idx) {
        est_spectrum <- method_data[[dataset_idx]][[model_idx]][[sim_idx]]
        rmse_per_freq <- sapply(seq_len(sel_freq), function(f) {
          sqrt(mean(Mod(est_spectrum[,,f] - true_spectrum[,,f])^2))
        })
        mean(rmse_per_freq)
      })
      rmse_table[model_idx, (dataset_idx - 1) * n_methods + method_idx] <- mean(rmse_freq_sim, na.rm = TRUE)
    }
  }
}

print(round(rmse_table * 100, 2))
