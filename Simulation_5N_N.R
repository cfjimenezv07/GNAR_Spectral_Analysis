# =====================================================================
# Simulation with more nodes
# =====================================================================

# Directory setup (choose one)
dir.pc <- "G:/My Drive/2024/York/"
dir.mac <- "~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/2024/York/"
dir <- dir.mac  # select active path

setwd(paste0(dir, "Rcodes"))
dir.r <- paste0(dir, "Results/")

# Libraries
library(igraph)
library(GNAR)

# =====================================================================
# Generate network and define parameters
# =====================================================================

set.seed(42)
gnar_network <- GNAR::fiveNet
plot(gnar_network)

num_simulations <- 500
n_nodes <- 5
n_models <- 5

# True model parameters (stationary)
params_list <- list(
  "GNAR(2,[1,1])" = list(
    alpha = list(rep(0.2, n_nodes), rep(0.2, n_nodes)),
    beta  = list(c(0.2), c(0.1))
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
    beta  = list(c(0.05, 0.05, 0.05), c(0.05, 0.05, 0.05), c(0.05, 0.05, 0.05))
  )
)

# =====================================================================
# Generate simulations
# =====================================================================

results_100  <- list()
results_200  <- list()
results_500  <- list()
results_1000 <- list()

for (model_name in names(params_list)) {
  alpha_params <- params_list[[model_name]]$alpha
  beta_params  <- params_list[[model_name]]$beta
  
  simulations_100 <- replicate(num_simulations,
                               GNARsim(n = 100, net = gnar_network, alphaParams = alpha_params, betaParams = beta_params),
                               simplify = FALSE
  )
  simulations_200 <- replicate(num_simulations,
                               GNARsim(n = 200, net = gnar_network, alphaParams = alpha_params, betaParams = beta_params),
                               simplify = FALSE
  )
  simulations_500 <- replicate(num_simulations,
                               GNARsim(n = 500, net = gnar_network, alphaParams = alpha_params, betaParams = beta_params),
                               simplify = FALSE
  )
  simulations_1000 <- replicate(num_simulations,
                                GNARsim(n = 1000, net = gnar_network, alphaParams = alpha_params, betaParams = beta_params),
                                simplify = FALSE
  )
  
  results_100[[model_name]]  <- simulations_100
  results_200[[model_name]]  <- simulations_200
  results_500[[model_name]]  <- simulations_500
  results_1000[[model_name]] <- simulations_1000
}

Datasets_5_v1 <- list(results_100, results_200, results_500, results_1000)
data_l <- c(100, 200, 500, 1000)

# =====================================================================
# Parametric GNAR spectrum estimation
# =====================================================================

source("Aux_GNAR_spec.R")

model_specs <- list(
  list(alphaOrder = 2, betaOrder = c(1, 1)),
  list(alphaOrder = 2, betaOrder = c(1, 2)),
  list(alphaOrder = 2, betaOrder = c(2, 3)),
  list(alphaOrder = 3, betaOrder = c(1, 2, 3)),
  list(alphaOrder = 3, betaOrder = c(3, 3, 3))
)

estimated_gnar_spec_5_v1 <- vector("list", length(Datasets_5_v1))

for (dataset_idx in seq_along(Datasets_5_v1)) {
  data_list <- Datasets_5_v1[[dataset_idx]]
  estimated_gnar_spec_5_v1[[dataset_idx]] <- lapply(seq_along(model_specs), function(model_idx) {
    spec <- model_specs[[model_idx]]
    lapply(seq_len(num_simulations), function(sim_idx) {
      vts <- data_list[[model_idx]][[sim_idx]]
      GNARSpec(vts, gnar_network, spec$alphaOrder, spec$betaOrder, globalalpha = TRUE)
    })
  })
}

# =====================================================================
# Parametric VAR spectrum
# =====================================================================

model_orders <- c(2, 2, 2, 3, 3)

all_VAR_spectra_5 <- lapply(1:4, function(i) {
  lapply(1:5, function(j) {
    alphaOrder <- model_orders[j]
    lapply(1:500, function(k) {
      vts <- Datasets_5_v1[[i]][[j]][[k]]
      tryCatch({
        compute_VAR_spectrum_single(vts, alphaOrder)
      }, error = function(e) {
        message(sprintf("Error at i=%d j=%d k=%d: %s", i, j, k, e$message))
        NULL
      })
    })
  })
})

# =====================================================================
# Nonparametric (agnostic) spectrum
# =====================================================================

Np_GNAR_AN_5_v1 <- list()
for (ii in seq_along(Datasets_5_v1)) {
  results <- Datasets_5_v1[[ii]]
  leng <- data_l[[ii]]
  Np_GNAR_AN_5_v1[[ii]] <- lapply(seq_along(results), function(model_idx) {
    simulations <- results[[model_idx]]
    lapply(simulations, function(simulation_data) {
      nonparametric_gnar_spectrum(simulation_data, leng, regularize = FALSE)$spectrum
    })
  })
}

# =====================================================================
# True GNAR spectrum
# =====================================================================

true_params_list <- params_list
residual_cov_matrix <- diag(n_nodes)
net <- gnar_network
True_gnar_spec_5 <- list()

for (ii in 1:length(Datasets_5_v1)) {
  nf <- data_l[ii]
  ff <- c(0:(nf - 1)) / nf
  sel.f <- which(ff > 0 & ff < 0.5)
  freq <- ff[sel.f]
  True_gnar_spec_5[[ii]] <- lapply(true_params_list, function(params) {
    compute_gnar_spectrum(params, net, freq, sigma = residual_cov_matrix)
  })
}

# =====================================================================
# RMSE, Coherence, Partial Coherence evaluation
# =====================================================================

True <- True_gnar_spec_5
n_datasets <- 4
n_simulations <- 500
n_models <- length(True[[1]])
n_nodes <- dim(True[[1]][[1]])[1]

method_list <- list(
  estimated_gnar_spec_5_v1,
  all_VAR_spectra_5,
  Np_GNAR_AN_5_v1
)
n_methods <- length(method_list)

rmse_table <- matrix(0, nrow = n_models, ncol = n_datasets * n_methods)
coh_rmse_table <- rmse_table
pcoh_rmse_table <- rmse_table

for (dataset_idx in 1:n_datasets) {
  cat("Processing dataset", dataset_idx, "\n")
  for (model_idx in 1:n_models) {
    cat("  Model", model_idx, "\n")
    true_spectrum <- True[[dataset_idx]][[model_idx]]
    sel_freq <- dim(true_spectrum)[3]
    for (method_idx in 1:n_methods) {
      cat("    Method", method_idx, "\n")
      rmse_freq_sim <- coh_freq_sim <- pcoh_freq_sim <- numeric(n_simulations)
      method_data <- method_list[[method_idx]]
      for (sim_idx in 1:n_simulations) {
        est_spectrum <- method_data[[dataset_idx]][[model_idx]][[sim_idx]]
        rmse_per_freq <- sapply(1:sel_freq, function(f)
          sqrt(mean(Mod(est_spectrum[,,f] - true_spectrum[,,f])^2)))
        rmse_freq_sim[sim_idx] <- mean(rmse_per_freq)
      }
      rmse_table[model_idx, (dataset_idx - 1) * n_methods + method_idx] <- mean(rmse_freq_sim, na.rm = TRUE)
    }
  }
}

resuts_rmse_5 <- list(rmse_table, coh_rmse_table, pcoh_rmse_table)
print(round(resuts_rmse_5[[1]] * 100, 2))
print(round(resuts_rmse_5[[2]] * 100, 2))
print(round(resuts_rmse_5[[3]] * 100, 2))
