# =====================================================================
# Simulation and GNAR Spectrum Analysis for 5 Nodes
# =====================================================================

# ------------------------------
# Directory setup
# ------------------------------
dir.pc <- "G:/My Drive/2024/York/"
dir.mac <- "~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/2024/York/"
dir <- dir.mac  # Select active path

setwd(paste0(dir, "Rcodes"))
dir.r <- paste0(dir, "Results/")

# ------------------------------
# Load required libraries
# ------------------------------
library(igraph)
library(GNAR)
library(glasso)

# ------------------------------
# Network and simulation setup
# ------------------------------
set.seed(42)
gnar_network <- GNAR::fiveNet
plot(gnar_network)

num_simulations <- 500
n_nodes <- 5
n_models <- 5
data_l <- c(100, 200, 500, 1000)

# ------------------------------
# True GNAR parameters (stationary)
# ------------------------------
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

# ------------------------------
# Generate simulated GNAR datasets
# ------------------------------
results_list <- vector("list", length(data_l))
names(results_list) <- paste0("T", data_l)

for (t_idx in seq_along(data_l)) {
  n_time <- data_l[t_idx]
  cat("Simulating datasets for T =", n_time, "\n")
  results_t <- vector("list", n_models)
  names(results_t) <- names(params_list)
  
  for (model_name in names(params_list)) {
    alpha_params <- params_list[[model_name]]$alpha
    beta_params <- params_list[[model_name]]$beta
    results_t[[model_name]] <- replicate(
      num_simulations,
      GNARsim(n = n_time, net = gnar_network, alphaParams = alpha_params, betaParams = beta_params),
      simplify = FALSE
    )
  }
  results_list[[t_idx]] <- results_t
}

Datasets_5_v1 <- results_list

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
# Parametric VAR spectrum estimation
# =====================================================================
model_orders <- c(2, 2, 2, 3, 3)

all_VAR_spectra_5 <- lapply(1:length(Datasets_5_v1), function(i) {
  lapply(1:n_models, function(j) {
    alphaOrder <- model_orders[j]
    lapply(1:num_simulations, function(k) {
      vts <- Datasets_5_v1[[i]][[j]][[k]]
      tryCatch({
        compute_VAR_spectrum_single(vts, alphaOrder)
      }, error = function(e) {
        message(sprintf("Error at dataset=%d, model=%d, sim=%d: %s", i, j, k, e$message))
        NULL
      })
    })
  })
})

# =====================================================================
# Nonparametric (agnostic) GNAR spectrum
# =====================================================================
Np_GNAR_AN_5_v1 <- lapply(seq_along(Datasets_5_v1), function(ii) {
  results <- Datasets_5_v1[[ii]]
  leng <- data_l[ii]
  lapply(seq_along(results), function(model_idx) {
    simulations <- results[[model_idx]]
    lapply(simulations, function(sim_data) {
      nonparametric_gnar_spectrum(sim_data, leng, regularize = FALSE)$spectrum
    })
  })
})

# =====================================================================
# Nonparametric GNAR spectrum penalized by adjacency
# =====================================================================
beta_param <- list(c(1,1), c(1,2), c(2,3), c(1,2,3), c(3,3,3))

log_progress <- function(stage, ii = NA, model = NA, sim = NA, freq = NA) {
  msg <- paste0("🔹 ", stage)
  if (!is.na(ii)) msg <- paste0(msg, " | Dataset: ", ii)
  if (!is.na(model)) msg <- paste0(msg, " | Model: ", model)
  if (!is.na(sim)) msg <- paste0(msg, " | Sim: ", sim)
  if (!is.na(freq)) msg <- paste0(msg, " | Freq: ", freq)
  message(msg)
}

construct_A_matrices <- function(r_max, use_tensor = TRUE) {
  if (use_tensor) {
    tensor <- get_k_stages_adjacency_tensor(as.matrix(gnar_network), 2 * r_max)
    S <- construct_S(tensor, 2 * r_max)
    A <- construct_A(S)
  } else {
    A <- construct_A(get.adjacency(GNARtoigraph(gnar_network)))
  }
  list(A = A, zero_matrix = as.matrix(find_zeros(A)))
}

run_simulation <- function(data_source, par_spectrum = FALSE, use_tensor = TRUE) {
  output_list <- vector("list", length(data_source))
  
  for (ii in seq_along(data_source)) {
    log_progress("Starting Dataset", ii)
    nf <- data_l[ii]
    ff <- seq(0, nf - 1) / nf
    sel.f <- which(ff > 0 & ff < 0.5)
    freq <- ff[sel.f]
    nfreq <- length(freq)
    
    data_set <- data_source[[ii]]
    Est_sdm <- vector("list", n_models)
    
    for (i in seq_len(n_models)) {
      log_progress("Processing Model", ii, i)
      r_max <- max(beta_param[[i]])
      A_info <- construct_A_matrices(r_max, use_tensor)
      A <- A_info$A
      zero_matrix <- A_info$zero_matrix
      
      Est_sdm[[i]] <- lapply(seq_len(num_simulations), function(j) {
        Np_sdm <- data_set[[i]][[j]]
        sdm_result <- array(0, dim = c(nrow(A)/2, nrow(A)/2, nfreq))
        
        for (k in seq_len(nfreq)) {
          C_k <- Re(Np_sdm[,,k])
          Q_k <- Im(Np_sdm[,,k])
          augmented_matrix <- 0.5 * rbind(cbind(C_k, -Q_k), cbind(Q_k, C_k))
          
          if (r_max < diameter(GNARtoigraph(gnar_network)) && dim(zero_matrix)[1] != 0) {
            Aug_cov <- glasso::glasso(augmented_matrix, rho = 0, zero = zero_matrix)$w
            CXX <- Aug_cov[1:(nrow(A)/2), 1:(nrow(A)/2)]
            CXY <- Aug_cov[1:(nrow(A)/2), (nrow(A)/2 + 1):nrow(A)]
            sdm_result[,,k] <- 2*CXX + 2i*CXY
          } else {
            sdm_result[,,k] <- Np_sdm[,,k]
          }
        }
        sdm_result
      })
    }
    output_list[[ii]] <- Est_sdm
    log_progress("Finished Dataset", ii)
  }
  return(output_list)
}

# Run penalized nonparametric simulations
Est_sdm_diff_l_A_5  <- run_simulation(Np_GNAR_AN_5_v1, par_spectrum = FALSE, use_tensor = TRUE)
Est_sdm_diff_l_A1_5 <- run_simulation(Np_GNAR_AN_5_v1, par_spectrum = FALSE, use_tensor = FALSE)


# Run penalized parametric simulations
Par_var_A_5  <- run_simulation(all_VAR_spectra_5, par_spectrum = FALSE, use_tensor = TRUE)
Par_var_A1_5 <- run_simulation(all_VAR_spectra_5, par_spectrum = FALSE, use_tensor = FALSE)

# =====================================================================
# True GNAR spectrum
# =====================================================================
residual_cov_matrix <- diag(n_nodes)
True_gnar_spec_5 <- lapply(seq_along(Datasets_5_v1), function(ii) {
  nf <- data_l[ii]
  ff <- seq(0, nf - 1) / nf
  sel.f <- which(ff > 0 & ff < 0.5)
  freq <- ff[sel.f]
  lapply(params_list, function(params) {
    compute_gnar_spectrum(params, gnar_network, freq, sigma = residual_cov_matrix)
  })
})

# =====================================================================
# RMSE, Coherence, Partial Coherence
# =====================================================================
n_datasets <- length(Datasets_5_v1)
n_simulations <- num_simulations
n_models <- length(params_list)
n_nodes <- dim(True_gnar_spec_5[[1]][[1]])[1]

method_list <- list(
  estimated_gnar_spec_5_v1,
  all_VAR_spectra_5,
  Par_var_A_5,
  Par_var_A1_5,
  Est_sdm_diff_l_A_5,
  Est_sdm_diff_l_A1_5,
  Np_GNAR_AN_5_v1
)
n_methods <- length(method_list)

rmse_table      <- matrix(0, nrow = n_models, ncol = n_datasets * n_methods)
coh_rmse_table  <- rmse_table
pcoh_rmse_table <- rmse_table

for (dataset_idx in 1:n_datasets) {
  cat("Processing dataset", dataset_idx, "\n")
  for (model_idx in 1:n_models) {
    cat("  Model", model_idx, "\n")
    true_spectrum <- True_gnar_spec_5[[dataset_idx]][[model_idx]]
    sel_freq <- dim(true_spectrum)[3]
    
    for (method_idx in 1:n_methods) {
      cat("    Method", method_idx, "\n")
      rmse_freq_sim <- coh_freq_sim <- pcoh_freq_sim <- numeric(n_simulations)
      method_data <- method_list[[method_idx]]
      
      for (sim_idx in 1:n_simulations) {
        est_spectrum <- method_data[[dataset_idx]][[model_idx]][[sim_idx]]
        
        # RMSE
        rmse_per_freq <- sapply(1:sel_freq, function(f) sqrt(mean(Mod(est_spectrum[,,f] - true_spectrum[,,f])^2)))
        rmse_freq_sim[sim_idx] <- mean(rmse_per_freq)
        
        # Coherence RMSE
        coh_per_freq <- numeric(sel_freq)
        for (f in 1:sel_freq) {
          coh_true <- coh_est <- matrix(1, n_nodes, n_nodes)
          for (i in 1:n_nodes) for (j in 1:n_nodes) if (i != j) {
            Sii_true <- Re(true_spectrum[i,i,f])
            Sjj_true <- Re(true_spectrum[j,j,f])
            Sii_est  <- Re(est_spectrum[i,i,f])
            Sjj_est  <- Re(est_spectrum[j,j,f])
            coh_true[i,j] <- Mod(true_spectrum[i,j,f])^2 / (Sii_true*Sjj_true)
            coh_est[i,j]  <- Mod(est_spectrum[i,j,f])^2 / (Sii_est*Sjj_est)
          }
          coh_per_freq[f] <- sqrt(mean((coh_est - coh_true)^2))
        }
        coh_freq_sim[sim_idx] <- mean(coh_per_freq)
        
        # Partial coherence RMSE
        pcoh_per_freq <- numeric(sel_freq)
        for (f in 1:sel_freq) {
          G_true <- tryCatch(solve(true_spectrum[,,f]), error = function(e) NULL)
          G_est  <- tryCatch(solve(est_spectrum[,,f]), error = function(e) NULL)
          if (is.null(G_true) || is.null(G_est)) {
            pcoh_per_freq[f] <- NA
            next
          }
          D_true <- diag(1/sqrt(abs(Re(diag(G_true)))))
          D_est  <- diag(1/sqrt(abs(Re(diag(G_est)))))
          Gamma_true <- -D_true %*% G_true %*% D_true
          Gamma_est  <- -D_est  %*% G_est  %*% D_est
          pcoh_per_freq[f] <- sqrt(mean((Mod(Gamma_est)^2 - Mod(Gamma_true)^2)^2))
        }
        pcoh_freq_sim[sim_idx] <- mean(pcoh_per_freq, na.rm = TRUE)
      }
      
      # Store average over simulations
      col_idx <- (dataset_idx - 1) * n_methods + method_idx
      rmse_table[model_idx, col_idx] <- mean(rmse_freq_sim, na.rm = TRUE)
      coh_rmse_table[model_idx, col_idx] <- mean(coh_freq_sim, na.rm = TRUE)
      pcoh_rmse_table[model_idx, col_idx] <- mean(pcoh_freq_sim, na.rm = TRUE)
    }
  }
}

resuts_rmse_5 <- list(rmse_table, coh_rmse_table, pcoh_rmse_table)

# Final Output (percentage)
print(round(resuts_rmse_5[[1]]*100,2))
print(round(resuts_rmse_5[[2]]*100,2))
print(round(resuts_rmse_5[[3]]*100,2))
