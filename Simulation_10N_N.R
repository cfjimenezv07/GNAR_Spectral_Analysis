################################################################################
# Simulation Study: GNAR Spectrum Estimation (10 Nodes)
# Author: Cristian Felipe Jimenez Varon
# University of York, Department of Mathematics
# ------------------------------------------------------------------------------
# This script simulates GNAR processes, estimates parametric and nonparametric
# spectra, including network-penalized estimators, and computes RMSE metrics.
################################################################################

# --- Libraries ----------------------------------------------------------------
library(igraph)
library(GNAR)
library(glasso)  # For network-penalized spectrum

# --- Network generation -------------------------------------------------------
generate_gnar_network <- function(n_nodes, p = 0.2) {
  g <- igraph::sample_gnp(n = n_nodes, p = p, directed = FALSE)
  while (!igraph::is_connected(g)) {
    g <- igraph::sample_gnp(n = n_nodes, p = p, directed = FALSE)
  }
  GNAR::as.GNARnet(as.matrix(igraph::as_adjacency_matrix(g)))
}

set.seed(42)
gnar_network <- generate_gnar_network(n_nodes = 10, p = 0.05)
plot(gnar_network, main = "10-Node GNAR Network")

# --- Simulation parameters ----------------------------------------------------
num_simulations <- 500
n_nodes         <- 10
n_models        <- 5
data_l          <- c(100, 200, 500, 1000)

# --- True parameter settings --------------------------------------------------
params_list <- list(
  "GNAR(2,[1,1])"   = list(alpha = list(rep(0.2, n_nodes), rep(0.2, n_nodes)),
                           beta  = list(0.2, 0.1)),
  "GNAR(2,[1,2])"   = list(alpha = list(rep(0.1, n_nodes), rep(0.1, n_nodes)),
                           beta  = list(c(0.075), c(0.05, 0.15))),
  "GNAR(2,[2,3])"   = list(alpha = list(rep(0.2, n_nodes), rep(0.1, n_nodes)),
                           beta  = list(c(0.075, 0.05), c(0.05, 0.05, 0.1))),
  "GNAR(3,[1,2,3])" = list(alpha = list(rep(0.1, n_nodes), rep(0.075, n_nodes), rep(0.05, n_nodes)),
                           beta  = list(c(0.1), c(0.075, 0.075), c(0.05, 0.05, 0.05))),
  "GNAR(3,[3,3,3])" = list(alpha = list(rep(0.15, n_nodes), rep(0.1, n_nodes), rep(0.05, n_nodes)),
                           beta  = list(rep(0.05,3), rep(0.05,3), rep(0.05,3)))
)

# --- Generate simulated datasets ---------------------------------------------
Datasets_10 <- lapply(data_l, function(T_len) {
  lapply(params_list, function(param) {
    replicate(num_simulations, GNARsim(T_len, gnar_network, param$alpha, param$beta), simplify = FALSE)
  })
})

################################################################################
# Parametric GNAR Spectrum Estimation
################################################################################
source("Aux_GNAR_spec.R")  # Contains GNARSpec, compute_VAR_spectrum_single, compute_nonparam_gnar_spec

model_specs <- list(
  list(alphaOrder = 2, betaOrder = c(1, 1)),
  list(alphaOrder = 2, betaOrder = c(1, 2)),
  list(alphaOrder = 2, betaOrder = c(2, 3)),
  list(alphaOrder = 3, betaOrder = c(1, 2, 3)),
  list(alphaOrder = 3, betaOrder = c(3, 3, 3))
)

estimated_gnar_spec_10 <- lapply(Datasets_10, function(data_set) {
  lapply(seq_along(model_specs), function(model_idx) {
    spec <- model_specs[[model_idx]]
    lapply(seq_len(num_simulations), function(sim_idx) {
      vts <- data_set[[model_idx]][[sim_idx]]
      GNARSpec(vts, gnar_network, spec$alphaOrder, spec$betaOrder, globalalpha = TRUE)
    })
  })
})

################################################################################
# Parametric VAR Spectrum Estimation (network-sparse)
################################################################################
model_orders <- c(2, 2, 2, 3, 3)

all_VAR_spectra_10 <- lapply(seq_along(Datasets_10), function(dataset_idx) {
  data_set <- Datasets_10[[dataset_idx]]
  lapply(seq_len(n_models), function(model_idx) {
    alphaOrder <- model_orders[model_idx]
    lapply(seq_len(num_simulations), function(sim_idx) {
      vts <- data_set[[model_idx]][[sim_idx]]
      tryCatch(compute_VAR_spectrum_single(vts, alphaOrder), error = function(e) NULL)
    })
  })
})

################################################################################
# Nonparametric GNAR Spectrum (network-agnostic)
################################################################################
Np_GNAR_AN_10_v1 <- lapply(seq_along(Datasets_10), function(dataset_idx) {
  data_set <- Datasets_10[[dataset_idx]]
  nf <- data_l[dataset_idx]
  lapply(seq_len(n_models), function(model_idx) {
    lapply(data_set[[model_idx]], function(sim_data) {
      nonparametric_gnar_spectrum(sim_data, nf, regularize = FALSE)$spectrum
    })
  })
})

################################################################################
# Nonparametric/Parametric Estimation Penalized by the Network
################################################################################
beta_param <- list(c(1,1), c(1,2), c(2,3), c(1,2,3), c(3,3,3))

log_progress <- function(stage, ii = NA, model = NA, sim = NA, freq = NA) {
  msg <- paste0("🔹 ", stage)
  if (!is.na(ii)) msg <- paste0(msg, " | Dataset: ", ii)
  if (!is.na(model)) msg <- paste0(msg, " | Model: ", model)
  if (!is.na(sim)) msg <- paste0(msg, " | Sim: ", sim)
  if (!is.na(freq)) msg <- paste0(msg, " | Freq: ", freq)
  message(msg)
}

safe_execute <- function(expr, dataset_id, model_id) {
  tryCatch(expr, error = function(e) {
    message("❌ STOPPED at Dataset ", dataset_id, ", Model ", model_id)
    stop(e)
  })
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

run_network_penalized <- function(data_source, par_spectrum = FALSE, use_tensor = TRUE) {
  output_list <- vector("list", length(data_source))
  
  for (ii in seq_along(data_source)) {
    log_progress("Starting Dataset", ii)
    nf <- data_l[ii]
    ff <- seq(0, nf-1)/nf
    sel.f <- which(ff>0 & ff<0.5)
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
      
      Est_sdm[[i]] <- safe_execute(
        lapply(seq_len(num_simulations), function(j) {
          Np_sdm <- data_set[[i]][[j]]
          sdm_result <- array(0, dim = c(nrow(A)/2, nrow(A)/2, nfreq))
          
          for (k in seq_len(nfreq)) {
            C_k <- Re(Np_sdm[,,k])
            Q_k <- Im(Np_sdm[,,k])
            augmented_matrix <- 0.5*rbind(cbind(C_k, -Q_k), cbind(Q_k, C_k))
            
            if (r_max < diameter(GNARtoigraph(gnar_network)) && dim(zero_matrix)[1]!=0) {
              Aug_cov <- glasso::glasso(augmented_matrix, rho = 0, zero = zero_matrix)$w
              CXX <- Aug_cov[1:(nrow(A)/2),1:(nrow(A)/2)]
              CXY <- Aug_cov[1:(nrow(A)/2),(nrow(A)/2 + 1):nrow(A)]
              sdm_result[,,k] <- 2*CXX + 2i*CXY
            } else {
              sdm_result[,,k] <- Np_sdm[,,k]
            }
          }
          sdm_result
        }), ii, i
      )
    }
    output_list[[ii]] <- Est_sdm
    log_progress("Finished Dataset", ii)
  }
  output_list
}

Est_sdm_diff_l_A_10  <- run_network_penalized(Np_GNAR_AN_10_v1, par_spectrum = FALSE, use_tensor = TRUE)
Est_sdm_diff_l_A1_10 <- run_network_penalized(Np_GNAR_AN_10_v1, par_spectrum = FALSE, use_tensor = FALSE)
Par_var_A_10         <- run_network_penalized(all_VAR_spectra_10, par_spectrum = TRUE, use_tensor = TRUE)
Par_var_A1_10        <- run_network_penalized(all_VAR_spectra_10, par_spectrum = TRUE, use_tensor = FALSE)

################################################################################
# True GNAR Spectrum Computation
################################################################################
residual_cov_matrix <- diag(n_nodes)
True_gnar_spec <- lapply(seq_along(Datasets_10), function(ii) {
  nf <- data_l[ii]
  ff <- seq(0, nf-1)/nf
  freq <- ff[ff>0 & ff<0.5]
  
  lapply(params_list, function(params) {
    compute_gnar_spectrum(params, gnar_network, freq, sigma = residual_cov_matrix)
  })
})

################################################################################
# RMSE Computation for all 7 Methods
################################################################################
method_list <- list(
  estimated_gnar_spec_10,
  all_VAR_spectra_10,
  Par_var_A_10,
  Par_var_A1_10,
  Est_sdm_diff_l_A_10,
  Est_sdm_diff_l_A1_10,
  Np_GNAR_AN_10_v1
)

n_methods <- length(method_list)
n_datasets <- length(Datasets_10)
n_models <- length(params_list)

rmse_table <- matrix(0, nrow = n_models, ncol = n_datasets * n_methods)

for (dataset_idx in seq_len(n_datasets)) {
  for (model_idx in seq_len(n_models)) {
    true_spectrum <- True_gnar_spec[[dataset_idx]][[model_idx]]
    sel_freq <- dim(true_spectrum)[3]
    
    for (method_idx in seq_len(n_methods)) {
      method_data <- method_list[[method_idx]]
      
      rmse_freq_sim <- sapply(seq_len(num_simulations), function(sim_idx) {
        est_spectrum <- method_data[[dataset_idx]][[model_idx]][[sim_idx]]
        rmse_per_freq <- sapply(seq_len(sel_freq), function(f) {
          sqrt(mean(Mod(est_spectrum[,,f] - true_spectrum[,,f])^2))
        })
        mean(rmse_per_freq)
      })
      rmse_table[model_idx, (dataset_idx-1)*n_methods + method_idx] <- mean(rmse_freq_sim, na.rm = TRUE)
    }
  }
}

print(round(rmse_table*100,2))
