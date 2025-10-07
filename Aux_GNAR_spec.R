################################################################################
# Auxiliary Functions for Parametric and Nonparametric GNAR Spectrum Estimation
# Author: Cristian Felipe Jimenez Varon
# University of York, Department of Mathematics
################################################################################

# --- Libraries ----------------------------------------------------------------
library(Matrix)
library(glasso)
library(GNAR)
library(astsa)
library(future.apply)
library(QZ)
library(parallel)
library(igraph)
library(vars)
library(rlist)

################################################################################
# I. Helper Functions for GNAR Spectrum
################################################################################

get_alpha <- function(n.nodes, j, coefs, globalalpha) {
  if (globalalpha) {
    alphastr <- paste0("dmatalpha", j)
  } else {
    alphastr <- paste(paste0("dmatalpha", j, "node"), 1:n.nodes, sep = "")
  }
  diag(coefs[alphastr], n.nodes)
}

get_betas <- function(j, sj, Ws, coefs) {
  betastr <- paste(paste0("dmatbeta", j, "."), 1:sj, sep = "")
  Reduce("+", mapply(function(beta, W) beta * W, as.list(coefs[betastr]), Ws, SIMPLIFY = FALSE))
}

get_phi <- function(n.nodes, j, sj, Ws, globalalpha, coefs) {
  get_alpha(n.nodes, j, coefs, globalalpha) + get_betas(j, sj, Ws, coefs)
}

get_var_mixing_matrices <- function(vts, net, alphaOrder, betaOrder, globalalpha) {
  n.nodes <- length(net$edges)
  r_max <- max(betaOrder)
  W <- weights_matrix(net, r_max)
  adj_matrix <- as.matrix(GNARtoigraph(net))
  stages_adj <- get_k_stages_adjacency_tensor(adj_matrix, r_max)
  
  fit <- GNARfit(vts, net, alphaOrder, betaOrder, globalalpha, fact.var = NULL)
  coefs <- coef(fit)
  
  prods <- lapply(stages_adj, function(adj) W * adj)
  
  lapply(seq_len(alphaOrder), function(j) {
    Reduce(`+`, mapply(`*`, coefs[paste0("dmatbeta", j, ".", seq_len(betaOrder[j]))],
                       prods[1:betaOrder[j]], SIMPLIFY = FALSE)) +
      diag(coefs[paste0("dmatalpha", j)], n.nodes)
  })
}

################################################################################
# II. Parametric GNAR Spectrum
################################################################################

GNARSpec <- function(vts, net, alphaOrder, betaOrder, globalalpha = TRUE) {
  VAR_coef <- get_var_mixing_matrices(vts, net, alphaOrder, betaOrder, globalalpha)
  nf <- nrow(vts)
  ff <- c(0:(nf - 1)) / nf
  freq <- ff[ff < 0.5 & ff > 0]
  
  nc <- ncol(vts)
  Sigma <- diag(nc)
  GNAR_spec <- array(NA, dim = c(nc, nc, length(freq)))
  
  for (v in seq_along(freq)) {
    U <- diag(nc) - Reduce(`+`, Map(`*`, lapply(VAR_coef, as.matrix),
                                    lapply(1:alphaOrder, function(r)
                                      complex(real = cos(2 * pi * r * freq[v]),
                                              imaginary = -sin(2 * pi * r * freq[v])))))
    GNAR_spec[,,v] <- solve(U) %*% Sigma %*% QZ::H(solve(U))
  }
  
  GNAR_spec
}

# ---------------------
# GNAR Spectrum Estimation for Global Bank Network_conn
# ---------------------

get_var_mixing_matrices_GBNC <- function(vts, net, alphaOrder, betaOrder, globalalpha, W = NULL) {
  n.nodes <- length(net$edges)
  r_max <- max(betaOrder)
  
  if (is.null(W)) {
    W <- weights_matrix(net, r_max)
  }
  
  adj_matrix <- as.matrix(GNARtoigraph(net))
  stages_adj <- get_k_stages_adjacency_tensor(adj_matrix, r_max)
  fit <- GNARfit(vts, net, alphaOrder, betaOrder, globalalpha, fact.var = NULL)
  coefs <- coef(fit)
  
  Ws <- lapply(stages_adj, function(adj) W * adj)
  
  lapply(seq_len(alphaOrder), function(j) {
    beta_contrib <- Reduce(`+`, mapply(`*`,
                                       coefs[paste0("dmatbeta", j, ".", seq_len(betaOrder[j]))],
                                       Ws[1:betaOrder[j]],
                                       SIMPLIFY = FALSE))
    alpha_contrib <- diag(coefs[paste0("dmatalpha", j)], n.nodes)
    beta_contrib + alpha_contrib
  })
}
residual_cov_matrix <- function(vts, net, alphaOrder, betaOrder, globalalpha, W = NULL) {
  # Rename and format the input data
  data <- vts
  colnames(data) <- paste0("Var", 1:dim(vts)[2])
  T <- nrow(data)
  P <- ncol(data)
  
  # Get VAR(alphaOrder) coefficient matrices from GNAR model
  VAR_coeff <- get_var_mixing_matrices_GBNC(vts, net, alphaOrder, betaOrder, globalalpha, W = W)
  
  # Check that VAR_coeff is a list of matrices [Phi_1, ..., Phi_p]
  p <- alphaOrder
  stopifnot(length(VAR_coeff) == p)
  
  # Stack the VAR coefficient matrices horizontally: B = [Phi_1 | ... | Phi_p]
  coef_matrix <- list.cbind(VAR_coeff)  # matrix of size (P x p*P)
  
  # Initialize matrix for fitted values
  fitted_values <- matrix(0, nrow = T, ncol = P)
  for (t in (p + 1):T) {
    # Create lagged vector of length p * P: [X_{t-1}; X_{t-2}; ...; X_{t-p}]
    lagged_vec <- unlist(lapply(1:p, function(k) data[t - k, ]))
    
    if (length(lagged_vec) != P * p) {
      stop(paste0("lagged_vec should be of length ", P * p, 
                  " but got ", length(lagged_vec)))
    }
    
    # Multiply and assign safely
    fitted_values[t, ] <- as.vector(coef_matrix %*% lagged_vec)
  }
  
  
  # Calculate residuals: ε̂_t = X_t - X̂_t
  residuals <- data - fitted_values
  
  # Estimate residual covariance matrix using residuals from t = p+1 to T
  residual_cov_matrix <- cov(residuals[(p + 1):T, ])
  
  return(residual_cov_matrix)
}

# GNAR Spectrum for GBNC
GNARSpec_GBNC <- function(vts, net, alphaOrder, betaOrder, globalalpha, W = NULL) {
  VAR_coef <- get_var_mixing_matrices_GBNC(vts, net, alphaOrder, betaOrder, globalalpha, W)
  nf <- dim(vts)[1]
  ff <- c(0:(nf - 1)) / nf
  sel.f <- which(ff < 0.5 & ff > 0)
  freq <- ff[sel.f]
  
  nc <- dim(VAR_coef[[1]])[1]
  Sigma <- residual_cov_matrix(vts, net, alphaOrder, betaOrder, globalalpha, W = W)
  # sigma=1
  # Sigma <- diag(nc)*sigma
  GNAR_spec <- array(data = NA, dim = c(nc, nc, length(freq)))
  
  GNAR_spec <- array(NA, dim = c(nc, nc, length(freq)))
  
  for (v in seq_along(freq)) {
    U <- diag(nc) - Reduce(`+`, Map(`*`, lapply(VAR_coef, as.matrix), 
                                    lapply(1:alphaOrder, function(r) complex(real = cos(2 * pi * r * freq[v]), imaginary = -sin(2 * pi * r * freq[v])))))
    GNAR_spec[,,v] <- solve(U) %*% Sigma %*% QZ::H(solve(U))
  }
  GNAR_spec
}


################################################################################
# III. Nonparametric GNAR Spectrum (Network Agnostic)
################################################################################

nonparametric_gnar_spectrum <- function(Y, leng, regularize = FALSE, epsilon = 1e-6) {
  nf <- leng
  ff <- c(0:(nf - 1)) / nf
  freq <- ff[ff > 0 & ff < 0.5]
  
  n <- nrow(Y)
  p <- ncol(Y)
  nfreq <- length(freq)
  
  optimal_span <- max(5, round(sqrt(n)))
  optimal_taper <- min(0.1, max(0.01, 1 / sqrt(n)))
  
  spectrum_obj <- astsa::mvspec(Y, plot = FALSE,
                                spans = c(optimal_span, optimal_span),
                                taper = optimal_taper)
  xfft <- spectrum_obj$fxx[,,ff > 0 & ff < 0.5]
  
  if (regularize) {
    for (i in seq_len(nfreq)) {
      eig_vals <- eigen(xfft[,,i], only.values = TRUE)$values
      if (min(Re(eig_vals)) < epsilon) {
        xfft[,,i] <- xfft[,,i] + diag(epsilon, p)
      }
    }
  }
  
  list(
    spectrum = xfft,
    spans = c(optimal_span, optimal_span),
    taper = optimal_taper,
    frequencies = freq
  )
}

################################################################################
# IV. True GNAR Spectrum
################################################################################

compute_weights <- function(network, r_max) {
  lapply(1:r_max, function(x) weights_matrix(network, x))
}

compute_gnar_spectrum <- function(params, net, freq, sigma) {
  n_nodes <- length(net$edges)
  alphaOrder <- length(params$alpha)
  
  Beta <- sapply(params$beta, length)
  r_max <- max(Beta)
  weights <- weights_matrix(net, r_max)
  
  spectrum_array <- array(0, dim = c(n_nodes, n_nodes, length(freq)))
  
  for (v in seq_along(freq)) {
    omega <- freq[v]
    U <- diag(n_nodes)
    
    for (j in 1:alphaOrder) {
      alpha_j <- diag(params$alpha[[j]], n_nodes, n_nodes)
      beta_j <- matrix(0, n_nodes, n_nodes)
      
      for (r in seq_along(params$beta[[j]])) {
        s_r <- get_k_stages_adjacency_tensor(as.matrix(net), r)
        beta_j <- beta_j + params$beta[[j]][r] * weights * s_r[[r]]
      }
      U <- U - (alpha_j + beta_j) * exp(-1i * 2 * pi * j * omega)
    }
    
    U_inv <- solve(U)
    spectrum_array[,,v] <- U_inv %*% sigma %*% Conj(t(U_inv))
  }
  
  spectrum_array
}

################################################################################
# V. Regularization Utility
################################################################################

make_invertible_spectral_matrix <- function(Np_sdm, lambda = 1e-6) {
  n_nodes <- nrow(Np_sdm[,,1])
  n_freq <- dim(Np_sdm)[3]
  
  Np_sdm_invertible <- array(0, dim = c(n_nodes, n_nodes, n_freq))
  
  for (k in 1:n_freq) {
    spectral_matrix <- (Np_sdm[,,k] + t(Conj(Np_sdm[,,k]))) / 2
    if (any(eigen(spectral_matrix)$values <= 0)) {
      spectral_matrix <- spectral_matrix + lambda * diag(n_nodes)
    }
    Np_sdm_invertible[,,k] <- spectral_matrix
  }
  
  Np_sdm_invertible
}

################################################################################
# VI. VAR Spectrum (for Benchmark Comparison)
################################################################################

VAR_coeff <- function(vts, alphaOrder) {
  nc <- ncol(vts)
  VAR_fit <- suppressWarnings(vars::VAR(vts, p = alphaOrder, type = 'none'))
  VAR_coef <- coef(VAR_fit)
  
  VAR_coef_all <- matrix(NA, nrow = nc, ncol = alphaOrder * nc)
  for (j in 1:nc) {
    VAR_coef_all[j, ] <- VAR_coef[[j]][, 1]
  }
  
  VAR_mat_coef <- vector("list", length = alphaOrder)
  for (r in 1:alphaOrder) {
    start_col <- (r - 1) * nc + 1
    end_col <- r * nc
    VAR_mat_coef[[r]] <- VAR_coef_all[, start_col:end_col]
  }
  
  VAR_mat_coef
}

compute_VAR_spectrum_single <- function(vts, alphaOrder) {
  T_len <- nrow(vts)
  ff <- c(0:(T_len - 1)) / T_len
  freq <- ff[ff < 0.5 & ff > 0]
  
  nc <- ncol(vts)
  Sigma <- diag(nc)
  VAR_coef <- VAR_coeff(vts, alphaOrder)
  
  VAR_spec <- array(NA, dim = c(nc, nc, length(freq)))
  
  for (v in seq_along(freq)) {
    u <- matrix(0, nc, nc)
    for (r in 1:alphaOrder) {
      u <- u + VAR_coef[[r]] * complex(real = cos(2 * pi * r * freq[v]),
                                       imaginary = -sin(2 * pi * r * freq[v]))
    }
    U <- diag(nc) - u
    VAR_spec[,,v] <- solve(U) %*% Sigma %*% QZ::H(solve(U))
  }
  
  VAR_spec
}

################################################################################
# VII. Model Misspecification
################################################################################

select_gnar_model <- function(vts, net, max_alpha, max_beta) {
  min_BIC <- Inf
  best_model <- list(alphaOrder = NULL, betaOrder = NULL)
  
  for (alpha in 1:max_alpha) {
    beta_combinations <- expand.grid(replicate(alpha, 1:max_beta, simplify = FALSE))
    
    for (i in seq_len(nrow(beta_combinations))) {
      beta_vector <- as.numeric(beta_combinations[i, ])
      model <- GNARfit(vts, net, alpha, beta_vector)
      
      if (!is.null(model) && inherits(model, "GNARfit")) {
        bic_value <- BIC(model)
        if (bic_value < min_BIC) {
          min_BIC <- bic_value
          best_model <- list(alphaOrder = alpha, betaOrder = beta_vector)
        }
      }
    }
  }
  
  best_model
}

GNARSpec_MM <- function(vts, net, max_alpha, max_beta, globalalpha = TRUE) {
  model <- select_gnar_model(vts, net, max_alpha, max_beta)
  alphaOrder <- model$alphaOrder
  betaOrder <- model$betaOrder
  
  VAR_coef <- get_var_mixing_matrices(vts, net, alphaOrder, betaOrder, globalalpha)
  Sigma <- diag(length(net$edges))
  
  nf <- nrow(vts)
  freq <- (0:(nf - 1)) / nf
  freq <- freq[freq < 0.5 & freq > 0]
  
  nc <- ncol(vts)
  GNAR_spec <- array(NA, dim = c(nc, nc, length(freq)))
  
  for (v in seq_along(freq)) {
    U <- diag(nc) - Reduce(`+`, Map(`*`, lapply(VAR_coef, as.matrix),
                                    lapply(1:alphaOrder, function(r)
                                      complex(real = cos(2 * pi * r * freq[v]),
                                              imaginary = -sin(2 * pi * r * freq[v])))))
    GNAR_spec[,,v] <- solve(U) %*% Sigma %*% QZ::H(solve(U))
  }
  
  GNAR_spec
}
