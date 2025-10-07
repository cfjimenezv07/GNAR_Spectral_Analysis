################################################################################
# Title: Plotting Spectral Estimation Methods (Results from Simulation_10N_N.R)
# Author: Cristian Felipe Jimenez Varon
# Affiliation: University of York, Department of Mathematics
# Description:
#   This script generates plots comparing different spectrum estimation methods.
#   It uses results obtained from Simulation_10N_N.R and produces coherence and
#   partial coherence plots for multiple node pairs and GNAR model configurations.
################################################################################

# ============================ #
#          Directories         #
# ============================ #

dir.pc  <- "G:/My Drive/2024/York/"
dir.mac <- "~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/2024/York/"

dir <- dir.mac  # <- Set active directory here (change to dir.pc if on PC)
setwd(paste0(dir, "Rcodes"))

dir.r <- paste0(dir, "Results/")
dir.p <- paste0(dir, "Plots_paper_spec/")


# ============================ #
#        Helper Function       #
# ============================ #

save_plot <- function(filename, plot_obj) {
  pdf(file = file.path(dir.p, filename), width = 8, height = 6)
  print(plot_obj)
  dev.off()
}


# ============================ #
#          Libraries           #
# ============================ #

library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)


# ============================ #
#         Load Results         #
# ============================ #
# NOTE: All results below are obtained from Simulation_10N_N.R

estimated_gnar_spec_10 <- readRDS(paste0(dir.r, "estimated_gnar_spec_10.rds"))
all_VAR_spectra_10     <- readRDS(paste0(dir.r, "all_VAR_spectra_10.rds"))
Par_var_A_10           <- readRDS(paste0(dir.r, "Par_var_A_10.rds"))
Par_var_A1_10          <- readRDS(paste0(dir.r, "Par_var_A1_10.rds"))
Est_sdm_diff_l_A_10    <- readRDS(paste0(dir.r, "Est_sdm_diff_l_A_10.rds"))
Est_sdm_diff_l_A1_10   <- readRDS(paste0(dir.r, "Est_sdm_diff_l_A1_10.rds"))
Np_GNAR_AN_10_v1       <- readRDS(paste0(dir.r, "Np_GNAR_AN_10_v1.rds"))
True_gnar_spec         <- readRDS(paste0(dir.r, "True_gnar_spec.rds"))


# ============================ #
#          Methods List        #
# ============================ #

estimated_methods <- list(
  EM1 = estimated_gnar_spec_10,
  EM2 = all_VAR_spectra_10,
  EM3 = Par_var_A_10,
  EM4 = Par_var_A1_10,
  EM5 = Est_sdm_diff_l_A_10,
  EM6 = Est_sdm_diff_l_A1_10,
  EM7 = Np_GNAR_AN_10_v1
)


# ============================ #
#         Configuration        #
# ============================ #

n_nodes      <- 10
dataset_list <- c(4)
model_list   <- c(1, 2, 4)
sim_idx      <- 1
true         <- True_gnar_spec
results_list <- list()


# ============================ #
#        Coherence Setup       #
# ============================ #

for (dataset_id in dataset_list) {
  results_list[[paste0("Dataset_", dataset_id)]] <- list()
  
  for (method_name in names(estimated_methods)) {
    method_spectrum <- estimated_methods[[method_name]]
    
    for (data_index in dataset_list) {
      for (model_id in model_list) {
        true_spectrum <- true[[data_index]][[model_id]]
        sel_freq <- dim(true_spectrum)[3]
        est_spectrum <- method_spectrum[[data_index]][[model_id]][[sim_idx]]
        
        # Skip if missing
        if (is.null(est_spectrum)) next
        
        model_results <- list(coherence = list(), partial_coherence = list())
        
        # --- Coherence ---
        coh_true <- array(0, dim = c(n_nodes, n_nodes, sel_freq))
        coh_est  <- array(0, dim = c(n_nodes, n_nodes, sel_freq))
        
        for (freq_idx in 1:sel_freq) {
          for (node_a in 1:n_nodes) {
            for (node_b in 1:n_nodes) {
              if (node_a != node_b) {
                Sii_true <- Re(true_spectrum[node_a, node_a, freq_idx])
                Sjj_true <- Re(true_spectrum[node_b, node_b, freq_idx])
                Sii_est  <- Re(est_spectrum[node_a, node_a, freq_idx])
                Sjj_est  <- Re(est_spectrum[node_b, node_b, freq_idx])
                
                coh_true[node_a, node_b, freq_idx] <- Mod(true_spectrum[node_a, node_b, freq_idx])^2 / (Sii_true * Sjj_true)
                coh_est[node_a, node_b, freq_idx]  <- Mod(est_spectrum[node_a, node_b, freq_idx])^2 / (Sii_est * Sjj_est)
              }
            }
          }
        }
        
        # --- Partial Coherence ---
        pcoh_true <- array(0, dim = c(n_nodes, n_nodes, sel_freq))
        pcoh_est  <- array(0, dim = c(n_nodes, n_nodes, sel_freq))
        
        for (freq_idx in 1:sel_freq) {
          G_true <- solve(true_spectrum[,,freq_idx])
          G_est  <- solve(est_spectrum[,,freq_idx])
          
          D_true <- diag(1 / sqrt(abs(Re(diag(G_true)))))
          D_est  <- diag(1 / sqrt(abs(Re(diag(G_est)))))
          
          Gamma_true <- -D_true %*% G_true %*% D_true
          Gamma_est  <- -D_est %*% G_est %*% D_est
          
          pcoh_true[,,freq_idx] <- Mod(Gamma_true)^2
          pcoh_est[,,freq_idx]  <- Mod(Gamma_est)^2
        }
        
        model_results$coherence$true <- coh_true
        model_results$coherence$est  <- coh_est
        model_results$partial_coherence$true <- pcoh_true
        model_results$partial_coherence$est  <- pcoh_est
        
        results_list[[paste0("Dataset_", data_index)]][[paste0("Method_", method_name)]][[paste0("Model_", model_id)]] <- model_results
      }
    }
  }
}


# ============================ #
#       Plot Configuration     #
# ============================ #

freq_grid <- seq(0, 0.5, length.out = 499)

method_palette <- c(
  "True" = "red",
  "EM1"  = "#FFD700",
  "EM2"  = "#90EE90",
  "EM3"  = "#FFA07A",
  "EM4"  = "#ADD8E6",
  "EM5"  = "blue",
  "EM6"  = "#D3D3D3",
  "EM7"  = "#D8BFD8"
)

line_widths <- c(
  "True" = 1.5,
  "EM1"  = 1.2,
  "EM2"  = 1,
  "EM3"  = 1,
  "EM4"  = 1,
  "EM5"  = 0.8,
  "EM6"  = 0.8,
  "EM7"  = 0.8
)


# ============================ #
#       Plotting Function      #
# ============================ #

plot_spectrum_gg <- function(quantity, model, i, j, dataset = 4, save = TRUE, methods = NULL) {
  
  ylab <- ifelse(quantity == "coherence", "Coherence", "Partial Coherence")
  df_list <- list()
  dataset_key <- paste0("Dataset_", dataset)
  model_key <- paste0("Model_", model)
  full_order <- c("EM2", "EM4", "EM3", "EM6", "EM5", "EM7", "EM1", "True")
  method_order <- if (is.null(methods)) full_order else methods
  
  # True spectrum
  if ("True" %in% method_order) {
    true_path <- results_list[[dataset_key]]$Method_EM1[[model_key]][[quantity]]$true[i, j, ]
    df_list[[length(df_list) + 1]] <- data.frame(Frequency = freq_grid, Value = true_path, Method = "True")
  }
  
  # Estimated spectra
  for (method in method_order[method_order != "True"]) {
    method_key <- paste0("Method_", method)
    est_path <- results_list[[dataset_key]][[method_key]][[model_key]][[quantity]]$est[i, j, ]
    df_list[[length(df_list) + 1]] <- data.frame(Frequency = freq_grid, Value = est_path, Method = method)
  }
  
  # Combine
  plot_data <- do.call(rbind, df_list)
  plot_data <- plot_data %>% filter(!is.na(Value))
  plot_data$Method <- factor(plot_data$Method, levels = method_order)
  
  # Plot
  p <- ggplot(plot_data, aes(x = Frequency, y = Value, color = Method, linetype = Method, size = Method)) +
    geom_line() +
    geom_hline(yintercept = 0, color = "black", size = 0.3) +
    geom_vline(xintercept = 0, color = "black", size = 0.3) +
    scale_color_manual(values = method_palette) +
    scale_size_manual(values = line_widths) +
    labs(y = ylab, x = "Frequency") +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      axis.title = element_text(face = "bold")
    )
  
  if (save) {
    method_tag <- if (is.null(methods)) "All" else paste0(methods, collapse = "")
    filename <- paste0(quantity, "_model_", model, "_nodes_", i, "_", j, "_", method_tag, ".pdf")
    save_plot(filename, p)
  }
  
  print(p)
}


# ============================ #
#      Generate All Plots      #
# ============================ #

plot_cases <- list(
  list(model = 1, i = 1, j = 8),
  list(model = 2, i = 1, j = 2),
  list(model = 4, i = 1, j = 9),
  list(model = 1, i = 1, j = 10),
  list(model = 2, i = 1, j = 5),
  list(model = 4, i = 1, j = 3)
)

# --- All methods ---
for (case in plot_cases) {
  plot_spectrum_gg("coherence", model = case$model, i = case$i, j = case$j)
  plot_spectrum_gg("partial_coherence", model = case$model, i = case$i, j = case$j)
}

# --- Subset: EM5, EM6, EM7 ---
for (case in plot_cases[1:3]) {
  plot_spectrum_gg("coherence", model = case$model, i = case$i, j = case$j,
                   methods = c("EM5", "EM6", "EM7", "True"))
  plot_spectrum_gg("partial_coherence", model = case$model, i = case$i, j = case$j,
                   methods = c("EM5", "EM6", "EM7", "True"))
}

# --- Subset: EM1, EM2, EM3, EM4 ---
for (case in plot_cases[1:3]) {
  plot_spectrum_gg("coherence", model = case$model, i = case$i, j = case$j,
                   methods = c("EM1", "EM2", "EM3", "EM4", "True"))
  plot_spectrum_gg("partial_coherence", model = case$model, i = case$i, j = case$j,
                   methods = c("EM1", "EM2", "EM3", "EM4", "True"))
}

################################################################################
# End of Script
# Notes:
#   - Results and spectra were computed in Simulation_10N_N.R.
#   - Figures are saved in the directory defined by `dir.p`.
################################################################################
