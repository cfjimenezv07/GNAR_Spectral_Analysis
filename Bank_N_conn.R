###############################################################################
# GLOBAL BANK NETWORK CONNECTEDNESS AND SPECTRAL FEATURES
###############################################################################

# This dataset was obtained from the supplementary material of:
# Demirer et al. (2018, Journal of Applied Econometrics)
library(data.table)
library(rugarch)
library(locits)
library(glmnet)
library(GNAR)
library(igraph)
library(RColorBrewer)
library(scales)
library(reshape2)
library(ggplot2)
library(stringr)
install.packages(file.path("packages", "fastCorbit_0.2.0.tar.gz"), 
                 repo = NULL, type = "source")
library(fastCorbit)

dir.p <- "~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/Financial_app/Plots/"
dir.r <- "~/Library/CloudStorage/GoogleDrive-cristian.jimenezvaron@york.ac.uk/My Drive/Financial_app/Rcodes/"
# =========================
# 0. SETUP BANK CODES
# =========================
bank_codes <- c(
  "hsba.gb", "mtbh.jp", "bnp.fr", "jpm.us", "dbk.de", "barc.gb", "aca.fr", "bac.us", "c.us",
  "mzh.jp", "gle.fr", "rbs.gb", "smtm.jp", "san.es", "wfc.us", "inga.nl", "lloy.gb", "ucg.it",
  "ubsn.ch", "csgn.ch", "gs.us", "nor.se", "isp.it", "ms.us", "td.ca", "ry.ca", "bbva.es",
  "cbk.de", "nab.au", "bns.ca", "cba.au", "stan.gb", "cmb1.cn", "anz.au", "wbc.au", "shgp.cn",
  "dan.dk", "sber.ru", "cmb2.cn", "bmo.ca", "itub4.br", "rsnh.jp", "nmrh.jp", "smtm.jp",
  "sbin.in", "dnb.no", "shba.se", "seba.se", "cm.ca", "bk.us", "usb.us", "bbdc4.br", "kbc.be",
  "pnc.us", "d05.sg", "pab.cn", "wrfh.kr", "dexb.be", "cof.us", "shf.kr", "swe.se", "hxb.cn",
  "ebs.at", "bmps.it", "stt.us", "sab.es", "uob.sg", "pop.es", "ibk.kr", "bbt.us", "bir.ie",
  "na.ca", "sti.us", "bp.it", "may.my", "aib.ie", "sbk.za", "axp.us", "ete.gr", "mqg.au", "ffg.jp",
  "boy.jp", "poh.fi", "fitb.us", "rf.us", "cbb.jp", "uni.it", "bcp.pr", "cimb.my", "bob.in",
  "isctr.tr", "bes.pr", "hkf.jp", "shzb.jp", "mb.it", "yfg.jp"
)

# =========================
# 1. LOAD AND FORMAT DATA
# =========================
ddly <- fread(paste0(dir.r,"ddly-data.csv"))
setnames(ddly, "dd/mm/yy", "Date")
ddly[, Date := as.Date(Date, format = "%d/%m/%y")]

# Extract time series matrix (T x p)
vol_ts <- as.matrix(ddly[, 2:97, with = FALSE])
dates <- ddly$Date
rownames(vol_ts) <- as.character(dates)
colnames(vol_ts) <- bank_codes

# =========================
# 2. PREPROCESSING
# =========================
# Safe log-transform
safe_log <- function(x) log(pmax(x, .Machine$double.eps))
log_vol_ts <- apply(vol_ts, 2, safe_log)

# Mean-zero centering
log_vol_ts_centered <- scale(log_vol_ts, center = TRUE, scale = FALSE)

# =========================
# 3. STATIONARITY TEST
# =========================
IsPowerOfTwo <- function(n) {
  if (n <= 0) return(FALSE)
  return((n & (n - 1)) == 0)
}

hwtos_func <- function(x) {
  if (IsPowerOfTwo(length(x))) locits:::hwtos2(x) else locits:::hwtos(x)
}

locits_results <- apply(log_vol_ts_centered, 2, function(x) {
  x <- na.omit(x)
  st_test <- hwtos_func(x)
  list(nreject = st_test$nreject, min_adj_pval = min(st_test$allpvals, na.rm = TRUE), st_test = st_test)
})

rej_flags <- sapply(locits_results, function(res) res$nreject > 0)
cat("Number of series rejecting second-order stationarity:", sum(rej_flags), "\n")
cat("Number of series accepted as second-order stationary:", sum(!rej_flags), "\n")

# Extract stationary series
stationary_mat <- log_vol_ts_centered[, !rej_flags, drop = FALSE]
rownames(stationary_mat) <- as.character(dates)
stationary_bank_codes <- colnames(stationary_mat)
cat("Bank codes passing the stationarity test:\n")
print(stationary_bank_codes)

# =========================
# 4. ESTIMATE LASSO-VAR
# =========================
set.seed(123)
data <- stationary_mat
colnames(data) <- stationary_bank_codes
p <- 4
K <- ncol(data)

YX <- embed(data, p + 1)
X <- YX[, (K + 1):ncol(YX)]
resids <- NULL
Bmat <- NULL

for (i in 1:K) {
  yi <- YX[, i]
  cvfit <- cv.glmnet(X, yi, alpha = 1, standardize = TRUE)
  coef_i <- as.numeric(coef(cvfit, s = "lambda.min")[-1])
  Bmat <- rbind(Bmat, coef_i)
  pred <- predict(cvfit, newx = X, s = "lambda.min")
  resids <- cbind(resids, yi - pred)
}

Sigma_hat <- crossprod(resids) / nrow(resids)

# =========================
# 5. BIC SELECTION FOR LAG ORDER
# =========================
max_p <- 10
bic_values <- numeric(max_p)
n <- nrow(stationary_mat)
K <- ncol(stationary_mat)

for (p in 1:max_p) {
  YX <- embed(stationary_mat, p + 1)
  X <- YX[, (K + 1):ncol(YX)]
  resids <- NULL
  Bmat <- NULL
  
  for (i in 1:K) {
    yi <- YX[, i]
    cvfit <- cv.glmnet(X, yi, alpha = 1, standardize = TRUE)
    coef_i <- coef(cvfit, s = "lambda.min")[-1]
    Bmat <- rbind(Bmat, as.numeric(coef_i))
    pred <- predict(cvfit, newx = X, s = "lambda.min")
    resids <- cbind(resids, yi - pred)
  }
  
  Sigma_hat <- crossprod(resids) / nrow(resids)
  logdet <- determinant(Sigma_hat, logarithm = TRUE)$modulus
  df_p <- sum(Bmat != 0)
  bic_values[p] <- as.numeric(logdet) + (log(n) * df_p / n)
}

best_p <- which.min(bic_values)
print(best_p)

# =========================
# 6. COMPUTE MA REPRESENTATION Φ_h
# =========================
Phi_array <- function(B, H, K, p) {
  A <- array(0, dim = c(K, K, p))
  for (i in 1:p) {
    A[, , i] <- matrix(B[, ((i - 1) * K + 1):(i * K)], K, K, byrow = FALSE)
  }
  Phi <- array(0, dim = c(K, K, H + 1))
  Phi[, , 1] <- diag(K)
  for (h in 1:H) {
    for (j in 1:min(h, p)) {
      Phi[, , h + 1] <- Phi[, , h + 1] + Phi[, , h - j + 1] %*% A[, , j]
    }
  }
  Phi
}

H <- 10
Phi <- Phi_array(Bmat, H = H, K = K, p = p)

# =========================
# 7. GENERALIZED FEVD
# =========================
fevd <- matrix(0, K, K)
Dinv <- diag(1 / diag(Sigma_hat + 1e-6 * diag(K)))

for (i in 1:K) {
  for (j in 1:K) {
    num <- 0
    denom <- 0
    for (h in 1:H) {
      Ah <- Phi[, , h + 1]
      row_i <- matrix(Ah[i, ], nrow = 1)
      e_j <- rep(0, K); e_j[j] <- 1
      num <- num + (row_i %*% Sigma_hat %*% Dinv %*% e_j)^2
      denom <- denom + (row_i %*% Sigma_hat %*% t(row_i))^2
    }
    fevd[i, j] <- num / denom
  }
}

# ========================================
# Normalize FEVD and Setup
# ========================================
fevd <- fevd / rowSums(fevd)
rownames(fevd) <- stationary_bank_codes
colnames(fevd) <- stationary_bank_codes

K <- nrow(fevd)

# === Percentile-Based Threshold Selection ===
# For threshold selection. Please run Threshold_selecting.R
threshold <- 0.0121
cat("Threshold set:", round(threshold, 4), "\n")

# ========================================
# Build Undirected Adjacency Matrix
# ========================================
adj_mat_undir <- matrix(0, K, K, dimnames = list(stationary_bank_codes, stationary_bank_codes))

for (i in 1:K) {
  for (j in 1:K) {
    if (i != j && (fevd[i, j] >= threshold || fevd[j, i] >= threshold)) {
      adj_mat_undir[i, j] <- 1
      adj_mat_undir[j, i] <- 1
    }
  }
}
diag(adj_mat_undir) <- 0

# ========================================
# Build Symmetric Weight Matrix
# ========================================
weight_mat <- matrix(0, K, K, dimnames = list(stationary_bank_codes, stationary_bank_codes))

for (i in 1:K) {
  for (j in 1:K) {
    if (i != j && adj_mat_undir[i, j] == 1) {
      avg_val <- mean(c(fevd[i, j], fevd[j, i]))
      weight_mat[i, j] <- avg_val
      weight_mat[j, i] <- avg_val
    }
  }
}
diag(weight_mat) <- 0

cat("Network density:", mean(adj_mat_undir), "\n")

# ========================================
# Visualize Network with igraph
# ========================================

igraph_net <- graph_from_adjacency_matrix(adj_mat_undir, mode = "undirected", diag = FALSE)
V(igraph_net)$label <- stationary_bank_codes
V(igraph_net)$name <- stationary_bank_codes

# Function to color nodes by region
get_label_color <- function(code) {
  if (grepl("\\.us$|\\.ca$|\\.br$", code)) return("red")        # Americas
  if (grepl("\\.es$|\\.fr$|\\.it$|\\.ch$|\\.gb$|\\.se$|\\.be$|\\.dk$|\\.no$", code)) return("darkblue")  # Europe
  if (grepl("\\.tr$|\\.in$|\\.kr$|\\.cn$|\\.jp$", code)) return("darkgreen") # Asia
  if (grepl("\\.au$", code)) return("darkviolet")               # Australia
  return("gray")
}

V(igraph_net)$color <- "black"
V(igraph_net)$label.color <- sapply(V(igraph_net)$name, get_label_color)

set.seed(42)
layout_net <- layout_with_fr(igraph_net)

# Plot network
pdf(file = paste0(dir.p, "Network.pdf"), width = 12, height = 10)
plot(igraph_net,
     layout = layout_net,
     vertex.size = 4,
     vertex.shape = "circle",
     vertex.label.cex = 0.7,
     vertex.label.family = "sans")
legend("topright",
       legend = c("Americas", "Europe", "Asia", "Australia"),
       text.col = c("red", "darkblue", "darkgreen", "darkviolet"),
       pch = NA, bty = "n", cex = 0.9)
dev.off()

# Filter out Japan and China nodes
keep_nodes <- !str_detect(V(igraph_net)$name, "\\.jp$|\\.cn$")
igraph_net_filtered <- induced_subgraph(igraph_net, vids = V(igraph_net)[keep_nodes])
layout_filtered <- layout_with_fr(igraph_net_filtered, niter = 3000) * 3

pdf(file = paste0(dir.p, "Network2.pdf"), width = 12, height = 10)
plot(igraph_net_filtered,
     layout = layout_filtered,
     vertex.size = 4,
     vertex.shape = "circle",
     vertex.label.cex = 0.7,
     vertex.label.family = "sans",
     edge.width = 0.5,
     edge.color = "gray70")
legend("topright",
       legend = c("Americas", "Europe", "Asia", "Australia"),
       text.col = c("red", "darkblue", "darkgreen", "darkviolet"),
       pch = NA, bty = "n", cex = 0.9)
dev.off()

cat("Graph diameter:", diameter(igraph_net, directed = FALSE), "\n")
cat("Number of edges:", gsize(igraph_net), "\n")

# ========================================
# Compute shortest paths / hops
# ========================================
hop_mat <- distances(igraph_net)
diag(hop_mat) <- NA



hop_df <- melt(hop_mat)
colnames(hop_df) <- c("From", "To", "Hops")
all_nodes <- sort(unique(c(hop_df$From, hop_df$To)))
hop_df$From <- factor(hop_df$From, levels = all_nodes)
hop_df$To <- factor(hop_df$To, levels = all_nodes)
hop_df <- hop_df[hop_df$Hops <= 5, ]
hop_df$Hops_factor <- cut(hop_df$Hops,
                          breaks = c(-Inf, 1, 2, 3, 4, 5),
                          labels = c("1 hop", "2 hops", "3 hops", "4 hops", "5 hops"),
                          right = TRUE)
hop_df <- na.omit(hop_df)

custom_colors <- c("1 hop" = "#00008B", "2 hops" = "#006400", "3 hops" = "#FF0000",
                   "4 hops" = "darkviolet", "5 hops" = "darkorange")

hop_heatmap <- ggplot(hop_df, aes(x = To, y = From, fill = Hops_factor)) +
  geom_tile(color = "lightgray") +
  scale_fill_manual(values = custom_colors, name = "Hops", drop = FALSE) +
  theme_minimal() +
  labs(title = "Number of Hops Between Nodes", x = "To", y = "From") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text = element_text(size = 6),
        panel.grid = element_blank())

print(hop_heatmap)
ggsave(filename = paste0(dir.p, "hops_heatmap.pdf"), plot = hop_heatmap, width = 8, height = 6)

# Function to compute shortest path string
get_shortest_path <- function(igraph_net, from_node, to_node) {
  stopifnot(from_node %in% V(igraph_net)$name,
            to_node %in% V(igraph_net)$name)
  sp <- shortest_paths(igraph_net, from = from_node, to = to_node, output = "vpath")$vpath[[1]]
  if (length(sp) == 0) return(paste("No path from", from_node, "to", to_node))
  paste(V(igraph_net)$name[sp], collapse = " -> ")
}

# Example usage
get_shortest_path(igraph_net, "jpm.us", "mzh.jp")
get_shortest_path(igraph_net, "bnp.fr", "cmb1.cn")
get_shortest_path(igraph_net, "cmb1.cn", "ffg.jp")

# ========================================
# GNAR Spectral Analysis
# ========================================
source("Aux_GNAR_spec.R")
gnar_net <- as.GNARnet(igraph_net)

# Plot corbit
pdf(file = paste0(dir.p, "corbit.pdf"), width = 12, height = 10)
corbit_plot_cpp(vts = data, net = gnar_net, max_lag = 20, max_stage = 5, partial = "yes")
dev.off()

# GNAR parametric estimation
# Model order selected from Corbit-plot above
alphaOrder <- 2
betaOrder <- c(1,1)
globalalpha <- TRUE

gnar_spec_param <- GNARSpec_GBNC(vts = data,
                            net = gnar_net,
                            alphaOrder = alphaOrder,
                            betaOrder = betaOrder,
                            globalalpha = globalalpha,
                            W = weight_mat)

# Functions for coherence / partial coherence
get_country <- function(code) sub(".*\\.", "", code)

compute_coherence <- function(spec_array, i, j) {
  S_ii <- spec_array[i, i, ]
  S_jj <- spec_array[j, j, ]
  S_ij <- spec_array[i, j, ]
  Mod(S_ij)^2 / (Mod(S_ii) * Mod(S_jj))
}

compute_partial_coherence <- function(spec_array, i, j) {
  n_freq <- dim(spec_array)[3]
  pcoh <- numeric(n_freq)
  for (k in 1:n_freq) {
    S <- spec_array[,,k]
    G <- solve(S)
    D_inv <- sqrt(abs(Re(diag(G))))
    if (any(D_inv == 0)) { pcoh[k] <- NA; next }
    D <- diag(1 / D_inv)
    Gamma <- -D %*% G %*% D
    pcoh[k] <- Mod(Gamma[i, j])^2
  }
  pcoh
}

# ==============================
# === 0. Setup and Inputs ===
# ==============================
countries <- sapply(stationary_bank_codes, get_country)
n_nodes <- length(stationary_bank_codes)
n_freq <- dim(gnar_spec_param)[3]
freqs <- seq(0, 0.5, length.out = n_freq)

group_palette <- c(
  "Same Country - Connected" = "#B2182B",
  "Diff Country - Connected" = "#2166AC",
  "Disconnected"             = "#4D4D4D"
)
line_types <- c(1, 2, 3)  # Recycled across groups

# ==============================
# === 1. Function to Find Pairs ===
# ==============================
adj_matrix =as.matrix(gnar_net)
find_pair <- function(c1, c2, connected = TRUE) {
  idx1 <- which(countries == c1)
  idx2 <- which(countries == c2)
  for (i in idx1) {
    for (j in idx2) {
      if (i != j && ((connected && adj_matrix[i,j] == 1) || (!connected && adj_matrix[i,j] == 0))) {
        return(c(i, j))
      }
    }
  }
  return(NULL)
}

make_pairs_df <- function(pairs_list, group_name) {
  df <- do.call(rbind, lapply(pairs_list, function(p) data.frame(t(find_pair(p[1], p[2], ifelse(length(p) == 3, p[3], TRUE))))))
  df$group <- group_name
  colnames(df) <- c("i", "j", "group")
  return(df)
}

# ==============================
# === 2. Define Pairs ===
# ==============================
countries <- sapply(stationary_bank_codes, get_country)
n_nodes <- length(stationary_bank_codes)
n_freq <- dim(gnar_spec_param)[3]
freqs <- seq(0, 0.5, length.out = n_freq)

# === 1. Explicit Pair Selection ===
find_pair <- function(c1, c2, connected = TRUE) {
  idx1 <- which(countries == c1)
  idx2 <- which(countries == c2)
  for (i in idx1) {
    for (j in idx2) {
      if (i != j && (connected && adj_matrix[i,j] == 1 || !connected && adj_matrix[i,j] == 0)) {
        return(c(i, j))
      }
    }
  }
  return(NULL)
}

same_country_connected <- rbind(
  data.frame(t(find_pair("us", "us")), group = "Same Country - Connected"),
  data.frame(t(find_pair("fr", "fr")), group = "Same Country - Connected"),
  data.frame(t(find_pair("jp", "jp")), group = "Same Country - Connected")
)

diff_country_connected <- rbind(
  data.frame(t(find_pair("us", "gb")), group = "Diff Country - Connected"),
  data.frame(t(find_pair("fr", "es")), group = "Diff Country - Connected"),
  data.frame(t(find_pair("jp", "kr")), group = "Diff Country - Connected")
)

disconnected_pairs <- rbind(
  data.frame(t(find_pair("it", "jp", connected = FALSE)), group = "Disconnected"),
  data.frame(t(find_pair("fr", "cn", connected = FALSE)), group = "Disconnected"),
  data.frame(t(find_pair("cn", "jp", connected = FALSE)), group = "Disconnected")
)

colnames(same_country_connected) <- colnames(diff_country_connected) <- colnames(disconnected_pairs) <- c("i", "j", "group")

# ==============================
# === 3. Assign Colors & Line Types ===
# ==============================
all_pairs <- rbind(same_country_connected, diff_country_connected, disconnected_pairs)
all_pairs$i <- as.integer(all_pairs$i)
all_pairs$j <- as.integer(all_pairs$j)
all_pairs$pair_label <- paste0(stationary_bank_codes[all_pairs$i], "_", stationary_bank_codes[all_pairs$j])

group_palette <- list(
  "Same Country - Connected" = "#B2182B",   # Dark Red
  "Diff Country - Connected" = "#2166AC",   # Dark Blue
  "Disconnected"             = "#4D4D4D"    # Dark Gray
)
line_types <- c(1, 2, 3)

all_pairs$Color <- NA
all_pairs$Lty <- NA

for (grp in unique(all_pairs$group)) {
  idxs <- which(all_pairs$group == grp)
  all_pairs$Color[idxs] <- group_palette[[grp]]
  all_pairs$Lty[idxs] <- line_types[seq_along(idxs)]
}
# ==============================
# === 4. Compute Spectral Quantities ===
# ==============================
mod_df <- data.frame()
phase_df <- data.frame()
pcoh_df <- data.frame()
coh_df <- data.frame()

for (row in 1:nrow(all_pairs)) {
  i <- all_pairs$i[row]
  j <- all_pairs$j[row]
  pair_label <- all_pairs$pair_label[row]
  group <- all_pairs$group[row]
  col <- all_pairs$Color[row]
  lty <- all_pairs$Lty[row]
  
  spec_ij <- gnar_spec_param[i, j, ]
  mod_vals <- Mod(spec_ij)
  phase_vals <- Arg(spec_ij)
  pcoh_vals <- compute_partial_coherence(gnar_spec_param, i, j)
  coh_vals <- compute_coherence(gnar_spec_param, i, j)
  
  mod_df <- rbind(mod_df, data.frame(Frequency = freqs, Value = mod_vals,
                                     Pair = pair_label, Group = group, Color = col, Lty = lty))
  phase_df <- rbind(phase_df, data.frame(Frequency = freqs, Value = phase_vals,
                                         Pair = pair_label, Group = group, Color = col, Lty = lty))
  pcoh_df <- rbind(pcoh_df, data.frame(Frequency = freqs, Value = pcoh_vals,
                                       Pair = pair_label, Group = group, Color = col, Lty = lty))
  coh_df <- rbind(coh_df, data.frame(Frequency = freqs, Value = coh_vals,
                                     Pair = pair_label, Group = group, Color = col, Lty = lty))
}

# ==============================
# === 5. Plotting Functions ===
# ==============================
plot_legend_inside <- function(colors, ltys, labels, x_pos = 0.95, y_pos_user, n_rows = 3, cex_legend = 1.2) {
  n_items <- length(labels)
  n_cols <- ceiling(n_items / n_rows)
  
  usr <- par("usr")
  x_legend <- grconvertX(x_pos, from = "npc", to = "user")
  y_legend <- y_pos_user
  
  line_segment_width <- 0.03 * (usr[2] - usr[1])
  label_widths <- strwidth(labels, units = "user", cex = cex_legend)
  item_widths <- label_widths + line_segment_width + 0.015 * (usr[2] - usr[1])
  
  row_items <- split(seq_len(n_items), rep(1:n_rows, each = n_cols, length.out = n_items))
  row_widths <- sapply(row_items, function(idx) sum(item_widths[idx]))
  box_width <- max(row_widths) + 0.02 * (usr[2] - usr[1])
  row_height <- 0.06 * (usr[4] - usr[3])
  box_height <- n_rows * row_height + 0.03 * (usr[4] - usr[3])
  
  x_legend <- x_legend - box_width  # right-aligned to corner
  
  # === Draw white rectangle behind legend ===
  rect(x_legend - 0.01 * (usr[2] - usr[1]), 
       y_legend - box_height,
       x_legend + box_width,
       y_legend,
       col = "white", border = "black", lwd = 1)
  
  # === Draw legend entries ===
  draw_row <- function(indices, y_pos_row) {
    x_pos_current <- x_legend
    for (i in indices) {
      segments(x_pos_current, y_pos_row,
               x_pos_current + line_segment_width, y_pos_row,
               col = colors[i], lwd = 3, lty = ltys[i])
      text(x_pos_current + line_segment_width + 0.007 * (usr[2] - usr[1]),
           y_pos_row,
           labels[i], adj = 0, cex = cex_legend)
      x_pos_current <- x_pos_current + item_widths[i]
    }
  }
  
  y_positions <- y_legend - (seq_len(n_rows) - 0.5) * row_height
  for (r in seq_len(n_rows)) {
    if (length(row_items[[r]]) > 0) {
      draw_row(row_items[[r]], y_positions[r])
    }
  }
}


plot_base <- function(df, main, ylab, n_legend_rows = 3, 
                      legend_cex = 1.2, cex_main = 2, cex_lab = 1.5, cex_axis = 1.2) {
  pairs <- unique(df$Pair)
  colors <- sapply(pairs, function(p) df$Color[df$Pair == p][1])
  ltys <- sapply(pairs, function(p) df$Lty[df$Pair == p][1])
  
  ymax <- max(df$Value, na.rm = TRUE)
  ymin <- min(df$Value, na.rm = TRUE)
  y_range <- ymax - ymin
  
  # Save user graphical settings
  op <- par(no.readonly = TRUE)
  
  # Set plot margins
  par(mar = c(4, 4, 4, 2))
  par(mgp = c(2, 0.7, 0))
  
  # Reserve space for the legend above the data
  row_height <- 0.06 * y_range
  legend_box_height <- n_legend_rows * row_height + 0.03 * y_range
  ylim_top <- ymax + legend_box_height * 2.0  # give ample space
  
  # === Base Plot ===
  plot(df$Frequency[df$Pair == pairs[1]], df$Value[df$Pair == pairs[1]], type = "n",
       ylim = c(ymin, ylim_top),
       xlab = "Frequency", ylab = ylab, main = main,
       cex.main = cex_main,
       cex.lab = cex_lab,
       cex.axis = cex_axis)
  
  # === Draw Lines ===
  for (i in seq_along(pairs)) {
    lines(df$Frequency[df$Pair == pairs[i]], df$Value[df$Pair == pairs[i]],
          col = colors[i], lwd = 2, lty = ltys[i])
  }
  
  # === Draw Legend in Upper Right ===
  plot_legend_inside(colors, ltys, pairs,
                     x_pos = 0.95,
                     y_pos_user = ylim_top - 0.05 * y_range,  # place near top
                     n_rows = n_legend_rows,
                     cex_legend = legend_cex)
  
  # Restore settings
  par(op)
}

# ==============================
# === 6. Plot All Spectral Quantities ===
# ==============================
plot_base(mod_df,
          main = expression(Modulus ~ "|" * S[ij](omega) * "|"),
          ylab = expression("|" * S[ij](omega) * "|"))


plot_base(phase_df,
          main = expression("Phase of " ~ S[ij](omega)),
          ylab = expression(Arg(S[ij](omega))))

plot_base(pcoh_df,
          main = expression("Partial Coherence " ~ gamma[ij]^2 * "(" * omega * ")"),
          ylab = expression(gamma[ij]^2 * "(" * omega * ")"))

plot_base(coh_df,
          main = expression("Coherence " ~ rho[ij]^2 * "(" * omega * ")"),
          ylab = expression(rho[ij]^2 * "(" * omega * ")"))



