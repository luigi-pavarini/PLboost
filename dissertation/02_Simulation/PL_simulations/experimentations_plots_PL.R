# Boosted PL regression

rm(list = ls())

# ==============================================================================
# SETUP
# ==============================================================================

# Needs to be settled!!!
base_path  <- "..."

# ---- Family selector --------------------------------------------------------
# Options: "NO", "TF", "LO", "SN", "PE", "Hyp", "SLASH"
fam <- "Hyp"

# ---- Image export settings --------------------------------------------------
img_width  <- 2000    # pixels
img_height <- 2000    # pixels
img_res    <- 330    # dpi

# ---- Derived paths ----------------------------------------------------------
fam_folder  <- switch(fam,
                      NO    = "1_NO",
                      TF    = "2_TF",
                      LO    = "3_LO",
                      SN    = "4_SN",
                      PE    = "5_PE",
                      Hyp   = "6_Hyp",
                      SLASH = "7_SLASH"
)

mstop_tag   <- if (fam == "SN") "mstop500" else "mstop100"
rdata_file  <- paste0("sim_PL_M100_", mstop_tag, "_", fam, ".RData")
save_path   <- file.path(base_path, "PL_simulations", fam_folder)

source(file.path(base_path, "essentials", "families_PL.R"))

M <- 100
n <- 800
p <- 1000

load(file.path(base_path, "PL_simulations", rdata_file))

# ==============================================================================
# VISUALIZATION
# Column order: (Intercept)[1] | V1[2]...V1000[1001]
# mu    true: V1[2], V2[3], V3[4], V4[5]           non-inf: [c(6:1001),]
# sigma true: V3[4], V4[5], V5[6], V6[7]           non-inf: [c(2:3, 8:1001),]
# ==============================================================================

# ---- Stopping iterations ----------------------------------------------------
png(file.path(save_path, paste0(fam, "_mstop.png")),
    width = img_width, height = img_height, res = img_res)

boxplot(sim$m_stop[1,], sim$m_stop[2,], sim$m_stop[3,],
        names = c("mu", "sigma", "lambda"),
        ylab = "Stopping iterations")

dev.off()

# ---- Mu coefficients --------------------------------------------------------
png(file.path(save_path, paste0(fam, "_mu_coefs.png")),
    width = img_width, height = img_height, res = img_res)

boxplot(sim$est_mu[2,], sim$est_mu[3,], sim$est_mu[4,],
        sim$est_mu[5,], c(sim$est_mu[c(6:1001),]),
        ylim = c(-1.2, 0.5), xaxt = "n", las = 1,
        ylab = expression("Coefficients for" ~ eta[mu]),
        main = ""
)
abline(h = 0, lty = 2, col = "gray")
axis(1, 1, expression(X[1]))
axis(1, 2, expression(X[2]))
axis(1, 3, expression(X[3]))
axis(1, 4, expression(X[4]))
axis(1, 5, "non-inf.")
segments(x0 = 0.5, x1 = 1.5, y0 = -0.99, y1 = -0.99, lty = 2, col = "red4")
segments(x0 = 1.5, x1 = 2.5, y0 =  0.49, y1 =  0.49, lty = 2, col = "red4")
segments(x0 = 2.5, x1 = 3.5, y0 = -1.17, y1 = -1.17, lty = 2, col = "red4")
segments(x0 = 3.5, x1 = 4.5, y0 =  0.38, y1 =  0.38, lty = 2, col = "red4")
segments(x0 = 4.5, x1 = 5.5, y0 =  0,    y1 =  0,    lty = 2, col = "red4")

dev.off()

# ---- Sigma coefficients -----------------------------------------------------
png(file.path(save_path, paste0(fam, "_sigma_coefs.png")),
    width = img_width, height = img_height, res = img_res)

boxplot(sim$est_sigma[4,], sim$est_sigma[5,], sim$est_sigma[6,],
        sim$est_sigma[7,], c(sim$est_sigma[c(2:3, 8:1001),]),
        ylim = c(-0.70, 1.00), xaxt = "n", las = 1,
        ylab = expression("Coefficients for" ~ eta[sigma]),
        main = ""
)
abline(h = 0, lty = 2, col = "gray")
axis(1, 1, expression(X[3]))
axis(1, 2, expression(X[4]))
axis(1, 3, expression(X[5]))
axis(1, 4, expression(X[6]))
axis(1, 5, "non-inf.")
segments(x0 = 0.5, x1 = 1.5, y0 = -0.41, y1 = -0.41, lty = 2, col = "red4")
segments(x0 = 1.5, x1 = 2.5, y0 = -0.24, y1 = -0.24, lty = 2, col = "red4")
segments(x0 = 2.5, x1 = 3.5, y0 =  0.81, y1 =  0.81, lty = 2, col = "red4")
segments(x0 = 3.5, x1 = 4.5, y0 =  0.70, y1 =  0.70, lty = 2, col = "red4")
segments(x0 = 4.5, x1 = 5.5, y0 =  0,    y1 =  0,    lty = 2, col = "red4")

dev.off()

# ---- Lambda -----------------------------------------------------------------
png(file.path(save_path, paste0(fam, "_lambda_coefs.png")),
    width = img_width, height = img_height, res = img_res)

boxplot(exp(sim$est_lambda[1,]), ylim = c(0, 10), xaxt = "n", las = 1,
        ylab = expression("Coefficients for" ~ lambda),
        main = ""
)
abline(h = 0, lty = 2, col = "gray")
axis(1, 1, "exp(intercept)")
segments(x0 = 0.5, x1 = 1.5, y0 = 1.2, y1 = 1.2, lty = 2, col = "red4")

dev.off()

# ==============================================================================
# LOLLIPOP PLOTS
# - Each column (run) is normalized to proportions first (don't mix samples)
# - Then average proportion across M runs is computed
# - Error bars show +/- sd_mult * sd across runs  (default: 1)
# - Stems: solid gray; error bars: dashed steel blue
# - Top 15 variables shown individually; rest collapsed into "others"
# ==============================================================================

lollipop_sel <- function(freq_mat, top_n = 15, param_label = "mu", sd_mult = 1) {
  
  # Step 1: normalize each run (column) to proportions — avoids mixing samples
  col_totals <- colSums(freq_mat)
  col_totals[col_totals == 0] <- 1  # guard against zero-sum columns
  prop_mat <- sweep(freq_mat, 2, col_totals, FUN = "/")
  
  # Drop intercept row
  prop_mat <- prop_mat[-1, ]
  
  # Step 2: average proportion and sd across M runs
  avg_prop <- rowSums(prop_mat) / M
  sd_prop  <- apply(prop_mat, 1, sd)
  
  # Step 3: select top_n by average proportion
  ord       <- order(avg_prop, decreasing = TRUE)[seq_len(top_n)]
  top_avg   <- avg_prop[ord]
  top_sd    <- sd_prop[ord]
  top_names <- names(top_avg)
  
  # Step 4: collapse remaining into "others"
  others_avg <- sum(avg_prop[-ord])
  others_sd  <- sqrt(sum(sd_prop[-ord]^2))  # pooled sd for "others"
  
  p_plot   <- c(top_avg, others_avg)
  sd_plot  <- c(top_sd,  others_sd)
  lbl_plot <- c(top_names, "others")
  n_plot   <- top_n + 1
  
  # Step 5: axis label
  xlab_expr <- switch(param_label,
                      mu    = expression("Avg. selection proportion across M runs" ~ eta[mu]),
                      sigma = expression("Avg. selection proportion across M runs" ~ eta[sigma])
  )
  
  # Step 6: plot
  old_par <- par(mar = c(4, 6, 3, 2))
  
  x_max <- max(p_plot + sd_mult * sd_plot, na.rm = TRUE) * 1.15
  
  plot(p_plot, seq_len(n_plot),
       ylim = c(n_plot, 1),
       xlim = c(0, x_max),
       yaxt = "n",
       xlab = xlab_expr,
       ylab = "",
       main = "",
       pch  = 19, col = "gray30", cex = 1.3)
  
  # Lollipop stems: light dashed gray
  segments(x0 = 0, x1 = p_plot,
           y0 = seq_len(n_plot), y1 = seq_len(n_plot),
           col = "gray80", lwd = 0.8, lty = 2)
  
  # Error bars: sd_mult * sd, solid steel blue
  x_lo <- pmax(p_plot - sd_mult * sd_plot, 0)
  x_hi <- p_plot + sd_mult * sd_plot
  cap  <- 0.18  # vertical cap half-height in plot units
  
  # Horizontal bar
  segments(x0 = x_lo, x1 = x_hi,
           y0 = seq_len(n_plot), y1 = seq_len(n_plot),
           col = "#4682B4", lwd = 1.4)
  # Left cap
  segments(x0 = x_lo, x1 = x_lo,
           y0 = seq_len(n_plot) - cap, y1 = seq_len(n_plot) + cap,
           col = "#4682B4", lwd = 1.4)
  # Right cap
  segments(x0 = x_hi, x1 = x_hi,
           y0 = seq_len(n_plot) - cap, y1 = seq_len(n_plot) + cap,
           col = "#4682B4", lwd = 1.4)
  
  axis(2, at = seq_len(n_plot), labels = lbl_plot, las = 1, cex.axis = 0.85)
  abline(v = 0, col = "gray60", lty = 2)
  par(old_par)
}

# ---- Lollipop: Mu -----------------------------------------------------------
png(file.path(save_path, paste0(fam, "_lollipop_mu_coefs.png")),
    width = img_width, height = img_height, res = img_res)

lollipop_sel(sim$freq_sel_mu,    top_n = 15, param_label = "mu",    sd_mult = 1)

dev.off()

# ---- Lollipop: Sigma --------------------------------------------------------
png(file.path(save_path, paste0(fam, "_lollipop_sigma_coefs.png")),
    width = img_width, height = img_height, res = img_res)

lollipop_sel(sim$freq_sel_sigma, top_n = 15, param_label = "sigma", sd_mult = 1)

dev.off()

message("All plots saved to: ", save_path)