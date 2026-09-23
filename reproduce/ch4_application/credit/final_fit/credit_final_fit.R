# ==============================================================================
# credit_final_fit.R
# Final IPL-TF fit on the full Vermont sample (n = 320) and quantile-residual
# diagnostics.
#
# Cross-validation fits 50 models and answers out-of-sample performance; model
# adequacy is a property of ONE fit on all the data, so this script refits once:
#   1. zeta chosen by Upsilon on the full sample
#   2. mstop chosen by cvrisk() on the full sample
#   3. final fit and randomized quantile residuals (Dunn & Smyth, 1996),
#      with the point mass at c = 1:
#        y_i = 1       -> Phi^{-1}(u_i),  u_i ~ U(1 - alpha_i, 1)
#        y_i in (0,1)  -> Phi^{-1}(F^{(c)}(y_i; alpha_i, mu_i, sigma_i, lambda))
#
# Outputs in $PLBOOST_OUT/final_fit/:
#   credit_qres_qq.png, credit_qres_index.png, credit_final_estimates.md,
#   credit_final_fit.RData (read by credit_qres_envelope.R)
# Supplementary diagnostic: not part of the text of dissertation v1.0.
# ==============================================================================

library(data.table); library(mboost); library(gamboostLSS)
library(gamlss.dist); library(EnvStats); library(parallel)

app_path <- Sys.getenv("PLBOOST_OUT", unset = ".")   # must contain data/
out_path <- file.path(app_path, "final_fit")
dir.create(out_path, showWarnings = FALSE)

library(PLboost); library(PLreg); library(gamboostLSS)   # PL / IPL families and distributions

# ---- dados, identico ao credit_run.r --------------------------------
load(file.path(app_path, "data", "df_lgd_all.RData"))
addr_state_cols <- grep("^addr_state_", names(df_lgd), value = TRUE)
df_lgd <- df_lgd[addr_state_VT == 1]
df_lgd[, (addr_state_cols) := NULL]
pred_cols <- setdiff(names(df_lgd), "lgd")
var_dict  <- data.table(new_name = paste0("V", seq_along(pred_cols)),
                        orig_name = pred_cols)
setnames(df_lgd, old = pred_cols, new = var_dict$new_name)
cat("n =", nrow(df_lgd), " p =", length(pred_cols),
    " | massa em 1:", sum(df_lgd$lgd == 1), "\n")

FAM <- "TF"; C_INFL <- 1
mstop_init <- 500; m_stop_cv <- 800; nu <- 0.001
zeta_grid <- c(3, 4, 5, 6, 8, 10, 15, 20)

grid_ipl <- make.grid(max = c(mu = m_stop_cv, sigma = m_stop_cv,
                              lambda = m_stop_cv, alpha = m_stop_cv),
                      min = mstop_init, length.out = 2, dense_mu_grid = TRUE)

upsilon_ipl <- function(y, mu, sigma, lambda, zeta, family) {
  cont <- y > 0 & y < 1
  fda  <- sort(pPL(y[cont], mu = mu[cont], sigma = sigma[cont],
                   lambda = lambda, zeta = zeta, family = family))
  n <- length(fda)
  mean(abs(qnorm(fda) - EnvStats::evNormOrdStats(n = n)))
}

# ATENCAO: `fam` precisa viver no ambiente global. O cvrisk() reavalia
# `families = fam` a partir da call guardada, noutro ambiente, e nao enxerga
# variavel local -- da "object 'fam' not found". O script de CV faz igual.
fit_ipl <- function(dat, zeta, mstop_vec, trace = FALSE) {
  assign("fam", IPL(family = FAM, stabilization = "MAD", c = C_INFL, zeta = zeta),
         envir = globalenv())
  glmboostLSS(formula = list(mu = lgd ~ ., sigma = lgd ~ .,
                             lambda = lgd ~ 1, alpha = lgd ~ .),
              data = dat, families = fam,
              control = boost_control(mstop = mstop_vec, nu = nu, trace = trace))
}

# ---- 1. zeta por Upsilon na base completa ------------------------------------
cat("\n[1] escolhendo zeta na base completa\n")
ups <- rep(NA_real_, length(zeta_grid))
for (z in seq_along(zeta_grid)) {
  m <- try(fit_ipl(df_lgd, zeta_grid[z], rep(mstop_init, 4)), silent = TRUE)
  if (inherits(m, "try-error")) next
  ups[z] <- upsilon_ipl(df_lgd$lgd,
    as.numeric(predict(m$mu,     newdata = df_lgd, type = "response")),
    as.numeric(predict(m$sigma,  newdata = df_lgd, type = "response")),
    as.numeric(predict(m$lambda, newdata = df_lgd, type = "response"))[1],
    zeta_grid[z], FAM)
  cat("   zeta =", zeta_grid[z], " Upsilon =", round(ups[z], 4), "\n")
  rm(m); gc()
}
zeta_best <- zeta_grid[which.min(ups)]
cat("   -> zeta escolhido:", zeta_best, "\n")

# ---- 2. mstop por CV na base completa ----------------------------------------
cat("\n[2] escolhendo mstop por cvrisk\n")
m0 <- fit_ipl(df_lgd, zeta_best, rep(mstop_init, 4))
set.seed(20260812L)
cvr  <- cvrisk(m0, grid = grid_ipl, folds = mboost::cv(model.weights(m0)),
               papply = mclapply)
opt  <- mstop(cvr)
print(opt)

# ---- 3. ajuste final ---------------------------------------------------------
cat("\n[3] ajuste final\n")
fit <- fit_ipl(df_lgd, zeta_best,
               c(opt["mu"], opt["sigma"], opt["lambda"], opt["alpha"]))

mu_h  <- as.numeric(predict(fit$mu,     newdata = df_lgd, type = "response"))
sg_h  <- as.numeric(predict(fit$sigma,  newdata = df_lgd, type = "response"))
lb_h  <- as.numeric(predict(fit$lambda, newdata = df_lgd, type = "response"))[1]
al_h  <- as.numeric(predict(fit$alpha,  newdata = df_lgd, type = "response"))

# ---- 4. residuo quantilico ---------------------------------------------------
y <- df_lgd$lgd
set.seed(20260812L)
rq <- numeric(length(y))
at_c <- (y == C_INFL)
rq[!at_c] <- qnorm(pIPL(y[!at_c], alpha = al_h[!at_c], mu = mu_h[!at_c],
                        sigma = sg_h[!at_c], lambda = lb_h, zeta = zeta_best,
                        family = FAM, c = C_INFL))
rq[at_c]  <- qnorm(runif(sum(at_c), 1 - al_h[at_c], 1))   # c = 1
rq <- pmin(pmax(rq, -8), 8)
cat("   r^q: media", round(mean(rq), 3), " sd", round(sd(rq), 3),
    " | fora de [-3,3]:", sum(abs(rq) > 3), "de", length(rq), "\n")

# ---- 5. graficos -------------------------------------------------------------
n <- length(rq); B <- 1000
env <- replicate(B, sort(rnorm(n)))
lo <- apply(env, 1, quantile, 0.025); hi <- apply(env, 1, quantile, 0.975)
mid <- apply(env, 1, median)
qt <- qnorm(ppoints(n))

png(file.path(out_path, "credit_qres_qq.png"), width = 1500, height = 1400, res = 240)
par(mar = c(4.2, 4.4, 1, 1))
plot(qt, sort(rq), pch = 3, cex = .55, xlab = "Theoretical quantiles",
     ylab = "Quantile residuals", ylim = range(c(rq, lo, hi)))
lines(qt, lo, lty = 1); lines(qt, hi, lty = 1); lines(qt, mid, lty = 2)
dev.off()

png(file.path(out_path, "credit_qres_index.png"), width = 1500, height = 1200, res = 240)
par(mar = c(4.2, 4.4, 1, 1))
plot(rq, pch = 3, cex = .55, xlab = "Observation indices", ylab = "Quantile residuals")
abline(h = c(-3, 0, 3), lty = c(2, 1, 2))
dev.off()

# ---- 6. tabela das estimativas ----------------------------------------------
coef_tab <- function(sub, nome) {
  cf <- coef(fit[[sub]], off2int = TRUE)
  if (!length(cf)) return(NULL)
  d <- data.table(Submodel = nome, Variable = names(cf), Estimate = as.numeric(cf))
  d[, Variable := ifelse(Variable %in% var_dict$new_name,
        var_dict$orig_name[match(Variable, var_dict$new_name)], Variable)]
  d[order(-abs(Estimate))]
}
tab <- rbindlist(list(coef_tab("mu","eta_mu"), coef_tab("sigma","eta_sigma"),
                      coef_tab("alpha","eta_alpha")), use.names = TRUE)

md <- c("# Ajuste final IPL-TF, base completa (n = 320)", "",
        sprintf("- zeta = %s, escolhido por Upsilon na base completa", zeta_best),
        sprintf("- mstop = mu %s, sigma %s, lambda %s, alpha %s",
                opt["mu"], opt["sigma"], opt["lambda"], opt["alpha"]),
        sprintf("- lambda_hat = %.4f", lb_h),
        sprintf("- Upsilon (base completa) = %.4f",
                upsilon_ipl(y, mu_h, sg_h, lb_h, zeta_best, FAM)),
        sprintf("- r^q: media %.3f, sd %.3f, |r| > 3 em %d de %d",
                mean(rq), sd(rq), sum(abs(rq) > 3), length(rq)), "",
        "| Submodel | Variable | Estimate |", "|---|---|---|",
        sprintf("| %s | %s | %.4f |", tab$Submodel, tab$Variable, tab$Estimate))
writeLines(md, file.path(out_path, "credit_final_estimates.md"))

save(fit, rq, mu_h, sg_h, lb_h, al_h, zeta_best, opt, var_dict, df_lgd,
     file = file.path(out_path, "credit_final_fit.RData"))
cat("\nPRONTO. Saidas em", out_path, "\n")
