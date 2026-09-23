# ==============================================================================
# credit_beoi_qres.R
# The same quantile-residual panels for the one-inflated beta baseline (BEOI),
# fitted by boosting on the full sample, for comparison with the IPL-TF panels.
#
# gamlss.dist::BEOI parameterization (from dBEOI):
#   log.lik = ifelse(x == 1, log(nu), log(1-nu) + dbeta(x, mu*sigma, (1-mu)*sigma))
# so nu is P(Y = 1) and F(1-) = 1 - nu.
#
#   (a) discrete   : u ~ U(1-nu, 1) if y = 1, U(0, 1-nu) otherwise
#   (b) continuous : Phi^{-1}[pbeta(y; mu*sigma, (1-mu)*sigma)] on y < 1
#   (c) joint      : Phi^{-1}[pBEOI(y)] on y < 1, u ~ U(1-nu, 1) at the mass
#
# Same recipe as the IPL: mstop by cvrisk() on the full sample, Atkinson
# envelope by refitting. B defaults to 200; override with B_ENV.
# Supplementary diagnostic: not part of the text of dissertation v1.0.
# ==============================================================================

library(data.table); library(mboost); library(gamboostLSS)
library(gamlss.dist); library(parallel)

app_path <- Sys.getenv("PLBOOST_OUT", unset = ".")   # must contain data/
out_path <- file.path(app_path, "final_fit")

# ---- dados, identico ao credit_run.r --------------------------------
load(file.path(app_path, "data", "df_lgd_all.RData"))
addr_state_cols <- grep("^addr_state_", names(df_lgd), value = TRUE)
df_lgd <- df_lgd[addr_state_VT == 1]
df_lgd[, (addr_state_cols) := NULL]
pred_cols <- setdiff(names(df_lgd), "lgd")
setnames(df_lgd, old = pred_cols, new = paste0("V", seq_along(pred_cols)))
y <- df_lgd$lgd; n <- length(y)
cat("n =", n, " p =", length(pred_cols), " | massa em 1:", sum(y == 1), "\n")

# Hiperparametros IDENTICOS aos do BEOI em credit_run.r, que e' a mesma origem de onde sairam os do IPL-TF
# usados em credit_final_fit.R. Nada aqui foi ajustado para este script.
LR <- 0.001                      # learning rate do boosting (nao confundir com o nu do BEOI)
mstop_init <- 500; m_stop_cv <- 800
B <- as.integer(Sys.getenv("B_ENV", "200"))

grid_beoi <- make.grid(max = c(mu = m_stop_cv, sigma = m_stop_cv, nu = m_stop_cv),
                       min = mstop_init, length.out = 2, dense_mu_grid = TRUE)

# `fam` no ambiente global: o cvrisk() reavalia a call noutro ambiente
fit_beoi <- function(dat, mstop_vec) {
  assign("fam", as.families(gamlss.dist::BEOI(), stabilization = "MAD"),
         envir = globalenv())
  glmboostLSS(formula = list(mu = lgd ~ ., sigma = lgd ~ ., nu = lgd ~ .),
              data = dat, families = fam,
              control = boost_control(mstop = mstop_vec, nu = LR, trace = FALSE))
}
qn <- function(p) qnorm(pmin(pmax(p, 1e-12), 1 - 1e-12))

res3 <- function(yy, mu, sg, nuv) {
  at1 <- (yy == 1)
  a <- mu * sg; b <- (1 - mu) * sg
  r_d <- qn(ifelse(at1, runif(length(yy), 1 - nuv, 1), runif(length(yy), 0, 1 - nuv)))
  r_c <- qn(pbeta(yy[!at1], shape1 = a[!at1], shape2 = b[!at1]))
  r_j <- numeric(length(yy))
  r_j[!at1] <- qn(pBEOI(yy[!at1], mu = mu[!at1], sigma = sg[!at1], nu = nuv[!at1]))
  r_j[at1]  <- qn(runif(sum(at1), 1 - nuv[at1], 1))
  lapply(list(disc = r_d, cont = r_c, joint = r_j),
         function(v) pmin(pmax(v, -8), 8))
}

cat("\n[1] mstop por cvrisk na base completa\n")
m0 <- fit_beoi(df_lgd, rep(mstop_init, 3))
set.seed(20260813L)
opt <- mstop(cvrisk(m0, grid = grid_beoi, folds = mboost::cv(model.weights(m0)),
                    papply = mclapply))
print(opt)

cat("\n[2] ajuste final\n")
fit  <- fit_beoi(df_lgd, c(opt["mu"], opt["sigma"], opt["nu"]))
mu_h <- as.numeric(predict(fit$mu,    newdata = df_lgd, type = "response"))
sg_h <- as.numeric(predict(fit$sigma, newdata = df_lgd, type = "response"))
nu_h <- as.numeric(predict(fit$nu,    newdata = df_lgd, type = "response"))
cat("   nu_hat (P(Y=1)) em [", round(min(nu_h),3), ",", round(max(nu_h),3),
    "] | media", round(mean(nu_h),3), " | observado", round(mean(y == 1),3), "\n")

set.seed(20260813L)
obs <- res3(y, mu_h, sg_h, nu_h)
for (k in names(obs))
  cat(sprintf("   %-6s media %6.3f  sd %5.3f\n", k, mean(obs[[k]]), sd(obs[[k]])))

cat("\n[3] envelope por reajuste, B =", B, "\n")
sim <- list(disc = matrix(NA_real_, n, B),
            cont = matrix(NA_real_, length(obs$cont), B),
            joint = matrix(NA_real_, n, B))
t0 <- Sys.time(); ok <- 0
for (bb in seq_len(B)) {
  ys <- rBEOI(n, mu = mu_h, sigma = sg_h, nu = nu_h)
  ys <- pmin(pmax(ys, 1e-6), 1)
  d2 <- copy(df_lgd); d2[, lgd := ys]
  m2 <- try(fit_beoi(d2, c(opt["mu"], opt["sigma"], opt["nu"])), silent = TRUE)
  if (inherits(m2, "try-error")) next
  r2 <- res3(ys,
             as.numeric(predict(m2$mu,    newdata = d2, type = "response")),
             as.numeric(predict(m2$sigma, newdata = d2, type = "response")),
             as.numeric(predict(m2$nu,    newdata = d2, type = "response")))
  sim$disc[, bb]  <- sort(r2$disc)
  sim$joint[, bb] <- sort(r2$joint)
  rc <- sort(r2$cont)
  sim$cont[, bb] <- if (length(rc) == nrow(sim$cont)) rc else
    approx(seq_along(rc), rc, xout = seq(1, length(rc), length.out = nrow(sim$cont)))$y
  ok <- ok + 1; rm(m2, d2); gc()
  if (bb %% 25 == 0) cat("   b", bb, "|", round(difftime(Sys.time(), t0, units="min"),1), "min\n")
}
cat("   reajustes ok:", ok, "de", B, "\n")

painel <- function(r, S, arq, tit) {
  r <- sort(r); m <- length(r)
  lo <- apply(S, 1, min, na.rm = TRUE); hi <- apply(S, 1, max, na.rm = TRUE)
  mi <- apply(S, 1, median, na.rm = TRUE); qt <- qnorm(ppoints(m))
  png(file.path(out_path, arq), width = 1500, height = 1400, res = 240)
  par(mar = c(4.2, 4.4, 1, 1))
  plot(qt, r, pch = 20, cex = .7, xlab = "Theoretical quantiles",
       ylab = "Quantile residuals", ylim = range(c(r, lo, hi), na.rm = TRUE))
  lines(qt, lo); lines(qt, hi); lines(qt, mi, lty = 2); dev.off()
  cat(sprintf("   %-12s fora do envelope: %d de %d\n", tit, sum(r < lo | r > hi), m))
}
cat("\n[4] paineis\n")
painel(obs$disc,  sim$disc,  "beoi_qres_a_discrete.png",   "(a) discreta")
painel(obs$cont,  sim$cont,  "beoi_qres_b_continuous.png", "(b) continua")
painel(obs$joint, sim$joint, "beoi_qres_c_joint.png",      "(c) conjunta")

# `fit` e' o modelo BEOI final, no mstop escolhido pelo cvrisk acima.
# Ele e' salvo para permitir diagnosticos posteriores (por exemplo, selecao
# de variaveis em fit$nu) sem reajustar o boosting.
save(fit, obs, sim, mu_h, sg_h, nu_h, opt,
     file = file.path(out_path, "credit_beoi_qres.RData"))
cat("\nPRONTO em", round(difftime(Sys.time(), t0, units="min"),1), "min\n")
