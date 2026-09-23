# ==============================================================================
# credit_qres_envelope.R
# Three quantile-residual panels for the final IPL-TF fit, each with a
# simulated envelope obtained by REFITTING (Atkinson's envelope: minimum and
# maximum over B replicates simulated from the fitted model).
#
#   (a) discrete part   : u ~ U(1-alpha,1) if y = c, U(0,1-alpha) otherwise
#   (b) continuous part : Phi^{-1}[F_PL(y)] on the observations with y < 1
#   (c) both            : residual of the mixture, u ~ U(1-alpha,1) at the mass
#
# The randomized quantile residual is used in the discrete part too: a Pearson
# residual takes two values per fitted alpha_i and its standardized version
# needs (Z'QZ)^{-1}, which does not exist with p = 348 > n = 320.
#
# Requires final_fit/credit_final_fit.RData (run credit_final_fit.R first).
# B defaults to 100; override with the B_ENV environment variable.
# Supplementary diagnostic: not part of the text of dissertation v1.0.
# ==============================================================================

library(data.table); library(mboost); library(gamboostLSS)
library(gamlss.dist); library(EnvStats); library(parallel)

app_path <- Sys.getenv("PLBOOST_OUT", unset = ".")   # must contain data/
out_path <- file.path(app_path, "final_fit")
library(PLboost); library(PLreg); library(gamboostLSS)   # PL / IPL families and distributions
load(file.path(out_path, "credit_final_fit.RData"))   # fit, mu_h, sg_h, lb_h, al_h, zeta_best, opt, df_lgd

FAM <- "TF"; C_INFL <- 1; nu <- 0.001
B <- as.integer(Sys.getenv("B_ENV", "100"))
y <- df_lgd$lgd
n <- length(y)
cat("n =", n, " B =", B, " zeta =", zeta_best, "\n")

fit_at <- function(dat) {
  assign("fam", IPL(family = FAM, stabilization = "MAD", c = C_INFL, zeta = zeta_best),
         envir = globalenv())
  glmboostLSS(formula = list(mu = lgd ~ ., sigma = lgd ~ .,
                             lambda = lgd ~ 1, alpha = lgd ~ .),
              data = dat, families = fam,
              control = boost_control(mstop = c(opt["mu"], opt["sigma"],
                                                opt["lambda"], opt["alpha"]),
                                      nu = nu, trace = FALSE))
}

# --- os tres residuos, dado um ajuste e um vetor de resposta ------------------
# qnorm() reclama se a probabilidade sair de [0,1] por arredondamento numerico
qn <- function(p) qnorm(pmin(pmax(p, 1e-12), 1 - 1e-12))

res3 <- function(yy, mu, sg, lb, al) {
  at_c <- (yy == C_INFL)
  # (a) discreto: quantilico aleatorizado da Bernoulli do indicador
  u_d <- ifelse(at_c, runif(length(yy), 1 - al, 1), runif(length(yy), 0, 1 - al))
  r_d <- qn(u_d)
  # (b) continuo: so as observacoes com y < 1
  r_c <- qn(pPL(yy[!at_c], mu = mu[!at_c], sigma = sg[!at_c],
                   lambda = lb, zeta = zeta_best, family = FAM))
  # (c) mistura
  r_j <- numeric(length(yy))
  r_j[!at_c] <- qn(pIPL(yy[!at_c], alpha = al[!at_c], mu = mu[!at_c],
                           sigma = sg[!at_c], lambda = lb, zeta = zeta_best,
                           family = FAM, c = C_INFL))
  r_j[at_c]  <- qn(runif(sum(at_c), 1 - al[at_c], 1))
  lapply(list(disc = r_d, cont = r_c, joint = r_j),
         function(v) pmin(pmax(v, -8), 8))
}

set.seed(20260812L)
obs <- res3(y, mu_h, sg_h, lb_h, al_h)
cat("observado | disc sd", round(sd(obs$disc),3),
    "| cont sd", round(sd(obs$cont),3),
    "| joint sd", round(sd(obs$joint),3), "\n")

# --- envelope por reajuste ----------------------------------------------------
sim <- list(disc = matrix(NA_real_, n, B),
            cont = matrix(NA_real_, length(obs$cont), B),
            joint = matrix(NA_real_, n, B))
t0 <- Sys.time(); ok <- 0
for (b in seq_len(B)) {
  ys <- rIPL(n = n, alpha = al_h, mu = mu_h, sigma = sg_h,
             lambda = lb_h, zeta = zeta_best, family = FAM, c = C_INFL)
  ys <- pmin(pmax(ys, 1e-6), 1)
  d2 <- copy(df_lgd); d2[, lgd := ys]
  m2 <- try(fit_at(d2), silent = TRUE)
  if (inherits(m2, "try-error")) { cat("  b", b, "falhou\n"); next }
  mu2 <- as.numeric(predict(m2$mu,     newdata = d2, type = "response"))
  sg2 <- as.numeric(predict(m2$sigma,  newdata = d2, type = "response"))
  lb2 <- as.numeric(predict(m2$lambda, newdata = d2, type = "response"))[1]
  al2 <- as.numeric(predict(m2$alpha,  newdata = d2, type = "response"))
  r2  <- res3(ys, mu2, sg2, lb2, al2)
  sim$disc[, b]  <- sort(r2$disc)
  sim$joint[, b] <- sort(r2$joint)
  rc <- sort(r2$cont)
  sim$cont[, b] <- if (length(rc) == nrow(sim$cont)) rc else
    approx(seq_along(rc), rc, xout = seq(1, length(rc),
           length.out = nrow(sim$cont)))$y
  ok <- ok + 1
  rm(m2, d2); gc()
  if (b %% 10 == 0) cat("  b", b, "de", B, "|",
      round(difftime(Sys.time(), t0, units = "min"), 1), "min\n")
}
cat("reajustes bem sucedidos:", ok, "de", B, "\n")

# --- graficos -----------------------------------------------------------------
painel <- function(r, S, arq, ylab) {
  r <- sort(r); m <- length(r)
  # envelope no sentido classico (Atkinson): minimo e maximo das replicas.
  # O intervalo pontual de 95% fica estreito demais com B grande e acusa
  # discrepancia onde nao ha.
  lo <- apply(S, 1, min, na.rm = TRUE)
  hi <- apply(S, 1, max, na.rm = TRUE)
  mi <- apply(S, 1, median, na.rm = TRUE)
  qt <- qnorm(ppoints(m))
  png(file.path(out_path, arq), width = 1500, height = 1400, res = 240)
  par(mar = c(4.2, 4.4, 1, 1))
  plot(qt, r, pch = 20, cex = .7, xlab = "Theoretical quantiles",
       ylab = ylab, ylim = range(c(r, lo, hi), na.rm = TRUE))
  lines(qt, lo); lines(qt, hi); lines(qt, mi, lty = 2)
  dev.off()
  sum(r < lo | r > hi)
}
f1 <- painel(obs$disc,  sim$disc,  "credit_qres_a_discrete.png",   "Quantile residuals")
f2 <- painel(obs$cont,  sim$cont,  "credit_qres_b_continuous.png", "Quantile residuals")
f3 <- painel(obs$joint, sim$joint, "credit_qres_c_joint.png",      "Quantile residuals")
cat("fora do envelope | (a)", f1, "de", n,
    "| (b)", f2, "de", length(obs$cont), "| (c)", f3, "de", n, "\n")

save(obs, sim, file = file.path(out_path, "credit_qres_envelope.RData"))
cat("PRONTO em", round(difftime(Sys.time(), t0, units = "min"), 1), "min\n")
