# ==============================================================================
# validation_multirep.R
#
# Multi-replication version of the Chapter 2 validation (boosting vs MLE).
# Each scenario is now repeated over M Monte Carlo replications and the tables
# report mean (standard deviation) across replications instead of a single fit.
#
# Same models, coefficient vectors and links as the single-sample validation;
# MLE with lambda held fixed at the true value. Seeded (set.seed(123)), so a
# rerun reproduces the dissertation tables exactly.
#
# Writes to $PLBOOST_OUT (default: the working directory):
#   val_pl.tex, val_ipl_continuous.tex, val_ipl_alpha.tex   (cells are "mean (sd)")
#   validation_results.md, validation_results.csv            (raw numbers)
#
# Requires: PLboost, gamboostLSS, mboost, faux, PLreg, Formula.
# Run:  Rscript validation_multirep.R
# ==============================================================================

suppressMessages({
  library(gamboostLSS); library(mboost); library(faux)
  library(PLreg); library(Formula)
})

# ---- output directory ---------------------------------------------------------
library(PLboost)
this_dir <- Sys.getenv("PLBOOST_OUT", unset = ".")   # tables and csv are written here
dir.create(this_dir, showWarnings = FALSE, recursive = TRUE)
cat("Output directory:", this_dir, "\n")

# ============================ CONFIG =========================================
M            <- 100        # Monte Carlo replications
fam          <- "NO"       # PL-NO throughout the validation
zeta         <- 2
lambda_true  <- 1.2
nu           <- 0.01

n_pl         <- 500        # PL validation sample size  (low-dim, MLE feasible)
mstop_pl     <- 800

n_ipl        <- 800        # IPL validation sample size.
                           #   800 so the
                           #   LOW-inflation scenario (alpha in [0.047,0.142]) has
                           #   enough boundary observations to fit the alpha submodel.
                           #   Lower it here if you want, but expect the low scenario
                           #   to fail in many reps below ~500.
mstop_ipl    <- 1000
c_infl       <- 1          # inflation at the upper boundary (mass at 1)
p_alpha      <- 2

set.seed(123)
# =============================================================================

fmt  <- function(m, s) if (is.na(m)) "---" else sprintf("$%.3f\\,(%.3f)$", m, s)
fmt1 <- function(x)   if (is.na(x)) "---" else sprintf("$%.3f$", x)
# A covariate that boosting never selected is absent from coef(): its estimate is
# ZERO (shrunk out), not missing. Returning NA here would average only over the
# replications in which it happened to be selected and overstate the estimate.
pick <- function(v, k, default = 0) if (k %in% names(v)) as.numeric(v[[k]]) else default
label_of <- function(k) {
  if (grepl("Intercept", k)) return("Intercept")
  if (grepl("^V[0-9]+$", k)) return(sprintf("$x_%s$", sub("^V", "", k)))
  if (grepl("^U[0-9]+$", k)) return(sprintf("$u_%s$", sub("^U", "", k)))
  k
}
normalise_intercept <- function(x) {
  names(x) <- gsub("^.*[Ii]ntercept.*$", "(Intercept)", names(x)); x
}
agg <- function(mat) list(mean = colMeans(mat, na.rm = TRUE),
                          sd   = apply(mat, 2, sd, na.rm = TRUE),
                          nok  = colSums(!is.na(mat)),
                          # proportion of replications with a non-zero estimate
                          # (for boosting this is the selection frequency)
                          selfreq = colMeans(abs(mat) > 1e-10, na.rm = TRUE))

# ======================= PL VALIDATION =======================================
mu_keys  <- c("(Intercept)", "V1", "V2", "V3", "V4")
sig_keys <- c("(Intercept)", "V3", "V4", "V5")

fit_pl_rep <- function(bmu, bsig) {
  X <- rnorm_multi(n_pl, vars = 5, mu = rep(0, 5), sd = rep(1, 5), r = 0.5,
                   varnames = paste0("V", 1:5), empirical = FALSE)
  Xm <- as.matrix(X)
  eta_mu  <- bmu["Intercept"]  + Xm %*% bmu[paste0("V", 1:5)]
  eta_sig <- bsig["Intercept"] + Xm %*% bsig[paste0("V", 1:5)]
  mu  <- pmax(pmin(plogis(eta_mu), 1 - 1e-6), 1e-6)
  sig <- exp(eta_sig)
  y   <- PLreg::rPL(n = n_pl, mu = mu, sigma = sig, lambda = lambda_true,
                    family = fam, zeta = zeta)
  y   <- pmax(pmin(y, 1 - 1e-6), 1e-6)
  dat <- data.frame(y = y, X)

  mb <- glmboostLSS(list(mu = y ~ V1 + V2 + V3 + V4, sigma = y ~ V3 + V4 + V5,
                         lambda = y ~ 1),
                    data = dat, families = PowerLogit(family = fam,
                    stabilization = "MAD", zeta = zeta),
                    control = boost_control(mstop = mstop_pl, nu = nu, trace = FALSE))
  b_mu  <- normalise_intercept(coef(mb, parameter = "mu",     off2int = TRUE))
  b_sig <- normalise_intercept(coef(mb, parameter = "sigma",  off2int = TRUE))
  b_lam <- as.numeric(unlist(coef(mb, parameter = "lambda",   off2int = TRUE)))

  # ATENCAO: o argumento de skewness NAO e' passado de proposito, para que o
  # PLreg estime o parametro em vez de receber o valor verdadeiro pronto.
  # Antes ele era informado, e por isso a coluna da MLE saia igual ao valor
  # verdadeiro, com desvio padrao zero: nao era estimativa, era a resposta
  # copiada, e a comparacao com o boosting nao dizia nada. O IPL da secao 2.5
  # sempre estimou, entao as duas validacoes estavam inconsistentes.
  mr <- PLreg::PLreg(y ~ V1 + V2 + V3 + V4 | V3 + V4 + V5, data = dat,
                     family = fam, zeta = zeta)
  ac <- coef(mr)
  m_mu  <- ac[!grepl("^\\(sigma\\)", names(ac))]
  m_sig <- setNames(ac[grepl("^\\(sigma\\)", names(ac))],
                    gsub("^\\(sigma\\)_", "", names(ac)[grepl("^\\(sigma\\)", names(ac))]))
  list(
    boost = c(setNames(sapply(mu_keys,  function(k) pick(b_mu,  k)), paste0("mu.",  mu_keys)),
              setNames(sapply(sig_keys, function(k) pick(b_sig, k)), paste0("sig.", sig_keys)),
              lambda = b_lam),
    mle   = c(setNames(sapply(mu_keys,  function(k) pick(m_mu,  k, NA_real_)), paste0("mu.",  mu_keys)),
              setNames(sapply(sig_keys, function(k) pick(m_sig, k, NA_real_)), paste0("sig.", sig_keys)),
              lambda = as.numeric(log(mr$coefficients$skewness)))
  )
}

run_pl_scenario <- function(bmu, bsig, tag) {
  cat(sprintf(">>> PL %s : %d reps\n", tag, M))
  B <- M2 <- NULL; ok <- 0; pml_fallback <- 0
  for (m in seq_len(M)) {
    # PLreg tenta verossimilhanca penalizada e, se nao converge, cai para ML.
    # Contamos quantas vezes isso acontece para poder relatar no texto.
    r <- withCallingHandlers(
           tryCatch(fit_pl_rep(bmu, bsig), error = function(e) NULL),
           warning = function(w) {
             if (grepl("pML", conditionMessage(w))) pml_fallback <<- pml_fallback + 1
             invokeRestart("muffleWarning")
           })
    if (is.null(r)) { next }
    B  <- rbind(B,  r$boost); M2 <- rbind(M2, r$mle); ok <- ok + 1
    if (m %% 20 == 0) cat("   ", m, "\n")
  }
  cat(sprintf("   %s: %d/%d reps ok | PLreg caiu de pML para ML em %d\n",
              tag, ok, M, pml_fallback))
  list(boost = agg(B), mle = agg(M2), ok = ok, pml_fallback = pml_fallback)
}

beta_mu_1  <- c(Intercept = 0.5, V1 = -0.99, V2 = 0.49, V3 = -1.17, V4 = 0.38, V5 = 0)
beta_sig_1 <- c(Intercept = 0.0, V1 = 0, V2 = 0, V3 = -0.41, V4 = -0.24, V5 = 0.81)
beta_mu_2  <- c(Intercept = 0.5, V1 = -0.50, V2 = 1.48, V3 = -1.45, V4 = 0.53, V5 = 0)
beta_sig_2 <- c(Intercept = 0.0, V1 = 0, V2 = 0, V3 = -0.69, V4 = -1.13, V5 = 0.91)
true_mu_1  <- c("(Intercept)"=0.5, V1=-0.99, V2=0.49, V3=-1.17, V4=0.38)
true_sig_1 <- c("(Intercept)"=0.0, V3=-0.41, V4=-0.24, V5=0.81)
true_mu_2  <- c("(Intercept)"=0.5, V1=-0.50, V2=1.48, V3=-1.45, V4=0.53)
true_sig_2 <- c("(Intercept)"=0.0, V3=-0.69, V4=-1.13, V5=0.91)

pl1 <- run_pl_scenario(beta_mu_1, beta_sig_1, "Scenario 1")
pl2 <- run_pl_scenario(beta_mu_2, beta_sig_2, "Scenario 2")

# ---- write val_pl.tex -------------------------------------------------------
cell <- function(sc, key) {
  if (!key %in% names(sc$mean)) return("---")
  fmt(sc$mean[[key]], sc$sd[[key]])
}
pl_row <- function(lead, key_disp, key_data, s1, s2)
  sprintf(" %s & %s & %s & %s & %s & & %s & %s & %s \\\\",
          lead, label_of(key_disp),
          fmt1(s1$true[[key_disp]]), cell(pl1$boost, key_data), cell(pl1$mle, key_data),
          fmt1(s2$true[[key_disp]]), cell(pl2$boost, key_data), cell(pl2$mle, key_data))

rows <- c()
for (i in seq_along(mu_keys)) {
  k <- mu_keys[i]; lead <- if (i==1) "$\\eta_\\mu$" else ""
  rows <- c(rows, pl_row(lead, k, paste0("mu.", k),
            list(true=true_mu_1), list(true=true_mu_2)))
}
rows <- c(rows, "\\midrule")
for (i in seq_along(sig_keys)) {
  k <- sig_keys[i]; lead <- if (i==1) "$\\eta_\\sigma$" else ""
  rows <- c(rows, pl_row(lead, k, paste0("sig.", k),
            list(true=true_sig_1), list(true=true_sig_2)))
}
rows <- c(rows, "\\midrule",
  sprintf("$\\log\\lambda$ & & %s & %s & %s & & %s & %s & %s \\\\",
          fmt1(log(lambda_true)), cell(pl1$boost,"lambda"), cell(pl1$mle,"lambda"),
          fmt1(log(lambda_true)), cell(pl2$boost,"lambda"), cell(pl2$mle,"lambda")))

writeLines(c(
  "% generated by validation_multirep.R -- do not edit by hand",
  "{\\footnotesize",
  "\\begin{tabular}{llrrrcrrr}", "\\toprule",
  "& & \\multicolumn{3}{c}{Scenario 1} & & \\multicolumn{3}{c}{Scenario 2} \\\\",
  "\\cmidrule(lr){3-5}\\cmidrule(lr){7-9}",
  "Submodel & Coefficient & True & Boosting & MLE & & True & Boosting & MLE \\\\",
  "\\midrule", rows, "\\bottomrule", "\\end{tabular}", "}"),
  file.path(this_dir, "val_pl.tex"))
cat("wrote val_pl.tex\n")

# ======================= IPL VALIDATION ======================================
family_ipl <- IPL(family = fam, stabilization = "MAD", c = c_infl, zeta = zeta)

IPLreg <- function(formula, data, family, zeta = 2, type = "ML", c) {
  of <- Formula(formula); mf <- model.frame(of, data = data)
  y <- model.response(mf); ind <- ifelse(y == c, 1, 0)
  dPL <- data[ind == 0, ]; dgl <- cbind(data, ind)
  fPL <- formula(of, lhs = 1, rhs = -3)
  fgl <- update(formula(of, lhs = 1, rhs = 3), ind ~ .)
  list(fitPL  = PLreg(fPL, data = dPL, family = family, zeta = zeta, type = type),
       fitglm = glm(fgl, family = binomial(), data = dgl))
}

al_keys <- c("(Intercept)", "U1", "U2")
beta_mu_ipl  <- c(-0.99, 0.49, -1.17, 0.38, 0)   # V1..V5
beta_sig_ipl <- c(0, 0, -0.41, -0.24, 0.81)
true_mu_ipl  <- c("(Intercept)"=0.5, V1=-0.99, V2=0.49, V3=-1.17, V4=0.38)
true_sig_ipl <- c("(Intercept)"=0.0, V3=-0.41, V4=-0.24, V5=0.81)

fit_ipl_rep <- function(kappa) {   # kappa = c(intercept, U1coef, U2coef)
  X <- rnorm_multi(n_ipl, vars = 5, mu = rep(0, 5), sd = rep(1, 5), r = 0.5,
                   varnames = paste0("V", 1:5), empirical = FALSE)
  Xm <- as.matrix(X)
  Xa <- matrix(runif(n_ipl * p_alpha), nrow = n_ipl); colnames(Xa) <- paste0("U", 1:p_alpha)
  eta_mu  <- 0.5 + Xm %*% beta_mu_ipl
  eta_sig <- 0.0 + Xm %*% beta_sig_ipl
  mu  <- pmax(pmin(plogis(eta_mu), 1 - 1e-6), 1e-6); sig <- exp(eta_sig)
  eta_a <- kappa[1] + Xa %*% kappa[-1]; alpha <- plogis(eta_a)
  y <- rIPL(n = n_ipl, mu = mu, sigma = sig, lambda = lambda_true,
            alpha = alpha, family = fam, zeta = zeta, c = c_infl)
  dat <- data.frame(y = y, Xa, X)

  mb <- glmboostLSS(list(mu = y ~ V1 + V2 + V3 + V4, sigma = y ~ V3 + V4 + V5,
                         lambda = y ~ 1, alpha = y ~ U1 + U2),
                    data = dat, families = family_ipl,
                    control = boost_control(mstop = mstop_ipl, nu = nu, trace = FALSE))
  b_mu  <- normalise_intercept(coef(mb, parameter = "mu",     off2int = TRUE))
  b_sig <- normalise_intercept(coef(mb, parameter = "sigma",  off2int = TRUE))
  b_lam <- as.numeric(unlist(coef(mb, parameter = "lambda",   off2int = TRUE)))
  b_al  <- normalise_intercept(coef(mb, parameter = "alpha",  off2int = TRUE))

  fit <- IPLreg(y ~ V1 + V2 + V3 + V4 | V3 + V4 + V5 | U1 + U2, data = dat,
                family = fam, zeta = zeta, type = "ML", c = c_infl)
  ac  <- coef(fit$fitPL)
  m_mu  <- ac[!grepl("^\\(sigma\\)", names(ac))]
  m_sig <- setNames(ac[grepl("^\\(sigma\\)", names(ac))],
                    gsub("^\\(sigma\\)_", "", names(ac)[grepl("^\\(sigma\\)", names(ac))]))
  # as.numeric() e' essencial: skewness vem nomeado "(lambda)", e sem isso o
  # c(lambda = m_lam) cria a chave "lambda.(lambda)" e a coluna some da tabela.
  m_lam <- as.numeric(log(fit$fitPL$coefficients$skewness))
  m_al  <- normalise_intercept(coef(fit$fitglm))
  list(
    boost = c(setNames(sapply(mu_keys,  function(k) pick(b_mu,  k)), paste0("mu.",  mu_keys)),
              setNames(sapply(sig_keys, function(k) pick(b_sig, k)), paste0("sig.", sig_keys)),
              lambda = b_lam,
              setNames(sapply(al_keys,  function(k) pick(b_al,  k)), paste0("al.",  al_keys))),
    mle   = c(setNames(sapply(mu_keys,  function(k) pick(m_mu,  k, NA_real_)), paste0("mu.",  mu_keys)),
              setNames(sapply(sig_keys, function(k) pick(m_sig, k, NA_real_)), paste0("sig.", sig_keys)),
              lambda = m_lam,
              setNames(sapply(al_keys,  function(k) pick(m_al,  k, NA_real_)), paste0("al.",  al_keys)))
  )
}

run_ipl_scenario <- function(kappa, tag) {
  cat(sprintf(">>> IPL %s : %d reps\n", tag, M))
  B <- M2 <- NULL; ok <- 0
  for (m in seq_len(M)) {
    r <- tryCatch(fit_ipl_rep(kappa), error = function(e) NULL)
    if (is.null(r)) next
    B <- rbind(B, r$boost); M2 <- rbind(M2, r$mle); ok <- ok + 1
    if (m %% 10 == 0) cat("   ", m, "\n")
  }
  cat(sprintf("   %s: %d/%d reps ok\n", tag, ok, M))
  list(boost = agg(B), mle = agg(M2), ok = ok)
}

kap_high <- c(-0.8, 0.4, -0.4)
kap_med  <- c(-1.5, 0.70, -0.3)
kap_low  <- c(-2.5, 0.7, -0.5)
true_al_high <- c("(Intercept)"=-0.8, U1=0.4, U2=-0.4)
true_al_med  <- c("(Intercept)"=-1.5, U1=0.70, U2=-0.3)
true_al_low  <- c("(Intercept)"=-2.5, U1=0.7, U2=-0.5)

ih <- run_ipl_scenario(kap_high, "High inflation")
im <- run_ipl_scenario(kap_med,  "Medium inflation")
il <- run_ipl_scenario(kap_low,  "Low inflation")

# ---- val_ipl_continuous.tex (mu, sigma, lambda; true common) ----------------
c3 <- function(sc, key) {
  if (!key %in% names(sc$mean)) return("---")
  fmt(sc$mean[[key]], sc$sd[[key]])
}
cont_row <- function(lead, kd, kk, tv)
  sprintf(" %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
          lead, label_of(kd), fmt1(tv),
          c3(ih$boost,kk), c3(ih$mle,kk), c3(im$boost,kk), c3(im$mle,kk),
          c3(il$boost,kk), c3(il$mle,kk))
rows <- c()
for (i in seq_along(mu_keys)) { k<-mu_keys[i]; rows<-c(rows,
  cont_row(if(i==1)"$\\eta_\\mu$" else "", k, paste0("mu.",k), true_mu_ipl[[k]])) }
rows <- c(rows, "\\midrule")
for (i in seq_along(sig_keys)) { k<-sig_keys[i]; rows<-c(rows,
  cont_row(if(i==1)"$\\eta_\\sigma$" else "", k, paste0("sig.",k), true_sig_ipl[[k]])) }
rows <- c(rows, "\\midrule",
  sprintf("$\\log\\lambda$ & & %s & %s & %s & %s & %s & %s & %s \\\\",
          fmt1(log(lambda_true)),
          c3(ih$boost,"lambda"), c3(ih$mle,"lambda"),
          c3(im$boost,"lambda"), c3(im$mle,"lambda"),
          c3(il$boost,"lambda"), c3(il$mle,"lambda")))
writeLines(c(
  "% generated by validation_multirep.R -- do not edit by hand",
  "{\\footnotesize",
  "\\begin{tabular}{llrrrrrrr}", "\\toprule",
  "& & & \\multicolumn{2}{c}{High} & \\multicolumn{2}{c}{Medium} & \\multicolumn{2}{c}{Low} \\\\",
  "\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\\cmidrule(lr){8-9}",
  "Submodel & Coefficient & True & Boost. & MLE & Boost. & MLE & Boost. & MLE \\\\",
  "\\midrule", rows, "\\bottomrule", "\\end{tabular}", "}"),
  file.path(this_dir, "val_ipl_continuous.tex"))
cat("wrote val_ipl_continuous.tex\n")

# ---- val_ipl_alpha.tex (row block per scenario) -----------------------------
al_block <- function(label, tv, sc) {
  out <- c()
  for (i in seq_along(al_keys)) {
    k <- al_keys[i]; lead <- if (i==1) label else ""
    out <- c(out, sprintf(" %s & %s & %s & %s & %s \\\\",
             lead, label_of(k), fmt1(tv[[k]]),
             c3(sc$boost, paste0("al.",k)), c3(sc$mle, paste0("al.",k))))
  }
  out
}
writeLines(c(
  "% generated by validation_multirep.R -- do not edit by hand",
  "{\\footnotesize",
  "\\begin{tabular}{llrrr}", "\\toprule",
  "Scenario & Coefficient & True & Boosting & MLE \\\\", "\\midrule",
  al_block("High inflation",   true_al_high, ih), "\\midrule",
  al_block("Medium inflation", true_al_med,  im), "\\midrule",
  al_block("Low inflation",    true_al_low,  il), "\\bottomrule",
  "\\end{tabular}", "}"),
  file.path(this_dir, "val_ipl_alpha.tex"))
cat("wrote val_ipl_alpha.tex\n")

# ---- machine-readable dump of every replication summary ---------------------
dump_rows <- function(tag, sc, truth) {
  ks <- names(sc$boost$mean)
  data.frame(
    scenario   = tag,
    coefficient= ks,
    true       = as.numeric(truth[ks]),
    boost_mean = as.numeric(sc$boost$mean[ks]),
    boost_sd   = as.numeric(sc$boost$sd[ks]),
    boost_selfreq = as.numeric(sc$boost$selfreq[ks]),
    mle_mean   = as.numeric(sc$mle$mean[ks]),
    mle_sd     = as.numeric(sc$mle$sd[ks]),
    reps_ok    = sc$ok,
    row.names  = NULL, stringsAsFactors = FALSE)
}
tv_pl <- function(tmu, tsig) {
  v <- c(setNames(as.numeric(tmu[mu_keys]),   paste0("mu.",  mu_keys)),
         setNames(as.numeric(tsig[sig_keys]), paste0("sig.", sig_keys)),
         lambda = log(lambda_true)); v
}
tv_ipl <- function(tal) {
  c(setNames(as.numeric(true_mu_ipl[mu_keys]),   paste0("mu.",  mu_keys)),
    setNames(as.numeric(true_sig_ipl[sig_keys]), paste0("sig.", sig_keys)),
    lambda = log(lambda_true),
    setNames(as.numeric(tal[al_keys]), paste0("al.", al_keys)))
}
allrows <- rbind(
  dump_rows("PL Scenario 1",        pl1, tv_pl(true_mu_1, true_sig_1)),
  dump_rows("PL Scenario 2",        pl2, tv_pl(true_mu_2, true_sig_2)),
  dump_rows("IPL High inflation",   ih,  tv_ipl(true_al_high)),
  dump_rows("IPL Medium inflation", im,  tv_ipl(true_al_med)),
  dump_rows("IPL Low inflation",    il,  tv_ipl(true_al_low)))
write.csv(allrows, file.path(this_dir, "validation_results.csv"), row.names = FALSE)
cat("wrote validation_results.csv\n")

# ---- human-readable markdown summary ----------------------------------------
md <- c(
  "# Validação multi-réplica — boosting vs máxima verossimilhança",
  "",
  sprintf("- Gerado em: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("- Réplicas Monte Carlo: M = %d", M),
  sprintf("- PL  (§2.4): n = %d, mstop = %d, nu = %.3f, família PL-%s", n_pl, mstop_pl, nu, fam),
  sprintf("- IPL (§2.5): n = %d, mstop = %d, nu = %.3f, massa em c = %d", n_ipl, mstop_ipl, nu, c_infl),
  sprintf("- λ verdadeiro = %.2f (log λ = %.3f); λ é ESTIMADO tanto no boosting quanto na MLE",
          lambda_true, log(lambda_true)),
  "",
  "Células: **média (desvio padrão)** entre as réplicas bem-sucedidas.",
  "`sel` = proporção de réplicas em que o boosting selecionou a covariável",
  "(covariável não selecionada entra como estimativa zero, não como faltante).",
  "")
fm <- function(x) if (is.na(x)) "---" else sprintf("%.3f", x)
for (tag in unique(allrows$scenario)) {
  d <- allrows[allrows$scenario == tag, ]
  md <- c(md, sprintf("## %s  (%d/%d réplicas ok)", tag, d$reps_ok[1], M), "",
          "| Coeficiente | Verdadeiro | Boosting (dp) | sel | MLE (dp) |",
          "|---|---|---|---|---|")
  for (i in seq_len(nrow(d))) {
    md <- c(md, sprintf("| `%s` | %s | %s (%s) | %.2f | %s (%s) |",
            d$coefficient[i], fm(d$true[i]),
            fm(d$boost_mean[i]), fm(d$boost_sd[i]), d$boost_selfreq[i],
            fm(d$mle_mean[i]),   fm(d$mle_sd[i])))
  }
  md <- c(md, "")
}
writeLines(md, file.path(this_dir, "validation_results.md"))
cat("wrote validation_results.md\n")

cat(sprintf("\n>>> done. M=%d | PL n=%d mstop=%d | IPL n=%d mstop=%d\n",
            M, n_pl, mstop_pl, n_ipl, mstop_ipl))
cat(sprintf("    reps ok: PL(%d,%d)  IPL(H=%d,M=%d,L=%d)\n",
            pl1$ok, pl2$ok, ih$ok, im$ok, il$ok))
