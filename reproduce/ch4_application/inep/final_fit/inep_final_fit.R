# ==============================================================================
# inep_final_fit.R
# Final PL-NO fit on the full sample (Medicina, Sao Paulo) and quantile-residual
# diagnostics, mirroring credit_final_fit.R.
#
# There is no inflation here: the response tx = Conceito Enade / 5 is strictly
# interior, so the residual is the plain quantile residual, no randomization:
#     r_i^q = Phi^{-1}[ F(y_i; mu_i, sigma_i, lambda) ]
# and a single panel instead of the three of the inflated case.
#
# Envelope by refitting (Atkinson): minimum and maximum over B replicates
# simulated from the fitted model. B defaults to 200; override with B_ENV.
# Outputs in $PLBOOST_OUT/final_fit/: inep_qres_qq.png,
# inep_final_estimates.md, inep_final_fit.RData.
# Supplementary diagnostic: not part of the text of dissertation v1.0.
# ==============================================================================

library(data.table); library(mboost); library(gamboostLSS)
library(gamlss.dist); library(EnvStats); library(parallel)

app_path <- Sys.getenv("PLBOOST_OUT", unset = ".")   # must contain essentials/ and data/
out_path <- file.path(app_path, "final_fit")
dir.create(out_path, showWarnings = FALSE)
source(file.path(app_path, "essentials", "families_PL.R"))
source(file.path(app_path, "essentials", "IPL_distributions.R"))

load(file.path(app_path, "data", "df_medicina_SP.RData"))
pred_cols <- setdiff(names(df_inep), "tx")
var_dict  <- data.table(new_name = paste0("V", seq_along(pred_cols)),
                        orig_name = pred_cols)
setnames(df_inep, old = pred_cols, new = var_dict$new_name)
y <- df_inep$tx; n <- length(y)
cat("n =", n, " p =", length(pred_cols),
    " | tx em [", round(min(y),4), ",", round(max(y),4), "]\n")

FAM <- "NO"; ZETA <- 2                     # PL-NO nao tem parametro extra
mstop_init <- 500; m_stop_cv <- 800; nu <- 0.001
B <- as.integer(Sys.getenv("B_ENV", "200"))

grid_pl <- make.grid(max = c(mu = m_stop_cv, sigma = m_stop_cv, lambda = m_stop_cv),
                     min = mstop_init, length.out = 2, dense_mu_grid = TRUE)

# `fam` no ambiente global: o cvrisk() reavalia a call noutro ambiente
fit_pl <- function(dat, mstop_vec) {
  assign("fam", PowerLogit(family = FAM, stabilization = "MAD", zeta = ZETA),
         envir = globalenv())
  glmboostLSS(formula = list(mu = tx ~ ., sigma = tx ~ ., lambda = tx ~ 1),
              data = dat, families = fam,
              control = boost_control(mstop = mstop_vec, nu = nu, trace = FALSE))
}
qn <- function(p) qnorm(pmin(pmax(p, 1e-12), 1 - 1e-12))

upsilon_pl <- function(yy, mu, sg, lb) {
  fda <- sort(pPL(yy, mu = mu, sigma = sg, lambda = lb, zeta = ZETA, family = FAM))
  mean(abs(qn(fda) - EnvStats::evNormOrdStats(n = length(fda))))
}

cat("\n[1] mstop por cvrisk na base completa\n")
m0 <- fit_pl(df_inep, rep(mstop_init, 3))
set.seed(20260812L)
opt <- mstop(cvrisk(m0, grid = grid_pl, folds = mboost::cv(model.weights(m0)),
                    papply = mclapply))
print(opt)

cat("\n[2] ajuste final\n")
fit  <- fit_pl(df_inep, c(opt["mu"], opt["sigma"], opt["lambda"]))
mu_h <- as.numeric(predict(fit$mu,     newdata = df_inep, type = "response"))
sg_h <- as.numeric(predict(fit$sigma,  newdata = df_inep, type = "response"))
lb_h <- as.numeric(predict(fit$lambda, newdata = df_inep, type = "response"))[1]

rq <- pmin(pmax(qn(pPL(y, mu = mu_h, sigma = sg_h, lambda = lb_h,
                       zeta = ZETA, family = FAM)), -8), 8)
cat("   r^q: media", round(mean(rq),3), " sd", round(sd(rq),3),
    " min", round(min(rq),3), " max", round(max(rq),3), "\n")

cat("\n[3] envelope por reajuste, B =", B, "\n")
S <- matrix(NA_real_, n, B); t0 <- Sys.time(); ok <- 0
for (b in seq_len(B)) {
  ys <- rPL(n = n, mu = mu_h, sigma = sg_h, lambda = lb_h, zeta = ZETA, family = FAM)
  ys <- pmin(pmax(ys, 1e-6), 1 - 1e-6)
  d2 <- copy(df_inep); d2[, tx := ys]
  m2 <- try(fit_pl(d2, c(opt["mu"], opt["sigma"], opt["lambda"])), silent = TRUE)
  if (inherits(m2, "try-error")) next
  S[, b] <- sort(pmin(pmax(qn(pPL(ys,
      mu = as.numeric(predict(m2$mu,     newdata = d2, type = "response")),
      sigma = as.numeric(predict(m2$sigma,  newdata = d2, type = "response")),
      lambda = as.numeric(predict(m2$lambda, newdata = d2, type = "response"))[1],
      zeta = ZETA, family = FAM)), -8), 8))
  ok <- ok + 1; rm(m2, d2); gc()
  if (b %% 25 == 0) cat("   b", b, "|", round(difftime(Sys.time(), t0, units="min"),1), "min\n")
}
cat("   reajustes ok:", ok, "de", B, "\n")

lo <- apply(S, 1, min, na.rm = TRUE); hi <- apply(S, 1, max, na.rm = TRUE)
mi <- apply(S, 1, median, na.rm = TRUE); qt <- qnorm(ppoints(n))
png(file.path(out_path, "inep_qres_qq.png"), width = 1500, height = 1400, res = 240)
par(mar = c(4.2, 4.4, 1, 1))
plot(qt, sort(rq), pch = 20, cex = .8, xlab = "Theoretical quantiles",
     ylab = "Quantile residuals", ylim = range(c(rq, lo, hi), na.rm = TRUE))
lines(qt, lo); lines(qt, hi); lines(qt, mi, lty = 2)
dev.off()
fora <- sum(sort(rq) < lo | sort(rq) > hi)
cat("   fora do envelope:", fora, "de", n, "\n")

cf <- function(sub, nome) {
  x <- coef(fit[[sub]], off2int = TRUE); if (!length(x)) return(NULL)
  d <- data.table(Submodel = nome, Variable = names(x), Estimate = as.numeric(x))
  d[, Variable := ifelse(Variable %in% var_dict$new_name,
        var_dict$orig_name[match(Variable, var_dict$new_name)], Variable)]
  d[order(-abs(Estimate))]
}
tab <- rbindlist(list(cf("mu","eta_mu"), cf("sigma","eta_sigma")), use.names = TRUE)
writeLines(c("# Ajuste final PL-NO, base completa", "",
  sprintf("- n = %d, mstop = mu %s, sigma %s, lambda %s", n, opt["mu"], opt["sigma"], opt["lambda"]),
  sprintf("- lambda_hat = %.4f", lb_h),
  sprintf("- Upsilon = %.4f", upsilon_pl(y, mu_h, sg_h, lb_h)),
  sprintf("- r^q: media %.3f, sd %.3f, min %.3f, max %.3f", mean(rq), sd(rq), min(rq), max(rq)),
  sprintf("- fora do envelope: %d de %d", fora, n), "",
  "| Submodel | Variable | Estimate |", "|---|---|---|",
  sprintf("| %s | %s | %.4f |", tab$Submodel, tab$Variable, tab$Estimate)),
  file.path(out_path, "inep_final_estimates.md"))

save(fit, rq, S, mu_h, sg_h, lb_h, opt, var_dict, df_inep,
     file = file.path(out_path, "inep_final_fit.RData"))
cat("\nPRONTO em", round(difftime(Sys.time(), t0, units="min"),1), "min\n")
