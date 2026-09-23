# ==============================================================================
# build_df_inep.r
# Builds the modeling data.table for the INEP application, now on the
# Educacao Superior CPC 2023 results (Conceito Enade continuo), scoped to
# MEDICINA / SP -- a single homogeneous course x state population.
#
# Motivation: the old Censo Escolar / ES version had n=1959 vs p=23 (3p=69),
# nowhere near the high-dimensional regime the method targets. This slice gives
# n=55 with p~19 predictors and 3 distributional parameters, so 3p=57 > n=55:
# squarely inside 3p > n, and still 3p > n_train (~44) inside each CV fold.
#
# Response: Conceito Enade (Continuo), a standardized score on [0,5].
#           tx := Conceito Enade / 5  -> unit interval. On MEDICINA/SP the score
#           never touches the boundaries (raw range [0.022, 4.970]), so tx is
#           strictly inside (0,1) and PL / Beta fit directly -- no boundary
#           adjustment needed. (A Smithson-Verkuilen squeeze would only be
#           required for a wider slice where the continuo hits exact 0 / 5.)
#
# Predictors (all non-circular w.r.t. Conceito Enade -- we deliberately DROP
# CPC and its own components as a response, and never use CPC continuo here):
#   counts        : Concluintes Inscritos, Concluintes Participantes
#   intake quality: Proporcao de concluintes com nota no Enem
#   corpo docente : Nota Padronizada Mestres / Doutores / Regime de Trabalho
#   estrutura     : Nota Padronizada Infraestrutura / Org. Didatico-Pedagogica /
#                   Oportunidade de Ampliacao
#   categoricals  : Organizacao Academica, Categoria Administrativa (one-hot).
#   Modalidade is constant (all Presencial in MEDICINA/SP) -> dropped.
# ==============================================================================

library(data.table)
library(readxl)
# --- resolucao de caminho relativa a este arquivo (nao depende de outras pastas) --
.script_dir <- function() {
  for (i in seq_len(sys.nframe())) {
    o <- sys.frame(i)$ofile
    if (!is.null(o)) return(dirname(normalizePath(o)))
  }
  a  <- commandArgs(trailingOnly = FALSE)
  sp <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(sp)) return(dirname(normalizePath(sp)))
  getwd()
}


app_path  <- Sys.getenv("PLBOOST_OUT", unset = ".")   # raw INEP files in, df_medicina_SP.RData out, both under data/
data_path <- file.path(app_path, "data")
dir.create(data_path, showWarnings = FALSE, recursive = TRUE)

AREA_FILTER <- "MEDICINA"
UF_FILTER   <- "SP"

# ── READ CPC 2023 results (one row per course) ────────────────────────────────
raw <- as.data.table(read_xlsx(file.path(data_path, "CPC_2023.xlsx"), sheet = 1))
# Normalise accented headers to ASCII so column selection is encoding-proof.
setnames(raw, make.unique(iconv(colnames(raw), to = "ASCII//TRANSLIT")))
g <- function(pat) grep(pat, colnames(raw), value = TRUE)[1]

col_area <- g("Area de Aval")
col_uf   <- g("Sigla da UF")

df <- raw[get(col_area) == AREA_FILTER & get(col_uf) == UF_FILTER]
cat("MEDICINA/SP rows:", nrow(df), "\n")

# ── RESPONSE: Conceito Enade continuo -> unit interval -> SV squeeze ───────────
col_resp <- g("Conceito Enade .Cont")
enade    <- suppressWarnings(as.numeric(df[[col_resp]]))
df <- df[!is.na(enade)]
enade <- enade[!is.na(enade)]

tx <- enade / 5                                      # [0,5] -> (0,1); strictly interior here
df[, tx := tx]

# ── JOIN Censo da Educacao Superior 2023 (course-level microdata) ─────────────
# Adds intake/financing/competition covariates so the regressed budget beats n:
# 19 CPC predictors + 14 census predictors = p = 33 -> 2p + 1 = 67 > n = 55.
# Join key: CO_CURSO (== "Codigo do Curso" in the CPC file). All non-circular
# w.r.t. Conceito Enade; counts are size, shares are composition.
df[, CO_CURSO := as.integer(get(g("odigo do Curso")))]

censo_csv <- file.path(data_path, "censo_superior_2023",
                       "microdados_censo_da_educacao_superior_2023", "dados",
                       "MICRODADOS_CADASTRO_CURSOS_2023.CSV")
censo_sel <- c("CO_CURSO", "QT_VG_TOTAL", "QT_ING", "QT_MAT", "QT_CONC",
               "QT_ING_FEM", "QT_ING_ENEM", "QT_ING_18_24", "QT_ING_BRANCA",
               "QT_ING_PRETA", "QT_ING_PARDA", "QT_CONC_FIES", "QT_CONC_PROUNII",
               "QT_CONC_PROUNIP", "QT_CONC_RESERVA_VAGA", "QT_CONC_PROCESCPUBLICA")
censo <- fread(censo_csv, sep = ";", encoding = "Latin-1", select = censo_sel)
censo[, CO_CURSO := as.integer(CO_CURSO)]

df <- merge(df, censo, by = "CO_CURSO", all.x = TRUE)
cat("Matched in census:", sum(!is.na(df$QT_ING)), "of", nrow(df), "\n")

# Composition shares (normalised so they carry mix, not size) + structural ratios.
sh <- function(a, b) round(a / pmax(b, 1), 4)
df[, `:=`(
  sh_ing_fem        = sh(QT_ING_FEM,   QT_ING),
  sh_ing_enem       = sh(QT_ING_ENEM,  QT_ING),
  sh_ing_jovem      = sh(QT_ING_18_24, QT_ING),
  sh_ing_branca     = sh(QT_ING_BRANCA, QT_ING),
  sh_ing_pretaparda = sh(QT_ING_PRETA + QT_ING_PARDA, QT_ING),
  sh_conc_fies      = sh(QT_CONC_FIES, QT_CONC),
  sh_conc_prouni    = sh(QT_CONC_PROUNII + QT_CONC_PROUNIP, QT_CONC),
  sh_conc_cotas     = sh(QT_CONC_RESERVA_VAGA, QT_CONC),
  sh_conc_escpub    = sh(QT_CONC_PROCESCPUBLICA, QT_CONC),
  concorrencia      = sh(QT_ING, QT_VG_TOTAL),      # ingressantes / vagas
  taxa_concl        = sh(QT_CONC, QT_ING)           # concluintes / ingressantes
)]
censo_num <- c("QT_VG_TOTAL", "QT_ING", "QT_MAT",
               "sh_ing_fem", "sh_ing_enem", "sh_ing_jovem", "sh_ing_branca",
               "sh_ing_pretaparda", "sh_conc_fies", "sh_conc_prouni",
               "sh_conc_cotas", "sh_conc_escpub", "concorrencia", "taxa_concl")

# ── PREDICTORS ────────────────────────────────────────────────────────────────
num_patterns <- c(
  "Concluintes Inscritos",
  "Concluintes Participantes",
  "Propor.*Enem",
  "Padronizada - Mestres",
  "Padronizada - Doutores",
  "Padronizada - Regime",
  "Padronizada - Infra",
  "Padronizada - Organiza",          # Nota Padronizada Org. Didatico-Pedagogica
  "Padronizada - Oportun"
)
num_cols <- c(vapply(num_patterns, g, character(1)), censo_num)
stopifnot(!any(is.na(num_cols)))
for (cc in num_cols) set(df, j = cc, value = suppressWarnings(as.numeric(df[[cc]])))

col_org <- g("Organiza.*Acad")       # Organizacao Academica
col_cat <- g("Categoria Admin")

# Clean, accent-free level labels. macOS iconv//TRANSLIT is unreliable and can
# leave R unicode escapes (e.g. "U00FA" for u-acute), so map those explicitly.
clean_lab <- function(x) {
  x <- as.character(x)
  esc <- c(U00E1="a", U00E2="a", U00E3="a", U00E9="e", U00EA="e", U00ED="i",
           U00F3="o", U00F4="o", U00F5="o", U00FA="u", U00E7="c",
           U00C1="A", U00C9="E", U00CD="I", U00D3="O", U00DA="U", U00C7="C")
  x <- tools::toTitleCase(x)
  x <- gsub("[^A-Za-z0-9]", "", x)          # "<U+00E1>" -> "U00E1"
  for (k in names(esc)) x <- gsub(k, esc[[k]], x, fixed = TRUE)
  x
}

# One-hot both categoricals (keep all levels; boosting simply won't select
# the dead columns, and full one-hot is what pushes p up to hit 3p > n).
df[, (col_org) := factor(get(col_org))]
df[, (col_cat) := factor(get(col_cat))]
dum_org <- model.matrix(~ get(col_org) - 1, data = df)
dum_cat <- model.matrix(~ get(col_cat) - 1, data = df)
colnames(dum_org) <- paste0("org_", clean_lab(gsub("get\\(col_org\\)", "", colnames(dum_org))))
colnames(dum_cat) <- paste0("cat_", clean_lab(gsub("get\\(col_cat\\)", "", colnames(dum_cat))))

df_model <- cbind(
  df[, c("tx", num_cols), with = FALSE],
  as.data.table(dum_org),
  as.data.table(dum_cat)
)

# Drop any constant columns (e.g. a Categoria level with a single course).
const_cols <- names(df_model)[sapply(df_model, function(x) length(unique(x)) <= 1)]
const_cols <- setdiff(const_cols, "tx")
if (length(const_cols) > 0) {
  cat("Dropping constant columns:", paste(const_cols, collapse = ", "), "\n")
  df_model[, (const_cols) := NULL]
}

# ── RENAME PREDICTORS TO V1..Vn + DICTIONARY ─────────────────────────────────
pred_cols <- setdiff(names(df_model), "tx")
new_names <- paste0("V", seq_along(pred_cols))
var_dict  <- data.table(new_name = new_names, orig_name = pred_cols)

setnames(df_model, old = pred_cols, new = new_names)
setcolorder(df_model, c("tx", new_names))
df_inep <- df_model

# ── REPORT ────────────────────────────────────────────────────────────────────
n_obs <- nrow(df_inep); p <- length(new_names)
cat("\nFinal df_inep:", n_obs, "rows x", p, "predictors\n")
# PL regresses mu and sigma on the p predictors; lambda ~ 1 (one scalar).
# Regression-coefficient budget = 2p + 1; nominal 3-parameter budget = 3p.
cat(sprintf("2p+1 = %d (regressed: mu,sigma + lambda scalar)   vs   n = %d   ->   %s\n",
            2 * p + 1, n_obs,
            if (2 * p + 1 > n_obs) "over-parameterised (need boosting)" else "2p+1 <= n"))
cat(sprintf("3p   = %d (nominal, 3 distributional params)      vs   n = %d\n", 3 * p, n_obs))
cat("tx range: [", round(min(df_inep$tx), 4), ",", round(max(df_inep$tx), 4), "]\n")
cat("tx mean :", round(mean(df_inep$tx), 4), "\n")
print(var_dict)

save(df_inep, var_dict, file = file.path(data_path, "df_medicina_SP.RData"))
cat("\nSaved:", file.path(data_path, "df_medicina_SP.RData"), "\n")
