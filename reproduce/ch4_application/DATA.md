# Data and scripts for the applications

Chapter 4 has two applications, each in its own subfolder:

- [`credit/`](credit) — **loss given default (LGD)** on Lending Club loans (inflated PL).
- [`inep/`](inep)   — **Conceito Enade** of Medicina courses in São Paulo (plain PL).

Every script resolves its working directory from the `PLBOOST_OUT` environment
variable (default: the current directory). That directory is expected to contain

```
$PLBOOST_OUT/
  essentials/   # families_PL.R, families_IPL.R, cv_load.R, cv_summary.R
  data/         # the prepared .RData for each application (see below)
  cross_val/    # created by the scripts for results and figures
```

The `essentials/cv_load.R`, `essentials/cv_summary.R` and `essentials/selection_table.R` helpers are shipped in
[`essentials/`](essentials); copy them (together with the family files) into
`$PLBOOST_OUT/essentials/`.

---

## Credit — Loss Given Default

The raw file is about 1.6 GB and is redistributed under Kaggle's terms, so it is
**not** included here. Download `accepted_2007_to_2018Q4.csv` from

<https://www.kaggle.com/datasets/wordsforthewise/lending-club>

(with the Kaggle CLI: `kaggle datasets download -d wordsforthewise/lending-club -f accepted_2007_to_2018Q4.csv --unzip`)
and build `data/df_lgd_all.RData`.

The analysed sample filters to defaulted loans in the state of Vermont, giving
`n = 320` loans and, after one-hot encoding, `p = 348` candidate covariates. The
response is the unrecovered fraction of the loan, restricted to `(0, 1]`; roughly
40% sits exactly at 1, the point mass the inflated model handles.

Predictive performance is assessed by **5-fold cross-validation repeated 10 times
(50 folds)**, stratified on the inflation indicator so every fold keeps the ~60/40 mix.

### Order of execution (`credit/`)

1. `credit_run.r`     — fits the six inflated power logit generators and the one-inflated beta baseline (BEOI) over the 50 folds
2. `credit_summary.r` — assembles the family-comparison table (MSE, MAE, Upsilon)
3. `credit_img_ipl.r` — coefficient, selection and diagnostic plots for the chosen IPL-t model (`fam <- "TF"`)
4. `credit_img_beoi.r` — the BEOI baseline plots
5. `../essentials/selection_table.R credit/credit_img_ipl.r` — the selection table (frequency, average coefficient, band avg ± 1.96 sd)

Supplementary, not in the text of dissertation v1.0 — diagnostics of a single fit on the
full sample, in `credit/final_fit/`: `credit_final_fit.R` (IPL-t final fit and quantile
residuals), then `credit_qres_envelope.R` (three residual panels with a simulated
envelope) and `credit_beoi_qres.R` (the same panels for BEOI).

---

## INEP — Conceito Enade (course quality)

The response is the Conceito Enade of the Medicina courses of São Paulo, rescaled to
the unit interval. It is built from two open INEP microdata sources of 2023, the
Conceito Preliminar de Curso (CPC) and the Censo da Educação Superior:

<https://www.gov.br/inep/pt-br/acesso-a-informacao/dados-abertos/microdados>

The cleaned sample has `n = 55` courses and `q = 33` predictors, with no boundary mass,
so the plain (non-inflated) power logit class is used. Same 50-fold cross-validation.

### Order of execution (`inep/`)

1. `extract_inep.r`, `build_df_inep.r` — download/parse the microdata into `data/df_medicina_SP.RData`
2. `inep_run.r`      — fits the six power logit generators and the beta baseline over the 50 folds
3. `inep_summary.r`, `inep_table.r` — family-comparison table and predictive summary
4. `inep_img.r`      — coefficient, selection and diagnostic plots for the PL-NO model
5. `../essentials/selection_table.R inep/inep_img.r` — the selection table

Supplementary, not in the text of dissertation v1.0: `inep/final_fit/inep_final_fit.R`,
the PL-NO final fit on the full sample with a quantile-residual envelope.
