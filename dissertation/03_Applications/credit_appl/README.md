# credit_appl

Application of the **PLboost** methodology to credit risk modeling — specifically, estimating **Loss Given Default (LGD)**: the fraction of a loan amount that is lost when a borrower defaults.

This application relates to **SDG 8** (Decent Work and Economic Growth), as accurate LGD estimation supports inclusive and resilient financial systems.

---

## Dataset

**Source:** [All Lending Club Loan Data — Kaggle](https://www.kaggle.com/datasets/wordsforthewise/lending-club)

The analysis uses the file **`accepted_2007_to_2018Q4.csv`**, which contains all accepted loan applications issued by Lending Club from 2007 through Q4 2018.

## What the scripts do

### `credit_application.r`
Full data preparation pipeline. Starting from the raw CSV, it:

1. **Filters** to defaulted loans only (`loan_status` ∈ `"Charged Off"`, `"Default"`)
2. **Constructs the LGD response variable** — the proportion of the loan amount not recovered after default: `LGD = 1 - recoveries / loan_amnt`, clipped to (0, 1)
3. **Engineers features** — credit age, days to last payment, co-borrower age
4. **Encodes categoricals** via one-hot encoding (grade, term, home ownership, purpose, state, etc.)
5. **Cleans and selects** the final predictor matrix, removing leakage variables (any information only available after the default event)

Key variables used:

| Variable | Description |
|---|---|
| `loan_amnt` | Total approved loan amount (USD) |
| `recoveries` | Amount recovered post charge-off — used to build LGD, then dropped |
| `int_rate` | Annual interest rate (%) |
| `dti` | Debt-to-income ratio |
| `fico_range_low/high` | FICO credit score range at origination |
| `grade` / `sub_grade` | Lending Club risk grade (A–G) |
| `annual_inc` | Borrower self-reported annual income |
| `revol_util` | Revolving credit line utilisation rate |

### `credit_application_mod.r` / `credit_application_mod_2.r`
Model fitting scripts. Applies the PL boosting families to the prepared LGD data, estimating location (μ), dispersion (σ), and skewness (λ) as functions of the selected covariates.

---

## Files in this folder

| File | Description |
|---|---|
| `credit_application.r` | Data preparation pipeline |
| `credit_application_mod.r` | Model fitting (version 1) |
| `credit_application_mod_2.r` | Model fitting (version 2) |
| `df_lgd_all.RData` | Pre-processed LGD dataset (ready for modelling) |
| `sim_results_credit_VT.RData` | Fitted model results |
| `accepted_2007_to_2018Q4.csv` | ⚠️ **Not included** — download from Kaggle (see above) |
| `essentials/` | Local copy of core PLboost functions |

---

## Dependencies

```r
install.packages(c(
  "data.table", "mltools",
  "PLreg", "mboost", "gamboostLSS", "MASS"
))
```
