# 03_Applications/itagrade_appl/

Application of the **PLboost** methodology to academic performance modeling — specifically, predicting the **GPA ratio**: each student's GPA divided by the maximum possible GPA, a continuous response strictly bounded in (0, 1).

This application relates to **SDG 4** (Quality Education), as understanding the environmental and demographic factors that drive academic performance supports equitable and effective education policy.

---

## Dataset

**Source:** [Zenodo — Student Academic Performance Dataset](https://zenodo.org/records/15423018)

The dataset contains academic and environmental data from **171 university students** across Milan and Rome (Mai, 2025).

### How to download

The script loads the data automatically — no manual download needed:

```r
url <- "https://zenodo.org/records/15423018/files/dataset.csv?download=1"
df  <- as.data.table(read.csv(url))
```

> The dataset is open access and licensed under CC0. It is small enough (~171 rows) that `df_itagrade.RData` — the pre-processed version — is included directly in this folder.

---

## Variables

### Response variable

| Variable | Description |
|---|---|
| `y` (GPA_ratio) | Student GPA divided by maximum possible GPA — bounded in (0, 1) |

A small epsilon correction (`1e-4`) is applied to avoid exact boundary values, though in this dataset `y` already lies in the range (0.667, 0.990).

### Predictors

**Demographic characteristics:**

| Variable | Description |
|---|---|
| `Gender` | Sex of the student |
| `Age` | Age in years (continuous) |
| `City` | City of study — Milano, Roma, Other |
| `Major` | Field of study — Economics, Engineering, Mathematics, Medicine, and others |

**Study environment conditions (binary):**

| Variable | Description |
|---|---|
| `study_time` | Period of study — Morning, Afternoon, Evening, Night |
| `study_location` | Place of study — Home, Library, Café, Other |
| `enough_sleep` | Whether the student sleeps enough |
| `natural_sun_exposure` | Whether the student has natural light exposure |
| `noisy_environment` | Whether the study environment is noisy |
| `heated_cooled` | Whether heating/cooling is adequate |
| `ventilated` | Whether the environment is well ventilated |
| `enough_desk_space` | Whether there is sufficient desk space |
| `often_distracted` | Whether the student is frequently distracted — **key predictor** (strong negative effect on GPA ratio) |
| `study_in_group` | Whether the student studies in groups |

---

## What the scripts do

### `itagrade_application.r`
Data preparation pipeline. Loads the raw CSV directly from Zenodo and:

1. **Builds the response** `y = GPA / max_GPA`, clipped to (ε, 1 − ε)
2. **Removes leakage** — drops the raw `GPA` and `GPA_ratio` columns after constructing `y`
3. **One-hot encodes** all categorical variables (City, Major, study_time, study_location, Gender)
4. **Saves** the cleaned dataset as `df_itagrade.RData`

### `itagrade_application_mod.r`
Model fitting and evaluation. Using the pre-processed dataset, it:

1. **Renames predictors** to `V1, V2, V3, ...` (blind fitting) and saves a variable dictionary
2. **Fits** a `glmboostLSS` model with the Power Logit Normal family, regressing μ and σ on all predictors and keeping λ fixed (`lambda ~ 1`)
3. **Tunes** the stopping iteration (`mstop`) via cross-validation on a grid up to 1500 iterations
4. **Evaluates** predictive performance across **100 train/test splits** (80/20), recording MSE per replicate
5. **Saves** all results including fitted coefficients and variable selection outcomes

---

## Files in this folder

| File | Description |
|---|---|
| `itagrade_application.r` | Data preparation pipeline |
| `itagrade_application_mod.r` | Model fitting and evaluation (100 replicates) |
| `df_itagrade.RData` | Pre-processed dataset — included, ready to use |
| `sim_results_itagrade.RData` | Fitted model results across all replicates |
| `essentials/` | Local copy of core PLboost functions |

---

## Dependencies

```r
install.packages(c(
  "data.table", "mltools",
  "PLreg", "mboost", "gamboostLSS", "MASS"
))
```
