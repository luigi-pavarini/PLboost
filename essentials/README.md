# essentials/

Core R functions shared across all modules of the **PLboost** project. Every script in `01_Testing/`, `02_Simulation/`, and `03_Applications/` depends on the files in this folder.

---

## Files

### `families_PL.R`
Implements the `gamboostLSS`-compatible boosting families for the **Power Logit (PL)** distribution. Defines three separate family objects — one per distributional parameter — that the boosting algorithm optimises independently:

| Family object | Parameter | Link function |
|---|---|---|
| `PLMu` | μ (location / median) | logit |
| `PLSigma` | σ (dispersion) | log |
| `PLLambda` | λ (skewness) | log — kept **fixed**, not regressed |

Supports all seven PL mixing families: `NO`, `TF`, `LO`, `PE`, `Hyp`, `SN`, `SLASH`.

**Dependencies:** `PLreg`, `mboost`, `gamboostLSS`, `MASS`

---

### `families_IPL.R`
Implements the `gamboostLSS`-compatible boosting families for the **Inflated Power Logit (IPL)** distribution — an extension of PL that accommodates probability mass at the boundaries (0 or 1). Defines four family objects:

| Family object | Parameter | Link function |
|---|---|---|
| `IPLMu` | μ (location / median) | logit |
| `IPLSigma` | σ (dispersion) | log |
| `IPLLambda` | λ (skewness) | log — kept **fixed** |
| `IPLAlpha` | α (inflation probability at boundary *c*) | logit |

The boundary `c` is set to either `0` or `1` at construction time.

**Dependencies:** `PLreg`, `mboost`, `gamboostLSS`, `MASS`, `VGAM`, `zipfR`, `gamlss.dist`, `GeneralizedHyperbolic`

---

### `IPL_distributions.R`
Provides the base distribution functions for the IPL model, following standard R conventions (`d`/`p`/`q`/`r`):

| Function | Description |
|---|---|
| `dIPL` | Density — returns `alpha` at the inflation point `c`, PL density elsewhere |
| `pIPL` | CDF |
| `qIPL` | Quantile function |
| `rIPL` | Random generation |

The IPL model is a mixture: with probability `alpha` the response equals `c` (point mass); with probability `1 - alpha` it follows a PL distribution.

**Dependencies:** `PLreg`, `VGAM`, `zipfR`

---

### `v_function.r`
Computes the **score auxiliary quantities** `z`, `v(z)`, and `v'(z)` used in the negative gradient calculation inside all boosting families. These are the key ingredients of the working scores for each distributional parameter.

| Argument | Description |
|---|---|
| `y` | Response vector in (0, 1) |
| `mu` | Fitted median vector in (0, 1) |
| `sigma` | Dispersion — scalar or vector |
| `lambda` | Power parameter — scalar |
| `family` | One of `"NO"`, `"TF"`, `"LO"`, `"SN"`, `"PE"`, `"Hyp"`, `"SLASH"` |
| `zeta` | Extra shape parameter (default = 2) |

Returns a named list: `z` (standardised latent residual), `v` (weight function), `dv` (derivative).

**Dependencies:** `VGAM`, `zipfR`

---

## Usage

Source all four files at the top of any script before fitting models:

```r
source(".../essentials/families_PL.R")
source(".../essentials/families_IPL.R")
source(".../essentials/IPL_distributions.R")
source(".../essentials/v_function.r")
```

> **Note:** `families_PL.R` and `families_IPL.R` internally source `v_function.r` and `IPL_distributions.R`. If you source the family files, the auxiliary files are loaded automatically — as long as all four files are in the same directory.

---

## Dependencies summary

```r
install.packages(c(
  "PLreg", "mboost", "gamboostLSS", "MASS",
  "VGAM", "zipfR", "gamlss.dist", "GeneralizedHyperbolic"
))
```
