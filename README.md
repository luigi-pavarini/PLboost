# PLboost

Gradient boosting families for the **power logit (PL)** class of distributions and its
adjusted extension (**APL**), for regression with a response bounded on the unit interval.

The package supplies `Families` objects compatible with
[**gamboostLSS**](https://cran.r-project.org/package=gamboostLSS), so the median,
dispersion, skewness and (in the adjusted case) boundary-probability submodels are
estimated jointly by component-wise gradient boosting. Variable selection and
regularization come out of the fitting procedure itself, which makes the approach usable
when the number of candidate covariates exceeds the number of observations.

It accompanies the master thesis *Boosting-based Inference and Variable Selection for
Modeling Bounded Continuous Data* (IME-USP). The power logit class is due to
Queiroz and Ferrari (2024), and the extension with a point mass at a boundary to
Queiroz and Ferrari (2024b), who call it the inflated power logit. Here it is called
the adjusted power logit, since 0 and 1 do not belong to the support of the PL
distribution and are incorporated by adjusting it.

## Installation

```r
# install.packages("remotes")
remotes::install_github("luigi-pavarini/PLboost")
```

`PLreg`, `gamboostLSS`, `mboost`, `VGAM` and `zipfR` are required and will be pulled in.

## Quick start

```r
library(PLboost)
library(gamboostLSS)

## three-parameter power logit, normal generator
fam <- PowerLogit(family = "NO", stabilization = "MAD")

fit <- glmboostLSS(
  formula  = list(mu = y ~ x1 + x2, sigma = y ~ x1 + x2, lambda = y ~ 1),
  data     = df,
  families = fam,
  control  = boost_control(mstop = 500, nu = 0.01)
)

coef(fit, parameter = "mu")
```

For a response with a point mass at a boundary:

```r
## four-parameter adjusted power logit, mass at one
fam <- APL(family = "TF", zeta = 3, c = 1, stabilization = "MAD")

fit <- glmboostLSS(
  formula  = list(mu = y ~ ., sigma = y ~ ., lambda = y ~ 1, alpha = y ~ .),
  data     = df,
  families = fam,
  control  = boost_control(mstop = 1000, nu = 0.001)
)
```

Runnable templates are in [`inst/examples/`](inst/examples).

## Constructors

| Function | Parameters | Use |
|---|---|---|
| `PowerLogit()` | mu, sigma, lambda | full three-parameter PL model |
| `GJS()` | mu, sigma | the special case lambda = 1 (generalized Johnson S_B) |
| `APL()` | mu, sigma, lambda, alpha | adjusted PL, mass at `c = 0` or `c = 1` |

All three take `family`, one of `"NO"`, `"TF"`, `"LO"`, `"PE"`, `"Hyp"`, `"SN"`,
`"SLASH"`, which selects the density generator function; `zeta`, the extra shape
parameter used by TF, PE, Hyp, SN and SLASH; and `stabilization`, one of `"none"`,
`"MAD"`, `"L2"`, for the gradient scaling.

The distribution functions `dAPL()`, `pAPL()`, `qAPL()` and `rAPL()` are exported too.

## Notes on the implementation

**Links.** The median uses a logit link, the dispersion and skewness a log link, and the
boundary probability a logit link.

**lambda is not regressed on covariates.** It is a scalar shared by all observations, so
it must enter the formula as `lambda = y ~ 1`. The families raise an error otherwise.

**Offsets.** The offset for the median is the logit of the sample median of the response.
For the dispersion it is obtained by a univariate optimization of the empirical risk over
`sigma` in `[0.1, 10]`, with the median fixed at the sample median and `lambda` at 1.
For the skewness it is obtained the same way, over `lambda` in `[0.1, 2]`, with the
median fixed at the sample median and the dispersion at a diffuse starting value.

**Internal calls.** The families call `PLreg:::dPL` and
`gamboostLSS:::stabilize_ngradient`, both internal to their packages. This is why PLboost
is distributed on GitHub rather than CRAN.

## Reproducing the master thesis

The [`reproduce_master_thesis/`](reproduce_master_thesis) directory holds the scripts behind every table and figure
of the master thesis, and is not part of the installed package.

| Directory | What it produces |
|---|---|
| `ch3_validation/` | boosting against `PLreg`/`IPLreg` maximum likelihood over 100 Monte Carlo replications per scenario (Sections 3.1 and 3.2): `validation_multirep.R` writes `validation_results.csv` and `validation_replications.csv`, `make_tables.py` turns the first into Tables 3.1 (PL) and 3.2 (APL), and `appendixA_boxplots.R` draws Figures A.1 to A.5 from the second. Seeded, so a rerun gives the published numbers |
| `ch3_simulation/` | the high-dimensional Monte Carlo study, `n = 800` with 1000 candidate covariates per submodel: `simulation_PL.R` (run once per generator, set `fam`) and `simulation_APL.R` store the results; the scripts in `figures/` draw every figure of Chapter 3 from them |
| `ch4_application/` | two applications, each under 50-fold cross-validation (5-fold repeated 10 times): `credit/` (loss given default, adjusted PL) and `inep/` (Conceito Enade course quality, plain PL); `essentials/selection_table.R` rebuilds the selection tables. See [`DATA.md`](reproduce_master_thesis/ch4_application/DATA.md) |

Scripts write their output to the directory given by the `PLBOOST_OUT` environment
variable, defaulting to the working directory:

```bash
PLBOOST_OUT=~/plboost-output Rscript reproduce_master_thesis/ch3_validation/validation_multirep.R
PLBOOST_OUT=~/plboost-output python3 reproduce_master_thesis/ch3_validation/make_tables.py
PLBOOST_OUT=~/plboost-output Rscript reproduce_master_thesis/ch3_validation/appendixA_boxplots.R
```

The Chapter 3 simulations are not seeded: a rerun gives statistically equivalent
results, not the same digits. The figure scripts reproduce the published figures exactly
from the stored `.RData` results.

The credit data is not redistributed here; see
[`reproduce_master_thesis/ch4_application/DATA.md`](reproduce_master_thesis/ch4_application/DATA.md).

## Residual diagnostics (experimental)

[`inst/diagnostics/`](inst/diagnostics) holds quantile residuals with simulated envelopes
for PL and APL fits, by boosting and by maximum likelihood (`residuals_boosting.R`,
function `envelope_PLboost`), and a worm plot of the quantile residuals
(`wp_PLboost.R`, function `plot_worm`). It is not part of the master thesis and is not
exported; `example_residuals.R` shows how to use it:

```r
source(system.file("diagnostics", "residuals_boosting.R", package = "PLboost"))
source(system.file("diagnostics", "wp_PLboost.R", package = "PLboost"))
```

## License

MIT. See [LICENSE](LICENSE).

## References

Queiroz, F. F. de and Ferrari, S. L. P. (2024). Power logit regression for modeling
bounded data. *Statistical Modelling* 24(1), 3–27.

Queiroz, F. F. de and Ferrari, S. L. P. (2024b). Modeling tropical tuna shifts: an
inflated power logit regression approach. *Biometrical Journal* 66(4).

Mayr, A., Fenske, N., Hofner, B., Kneib, T. and Schmid, M. (2012). Generalized additive
models for location, scale and shape for high dimensional data. *JRSS-C* 61(3), 403–427.

Hofner, B., Mayr, A. and Schmid, M. (2016). gamboostLSS: an R package for model building
and variable selection in the GAMLSS framework. *JSS* 74(1), 1–31.
