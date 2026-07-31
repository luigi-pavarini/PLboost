# PLboost

Gradient boosting families for the **power logit (PL)** class of distributions and its
inflated extension (**IPL**), for regression with a response bounded on the unit interval.

The package supplies `Families` objects compatible with
[**gamboostLSS**](https://cran.r-project.org/package=gamboostLSS), so the median,
dispersion, skewness and (in the inflated case) boundary-probability submodels are
estimated jointly by component-wise gradient boosting. Variable selection and
regularization come out of the fitting procedure itself, which makes the approach usable
when the number of candidate covariates exceeds the number of observations.

It accompanies the MSc dissertation *Boosting-based Inference and Variable Selection for
Modeling Bounded Continuous Data* (IME-USP). The power logit class is due to
Queiroz and Ferrari (2024), and the inflated extension to Queiroz and Ferrari (2024b).

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
## four-parameter inflated power logit, mass at one
fam <- IPL(family = "TF", zeta = 3, c = 1, stabilization = "MAD")

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
| `IPL()` | mu, sigma, lambda, alpha | inflated PL, mass at `c = 0` or `c = 1` |

All three take `family`, one of `"NO"`, `"TF"`, `"LO"`, `"PE"`, `"Hyp"`, `"SN"`,
`"SLASH"`, which selects the density generator function; `zeta`, the extra shape
parameter used by TF, PE, Hyp, SN and SLASH; and `stabilization`, one of `"none"`,
`"MAD"`, `"L2"`, for the gradient scaling.

The distribution functions `dIPL()`, `pIPL()`, `qIPL()` and `rIPL()` are exported too.

## Notes on the implementation

**Links.** The median uses a logit link, the dispersion and skewness a log link, and the
inflation parameter a logit link.

**lambda is not regressed on covariates.** It is a scalar shared by all observations, so
it must enter the formula as `lambda = y ~ 1`. The families raise an error otherwise.

**Offsets.** The offset for the median is the logit of the sample median of the response.
For the dispersion it is obtained by a univariate optimization of the empirical risk over
`sigma` in `[0.1, 10]`, with the median fixed at the sample median and `lambda` at 1.
For the skewness the offset is **fixed at `log(0.1)`**: `PLLambda()` calls `optimize()`
over `lambda` in `[0.1, 2]` but then overwrites the result, so the starting value is
always 0.1 regardless of the data. This is the behaviour that produced every result
reported in the dissertation, and it is kept here unchanged for reproducibility. Be aware
of it when interpreting the estimated skewness, which tends to drift upward from that
starting point.

**Internal calls.** The families call `PLreg:::dPL` and
`gamboostLSS:::stabilize_ngradient`, both internal to their packages. This is why PLboost
is distributed on GitHub rather than CRAN.

## Reproducing the dissertation

The [`reproduce/`](reproduce) directory holds the scripts behind the results, and is not
part of the installed package.

| Directory | What it produces |
|---|---|
| `ch2_validation/` | boosting against `PLreg`/`IPLreg` maximum likelihood, on simulated data where both are feasible (Tables 2.3 to 2.5) |
| `ch3_simulation/` | the high-dimensional Monte Carlo study, `n = 800` with 1000 candidate covariates per submodel |
| `ch4_application/` | two applications, each under 50-fold cross-validation (5-fold repeated 10 times): `credit/` (loss given default, inflated PL) and `inep/` (Conceito Enade course quality, plain PL) |

Scripts write their output to the directory given by the `PLBOOST_OUT` environment
variable, defaulting to the working directory:

```bash
PLBOOST_OUT=~/plboost-output Rscript reproduce/ch2_validation/validation_PL.R
```

The credit data is not redistributed here; see
[`reproduce/ch4_application/DATA.md`](reproduce/ch4_application/DATA.md).

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
