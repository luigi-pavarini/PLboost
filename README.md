# Dissertation code (IME-USP)

This repository contains the R code used in my dissertation work at **IME-USP**.
It is organized as a set of scripts for **testing**, **simulation**, and **applications**, plus a small folder of shared functions.

## Author and advisor

- **Author**: Luigi Pavarini de Lima  
  MSc student, Institute of Mathematics and Statistics and Computer Science, University of São Paulo (IME-USP)
- **Advisor**: Francisco Felipe de Queiroz (PhD)

São Paulo, 2026

## Folder map

- `01_Testing/`: scripts used to test IPL/PL setups
- `02_Simulation/`: scripts used to run simulations and produce plots/results
- `03_Applications/`: application scripts (e.g. credit and itagrade)
- `essentials/`: shared R functions used across folders
- `examples/`: small runnable examples (see below)

## How to run

### Requirements

- R (and optionally RStudio)

### Packages

Install the packages used by the scripts you plan to run. Some scripts use (at least):
`data.table`, `mltools`, `gamboostLSS`, `mboost`, `Formula`, `faux`, `PLreg`.

```r
install.packages(c("data.table", "mltools", "gamboostLSS", "mboost", "Formula", "faux", "PLreg"))
```

### Running a script
#### Example: PLboost (boosted PL regression)

There are very small templates under `examples/` (fill placeholders and run). For PLboost:

```bash
Rscript examples/plboost_template.R
```

#### Example: IPLboost (boosted inflated PL regression)

For one-inflated / zero-inflated responses:

```bash
Rscript examples/iplboost_template.R
```

## Notes about paths and data

- Some scripts currently define an **absolute** `base_path` (local path on my machine). If you are running this on a different computer, adjust `base_path` (or switch to relative paths).
-- The essentilas files like `families_IPL.R` and `families_PL.R` must have the source settled!!!

- Large datasets and generated outputs (e.g. `.RData`, large `.csv`, images, archives) are intentionally **ignored** via `.gitignore` to keep the repository lightweight.
