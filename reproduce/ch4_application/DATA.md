# Data for the credit application

The application models **loss given default** (LGD) on loans from the Lending Club
peer-to-peer platform. The raw file is about 1.6 GB and is redistributed under Kaggle's
terms, so it is **not** included in this repository.

## Getting the data

Download `accepted_2007_to_2018Q4.csv` from

<https://www.kaggle.com/datasets/wordsforthewise/lending-club>

and place it in this directory, or point the scripts at it:

```bash
mkdir -p reproduce/ch4_application/data
# move the downloaded csv to reproduce/ch4_application/data/
```

A Kaggle account is needed. With the Kaggle CLI:

```bash
kaggle datasets download -d wordsforthewise/lending-club -f accepted_2007_to_2018Q4.csv \
  -p reproduce/ch4_application/data --unzip
```

## What the scripts do with it

The sample analysed in the dissertation is built by filtering to loans that ended in
default and to the state of Vermont, which gives `n = 320` loans and, after one-hot
encoding of the borrower and loan characteristics, `p = 348` candidate covariates. The
response is the unrecovered fraction of the loan, restricted to `(0, 1]`; roughly 40% of
the sample sits exactly at 1, which is the point mass the inflated model is there to
handle.

Predictive performance is assessed over 100 independent random 80/20 train-test splits.
A 5-fold cross-validation was also run during the work and gave the same conclusions;
the split-based version is the one reported.

## Order of execution

1. `credit_run_2.1.r` — fits the power logit generators over the replications
2. `credit_run_2.2.r` — fits the one-inflated beta baseline (BEOI)
3. `credit_combine.r` — assembles the family comparison table
4. `credit_img_ipl_TF.r`, `credit_img_beoi.r` — coefficient and selection plots

Set `PLBOOST_OUT` to choose where results and figures are written.
